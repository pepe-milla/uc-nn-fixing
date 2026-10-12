% Parallel NN-assisted assessment, STRATEGY 3 (PHES case): the network
% input is augmented with the solution of the LP relaxation of the same
% instance, so each day requires TWO GAMS solves — first the relaxation
% (whose time is overhead of the strategy), then the fixed MIP.
%
% The three phases are timed separately so that the results can be reported
% both with and without the relaxation overhead, as the article does.
%
% Requirements:
%   - Parallel Computing Toolbox, and pinWorker.m in the path.
%   - N copies of the GAMS model folder, each containing BOTH
%     MIP_GCI_SIMPL.gms and RMIP_GCI_SIMPL.gms (plus data_gms/results_rmip.gms).
%   - NN_S3.mat with netTrained, ndMu, ndSigma (saved by define_NN_S3.m)
%   - Fixing_bands_S3.mat, calibrated with THAT network.

% ---- Machine / experiment parameters ----------------------------------
numWorkers = 4;
gamsExe    = uc_config().gamsExe;
workerRoot = uc_config().workerRoot;
runMIP     = 'MIP_GCI_SIMPL.gms';
runRMIP    = 'RMIP_GCI_SIMPL.gms';
logName    = 'MIP_GCI_SIMPL.log';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

% ---- Result containers -------------------------------------------------
Z = zeros(numDays,1);
U = zeros(numDays,16,24);   Ug = zeros(numDays,6,24);   Up = zeros(numDays,6,24);
DELTA = zeros(numDays,16,24); DELTAg = zeros(numDays,6,24); DELTAp = zeros(numDays,6,24);
POWER = zeros(numDays,16,24); POWERg = zeros(numDays,6,24); POWERp = zeros(numDays,6,24);
V = zeros(numDays,24);      WC = zeros(numDays,24);
UR = zeros(numDays,28,24);                       % relaxation actually used
COMP_TIME = zeros(numDays,1); OPT_GAP = zeros(numDays,1);
T = zeros(numDays,1);       TREL = zeros(numDays,1);  PT = zeros(numDays,1);
FIXED_GEN = zeros(numDays,1); FIXED_HGEN = zeros(numDays,1);
OK = false(numDays,1);

% ---- Parallel runs -----------------------------------------------------
p = gcp('nocreate');
if isempty(p) || p.NumWorkers ~= numWorkers
    delete(gcp('nocreate'));
    parpool('Processes', numWorkers);
end

parfor d = 1:numDays
    pinWorker();
    tk = getCurrentTask();
    if isempty(tk), wid = 1; else, wid = tk.ID; end
    wf = fullfile(workerRoot, sprintf('worker%d', wid));
    try
        writeNetDemand(fullfile(wf,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));

        % --- Phase 1: LP relaxation (overhead of Strategy 3) ------------
        relFile = fullfile(wf, 'result_rmip.txt');
        if isfile(relFile), delete(relFile); end
        t0 = tic;
        system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runRMIP));
        TREL(d) = toc(t0);
        rel = readRelaxation(relFile);
        Rday = [rel.u; rel.ug; rel.up];              % [28 x 24]
        Rday = sortGroupsDay(Rday, idGroups);        % canonical (descending)
        UR(d,:,:) = Rday;

        % --- Phase 2: prediction and fixing masks -----------------------
        t0 = tic;
        x = [ (ND(d,:)'-ndMu)/ndSigma , Rday' ];     % [24 x 29]
        y = extractdata(predict(netTrained, dlarray(single(x), 'TC')));
        if size(y,1) ~= numGen, y = y'; end
        u_values = double(y);                        % [28 x 24]

        fixON  = u_values > hiBand(:);
        fixOFF = u_values < loBand(:);
        for k = 1:numel(idGroups)
            idx = idGroups{k};
            for j = 2:numel(idx)
                fixON(idx(1:j-1),:)  = fixON(idx(1:j-1),:)  | fixON(idx(j),:);
                fixOFF(idx(j:end),:) = fixOFF(idx(j:end),:) | fixOFF(idx(j-1),:);
            end
        end
        conflict = fixON & fixOFF;
        fixON(conflict)  = false;
        fixOFF(conflict) = false;
        for t = 1:24
            nT = sum(fixON(idxTurb,t));  nP = sum(fixON(idxPump,t));
            while nT + nP > 6
                if nP == 0 || (nT > 0 && u_values(idxTurb(nT),t) <= u_values(idxPump(nP),t))
                    fixON(idxTurb(nT),t) = false;  nT = nT - 1;
                else
                    fixON(idxPump(nP),t) = false;  nP = nP - 1;
                end
            end
        end
        [nThermal, nHydro] = writeFixedVars(fullfile(wf,'data_gms'), ...
            fixON, fixOFF, idxTurb, idxPump, d);
        PT(d) = toc(t0);
        FIXED_GEN(d)  = nThermal;
        FIXED_HGEN(d) = nHydro;

        % --- Phase 3: fixed MIP -----------------------------------------
        resFile = fullfile(wf, 'result.txt');
        if isfile(resFile), delete(resFile); end
        t0 = tic;
        system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runMIP));
        T(d) = toc(t0);

        res = readResults(resFile);
        Z(d) = res.z;
        U(d,:,:) = res.u;         Ug(d,:,:) = res.ug;       Up(d,:,:) = res.up;
        DELTA(d,:,:) = res.delta; DELTAg(d,:,:) = res.deltag; DELTAp(d,:,:) = res.deltap;
        POWER(d,:,:) = res.P;     POWERg(d,:,:) = res.Pg;   POWERp(d,:,:) = res.Pp;
        V(d,:) = res.v;           WC(d,:) = res.wc;
        [OPT_GAP(d), COMP_TIME(d)] = readLog(fullfile(wf, logName));
        OK(d) = true;
        fprintf('Day %4d (worker %d): Z = %.2f, gap = %.4f | rel %.1f s + MIP %.1f s\n', ...
            d, wid, Z(d), OPT_GAP(d), TREL(d), T(d));
    catch err
        fprintf(2, 'Day %4d FAILED (worker %d): %s\n', d, wid, err.message);
    end
end

% ---- Summary and save --------------------------------------------------
nFail = nnz(~OK);
fprintf('\nCompleted: %d/%d days.\n', numDays-nFail, numDays);
fprintf('  MIP solo:                    media %.1f s  mediana %.1f s\n', ...
    mean(T(OK)), median(T(OK)));
fprintf('  relajacion (sobrecarga):     media %.1f s  (%.1f %% del total)\n', ...
    mean(TREL(OK)), 100*mean(TREL(OK))/mean(TREL(OK)+PT(OK)+T(OK)));
fprintf('  total extremo a extremo:     media %.1f s  mediana %.1f s\n', ...
    mean(TREL(OK)+PT(OK)+T(OK)), median(TREL(OK)+PT(OK)+T(OK)));
fprintf('  gap medio %.4f | fijadas %.0f/dia (%.1f %%)\n', mean(OPT_GAP(OK)), ...
    mean(FIXED_GEN(OK)+FIXED_HGEN(OK)), 100*mean(FIXED_GEN(OK)+FIXED_HGEN(OK))/672);
if nFail > 0
    fprintf(2, 'FAILED days: %s\n', mat2str(find(~OK)'));
end
save(outFile, 'ND','U','Ug','Up','UR','Z','DELTA','DELTAg','DELTAp', ...
     'POWER','POWERg','POWERp','V','WC','COMP_TIME','OPT_GAP', ...
     'T','TREL','PT','FIXED_GEN','FIXED_HGEN','OK');
fprintf('Saved results to %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function X = sortGroupsDay(X, groups)
% Descending sort within each group of identical units, per hour.
% X = [numGen x 24]. Continuous counterpart of the stacked encoding.
for k = 1:numel(groups)
    idx = groups{k};
    X(idx,:) = sort(X(idx,:), 1, 'descend');
end
end

function writeNetDemand(dataGmsFolder, DayDemand, DayWind)
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
fprintf(fid, '* Data reading from Matlab with assess_NN_S3_par.m\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /1000000/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
fclose(fid);
end

function [nThermal, nHydro] = writeFixedVars(dataGmsFolder, fixON, fixOFF, idxTurb, idxPump, d)
fid = fopen(fullfile(dataGmsFolder, 'fixed_variables.gms'), 'w');
fprintf(fid, '* Variables for d = %d (Strategy 3)\n', d);
nThermal = 0;
for i = 1:16
    for t = 1:24
        if fixON(i,t)
            fprintf(fid, ' u.fx(''%d'',''%d'') = 1;\n', i, t); nThermal = nThermal+1;
        elseif fixOFF(i,t)
            fprintf(fid, ' u.fx(''%d'',''%d'') = 0;\n', i, t); nThermal = nThermal+1;
        end
    end
end
nHydro = 0;
for j = 1:numel(idxTurb)
    i = idxTurb(j);  c = j;
    for t = 1:24
        if fixON(i,t)
            fprintf(fid, ' ug.fx(''%d'',''%d'') = 1;\n', c, t); nHydro = nHydro+1;
        elseif fixOFF(i,t)
            fprintf(fid, ' ug.fx(''%d'',''%d'') = 0;\n', c, t); nHydro = nHydro+1;
        end
    end
end
for j = 1:numel(idxPump)
    i = idxPump(j);  c = 7 - j;
    for t = 1:24
        if fixON(i,t)
            fprintf(fid, ' up.fx(''%d'',''%d'') = 1;\n', c, t); nHydro = nHydro+1;
        elseif fixOFF(i,t)
            fprintf(fid, ' up.fx(''%d'',''%d'') = 0;\n', c, t); nHydro = nHydro+1;
        end
    end
end
fclose(fid);
end

function res = readRelaxation(resFile)
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result_rmip.txt no encontrado');
cleaner = onCleanup(@() fclose(fid));
fgetl(fid);
res.z  = fscanf(fid, '%f', 1);
res.u  = readMat(fid, 'Solution: variable u(i=row, t=col)',  16);
res.ug = readMat(fid, 'Solution: variable ug(c=row, t=col)',  6);
res.up = readMat(fid, 'Solution: variable up(c=row, t=col)',  6);
end

function res = readResults(resFile)
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result.txt not found');
cleaner = onCleanup(@() fclose(fid));
fgetl(fid);
res.z = fscanf(fid, '%f', 1);
res.u      = readMat(fid, 'Solution: variable u(i=row, t=col)',      16);
res.ug     = readMat(fid, 'Solution: variable ug(c=row, t=col)',      6);
res.up     = readMat(fid, 'Solution: variable up(c=row, t=col)',      6);
res.delta  = readMat(fid, 'Solution: variable delta(i=row, t=col)',  16);
res.deltag = readMat(fid, 'Solution: variable deltag(c=row, t=col)',  6);
res.deltap = readMat(fid, 'Solution: variable deltap(c=row, t=col)',  6);
res.P      = readMat(fid, 'Solution: variable P(i=row, t=col)',      16);
res.Pg     = readMat(fid, 'Solution: variable Pg(c=row, t=col)',      6);
res.Pp     = readMat(fid, 'Solution: variable Pp(c=row, t=col)',      6);
res.v      = readMat(fid, 'Solution: variable v(t=col)',              1)';
res.wc     = readMat(fid, 'Solution: variable wc(t=col)',             1)';
end

function M = readMat(fid, header, nRows)
while true
    line = fgetl(fid);
    assert(ischar(line), 'Section not found: %s', header);
    if contains(line, header), break; end
end
M = fscanf(fid, '%f', [24, nRows])';
end

function [gap, csec] = readLog(logFile)
gap = NaN; csec = NaN;
txt = fileread(logFile);
tok = regexp(txt, 'elapsed\s+([\d:]+\.\d+)', 'tokens');
if ~isempty(tok)
    parts = sscanf(tok{end}{1}, '%d:%d:%f');
    csec = parts(1)*3600 + parts(2)*60 + parts(3);
end
tok = regexp(txt, 'Relative gap:\s+([\d\.]+)', 'tokens');
if ~isempty(tok)
    gap = str2double(tok{end}{1});
end
end
