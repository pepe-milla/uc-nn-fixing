% Parallel NN-assisted assessment (PHES case): solves every day of the
% test period with the calibrated fixing bands applied, using N concurrent
% GAMS workers. Parallel counterpart of assess_NN.m.
%
% Requirements:
%   - Parallel Computing Toolbox.
%   - N copies of the GAMS model folder: <workerRoot>\worker1..workerN.
%     Their cplex.op9 is overwritten with the settings below; their
%     fixed_variables.gms is rewritten per day by each worker.
%   - NN.mat (netTrained) and Fixing_bands.mat from the SAME training run.

% ---- Machine / experiment parameters ----------------------------------
numWorkers = 4;
gamsExe    = uc_config().gamsExe;
workerRoot = uc_config().workerRoot;
runFile    = 'MIP_GCI_SIMPL.gms';
logName    = 'MIP_GCI_SIMPL.log';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

% ---- NN predictions and fixing masks (serial: it is fast) --------------
% The net was trained on canonically encoded (stacked) labels: within each
% group of identical units, layer j means "at least j units are ON".
numGen = numel(loBand);
FIXON  = false(numGen, 24, numDays);
FIXOFF = false(numGen, 24, numDays);
PT = zeros(numDays,1);
for d = 1:numDays
    tstart = tic;
    % R2023a-compatible inference: in older releases dlnetwork/predict
    % requires a formatted dlarray input ('TC' = 24 time steps x 1
    % channel); numeric input is only accepted in newer ones. The output
    % orientation is normalized to [numGen x 24] afterwards.
    y = extractdata(predict(netTrained, dlarray(single(ND(d,:)'), 'TC')));
    if size(y,1) ~= numGen, y = y'; end
    u_values = double(y);                        % [28 x 24] layer probs

    fixON  = u_values > hiBand(:);
    fixOFF = u_values < loBand(:);

    % Stacked-group consistency: ON at layer j forces ON at shallower
    % layers; OFF at layer j forces OFF at deeper layers.
    for k = 1:numel(idGroups)
        idx = idGroups{k};
        for j = 2:numel(idx)
            fixON(idx(1:j-1),:) = fixON(idx(1:j-1),:) | fixON(idx(j),:);
            fixOFF(idx(j:end),:) = fixOFF(idx(j:end),:) | fixOFF(idx(j-1),:);
        end
    end
    conflict = fixON & fixOFF;
    fixON(conflict)  = false;
    fixOFF(conflict) = false;

    % Reversible-plant cap: at most 6 turbine+pump layers fixed ON per
    % hour. Release the least confident deep layers when violated.
    for t = 1:24
        nT = sum(fixON(idxTurb,t));
        nP = sum(fixON(idxPump,t));
        while nT + nP > 6
            if nP == 0 || (nT > 0 && u_values(idxTurb(nT),t) <= u_values(idxPump(nP),t))
                fixON(idxTurb(nT),t) = false;  nT = nT - 1;
            else
                fixON(idxPump(nP),t) = false;  nP = nP - 1;
            end
        end
    end
    PT(d) = toc(tstart);
    FIXON(:,:,d)  = fixON;
    FIXOFF(:,:,d) = fixOFF;
end

% ---- Result containers -------------------------------------------------
Z = zeros(numDays,1);
U = zeros(numDays,16,24);   Ug = zeros(numDays,6,24);   Up = zeros(numDays,6,24);
DELTA = zeros(numDays,16,24); DELTAg = zeros(numDays,6,24); DELTAp = zeros(numDays,6,24);
POWER = zeros(numDays,16,24); POWERg = zeros(numDays,6,24); POWERp = zeros(numDays,6,24);
V = zeros(numDays,24);      WC = zeros(numDays,24);
COMP_TIME = zeros(numDays,1); OPT_GAP = zeros(numDays,1); T = zeros(numDays,1);
FIXED_GEN = zeros(numDays,1); FIXED_HGEN = zeros(numDays,1);
OK = false(numDays,1);

% ---- Parallel GAMS runs ------------------------------------------------
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
        [nThermal, nHydro] = writeFixedVars(fullfile(wf,'data_gms'), ...
            FIXON(:,:,d), FIXOFF(:,:,d), idxTurb, idxPump, d);
        FIXED_GEN(d)  = nThermal;
        FIXED_HGEN(d) = nHydro;

        resFile = fullfile(wf, 'result.txt');
        if isfile(resFile), delete(resFile); end

        cmd = sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runFile);
        tstart = tic;
        system(cmd);
        T(d) = toc(tstart);

        res = readResults(resFile);
        Z(d) = res.z;
        U(d,:,:) = res.u;         Ug(d,:,:) = res.ug;       Up(d,:,:) = res.up;
        DELTA(d,:,:) = res.delta; DELTAg(d,:,:) = res.deltag; DELTAp(d,:,:) = res.deltap;
        POWER(d,:,:) = res.P;     POWERg(d,:,:) = res.Pg;   POWERp(d,:,:) = res.Pp;
        V(d,:) = res.v;           WC(d,:) = res.wc;

        [OPT_GAP(d), COMP_TIME(d)] = readLog(fullfile(wf, logName));
        OK(d) = true;
        fprintf('Day %4d done (worker %d): Z = %.2f, gap = %.4f, %.1f s, fixed %d+%d\n', ...
            d, wid, Z(d), OPT_GAP(d), T(d), FIXED_GEN(d), FIXED_HGEN(d));
    catch err
        fprintf(2, 'Day %4d FAILED (worker %d): %s\n', d, wid, err.message);
    end
end

% ---- Summary and save --------------------------------------------------
nFail = nnz(~OK);
fprintf('\nCompleted: %d/%d days. Mean gap %.4f, mean time %.1f s.\n', ...
    numDays-nFail, numDays, mean(OPT_GAP(OK)), mean(T(OK)));
fprintf('Mean fixed: %.1f thermal + %.1f hydro of 672 binaries/day (%.1f %%).\n', ...
    mean(FIXED_GEN(OK)), mean(FIXED_HGEN(OK)), ...
    100*mean(FIXED_GEN(OK)+FIXED_HGEN(OK))/672);
if nFail > 0
    fprintf(2, 'FAILED days (check for infeasibility, rerun): %s\n', mat2str(find(~OK)'));
end
save(outFile, 'ND','U','Ug','Up','Z','DELTA','DELTAg','DELTAp', ...
     'POWER','POWERg','POWERp','V','WC','COMP_TIME','OPT_GAP','T','PT', ...
     'FIXED_GEN','FIXED_HGEN','OK');
fprintf('Saved results to %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function writeNetDemand(dataGmsFolder, DayDemand, DayWind)
% Writes net_demand.gms with the same format as write_net_demand.m
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
fprintf(fid, '* Data reading from Matlab with assess_NN_par.m\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /1000000/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
fclose(fid);
end

function [nThermal, nHydro] = writeFixedVars(dataGmsFolder, fixON, fixOFF, idxTurb, idxPump, d)
% Writes fixed_variables.gms. Thermal layers map one-to-one; turbine
% layers map to physical units upwards (layer j -> unit j) and pump layers
% DOWNWARDS (layer j -> unit 7-j): together with the nT+nP <= 6 cap this
% guarantees no unit is fixed to turbine and pump in the same hour.
fid = fopen(fullfile(dataGmsFolder, 'fixed_variables.gms'), 'w');
fprintf(fid, '* Variables for d = %d\n', d);
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

function res = readResults(resFile)
% Reads result.txt (same section order as results.gms)
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result.txt not found');
cleaner = onCleanup(@() fclose(fid));
fgetl(fid);                                   % header line
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
% Extracts the last reported elapsed time and relative gap from the log
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

function pinWorker()
persistent done
if ~isempty(done), return; end
tk = getCurrentTask();
if isempty(tk), return; end
wid  = double(tk.ID);                          % 1..4
mask = bitshift(int64(5), 4*(wid-1));          % 5 = binario 101
proc = System.Diagnostics.Process.GetCurrentProcess();
proc.ProcessorAffinity = System.IntPtr(mask);
done = true;
fprintf('Worker %d anclado a logicos %d y %d\n', wid, 4*(wid-1), 4*(wid-1)+2);
end
