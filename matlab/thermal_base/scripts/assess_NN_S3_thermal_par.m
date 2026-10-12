% NN-assisted assessment for the thermal case, STRATEGY 3: the network
% input is augmented with the LP relaxation of the same instance, so each
% day needs TWO GAMS solves — first the relaxation (whose time is overhead
% of the strategy), then the fixed MIP.
%
% The phases are timed separately so that the results can be reported both
% with and without the relaxation overhead, as the article does. This
% matters much more here than in the PHES case: the relaxation costs about
% a second while the assisted MIP takes a few seconds.
%
% Requirements:
%   - Parallel Computing Toolbox and pinWorker.m in the path.
%   - N copies of the thermal model folder, each containing BOTH
%     MIP_GCI_SIMPL_sinPHSP.gms and RMIP_GCI_SIMPL_sinPHSP.gms, plus
%     data_gms/results_rmip.gms.
%   - NN_S3.mat with netTrained, ndMu, ndSigma (saved by define_NN_S3.m)
%   - Fixing_bands_S3.mat, calibrated with THAT network.

% ---- Machine / experiment parameters ----------------------------------
numWorkers  = 4;
gamsExe     = uc_config().gamsExe;
workerRoot  = uc_config().workerRoot;
workerPrefix= 'worker_th';
runMIP      = 'MIP_GCI_SIMPL_sinPHSP.gms';
runRMIP     = 'RMIP_GCI_SIMPL_sinPHSP.gms';
logName     = 'MIP_GCI_SIMPL_sinPHSP.log';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

% ---- Result containers -------------------------------------------------
Z = zeros(numDays,1);
U = zeros(numDays,16,24);  DELTA = zeros(numDays,16,24);  POWER = zeros(numDays,16,24);
WC = zeros(numDays,24);    UR = zeros(numDays,16,24);
COMP_TIME = zeros(numDays,1);  OPT_GAP = zeros(numDays,1);
T = zeros(numDays,1);  TREL = zeros(numDays,1);  PT = zeros(numDays,1);
FIXED_VAR = zeros(numDays,1);  OK = false(numDays,1);

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
    wf = fullfile(workerRoot, sprintf('%s%d', workerPrefix, wid));
    try
        th_writeNetDemand(fullfile(wf,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));

        % --- Phase 1: LP relaxation (overhead of Strategy 3) ------------
        relFile = fullfile(wf, 'result_rmip.txt');
        if isfile(relFile), delete(relFile); end
        t0 = tic;
        system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runRMIP));
        TREL(d) = toc(t0);
        Rday = readRelaxation(relFile);              % [16 x 24]
        Rday = sortGroupsDay(Rday, idGroups);        % canonical (descending)
        UR(d,:,:) = Rday;

        % --- Phase 2: prediction and fixing masks -----------------------
        t0 = tic;
        x = [ (ND(d,:)'-ndMu)/ndSigma , Rday' ];     % [24 x 17]
        y = extractdata(predict(netTrained, dlarray(single(x), 'TC')));
        if size(y,1) ~= numGen, y = y'; end
        u_values = double(y);                        % [16 x 24]

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
        FIXED_VAR(d) = th_writeFixedVars(fullfile(wf,'data_gms'), fixON, fixOFF, d);
        PT(d) = toc(t0);

        % --- Phase 3: fixed MIP -----------------------------------------
        resFile = fullfile(wf, 'result.txt');
        if isfile(resFile), delete(resFile); end
        t0 = tic;
        system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runMIP));
        T(d) = toc(t0);

        res = th_readResults(resFile);
        Z(d) = res.z;
        U(d,:,:) = res.u;  DELTA(d,:,:) = res.delta;  POWER(d,:,:) = res.P;
        WC(d,:) = res.wc;
        [OPT_GAP(d), COMP_TIME(d)] = th_readLog(fullfile(wf, logName));
        OK(d) = true;
        fprintf('Day %4d (worker %d): Z = %.2f, gap = %.5f | rel %.1f s + MIP %.1f s\n', ...
            d, wid, Z(d), OPT_GAP(d), TREL(d), T(d));
    catch err
        fprintf(2, 'Day %4d FAILED (worker %d): %s\n', d, wid, err.message);
    end
end

% ---- Summary and save --------------------------------------------------
nFail = nnz(~OK);
fprintf('\nCompleted: %d/%d days. Mean gap %.5f.\n', numDays-nFail, numDays, mean(OPT_GAP(OK)));
fprintf('  MIP solo:                media %6.2f s  mediana %6.2f s\n', ...
    mean(T(OK)), median(T(OK)));
fprintf('  relajacion (sobrecarga): media %6.2f s  = %.1f %% del total\n', ...
    mean(TREL(OK)), 100*mean(TREL(OK))/mean(TREL(OK)+PT(OK)+T(OK)));
fprintf('  total extremo a extremo: media %6.2f s  mediana %6.2f s\n', ...
    mean(TREL(OK)+PT(OK)+T(OK)), median(TREL(OK)+PT(OK)+T(OK)));
fprintf('  fijadas %.1f/dia (%.1f %%) | dias al limite: %d\n', ...
    mean(FIXED_VAR(OK)), 100*mean(FIXED_VAR(OK))/384, nnz(T(OK) >= 0.98*tilim));
if nFail > 0
    fprintf(2, 'FAILED days: %s\n', mat2str(find(~OK)'));
end
save(outFile, 'ND','U','DELTA','POWER','WC','UR','Z','COMP_TIME','OPT_GAP', ...
     'T','TREL','PT','FIXED_VAR','OK');
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

function u = readRelaxation(resFile)
% Reads result_rmip.txt: skips the objective value and returns the relaxed
% commitment matrix [16 x 24].
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result_rmip.txt no encontrado');
cleaner = onCleanup(@() fclose(fid));
fgetl(fid);
fscanf(fid, '%f', 1);
while true
    line = fgetl(fid);
    assert(ischar(line), 'Seccion de u no encontrada en result_rmip.txt');
    if contains(line, 'Solution: variable u(i=row, t=col)'), break; end
end
u = fscanf(fid, '%f', [24, 16])';
end
