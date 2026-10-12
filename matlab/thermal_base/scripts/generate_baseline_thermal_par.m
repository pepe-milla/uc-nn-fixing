% Baseline campaign for the thermal case (no PSHP), run in parallel.
% Solves the unassisted UC for every day of the chosen period, with the
% fixed_variables.gms of each worker emptied beforehand so that no residual
% fixing can contaminate the reference results.
%
% Requirements:
%   - Parallel Computing Toolbox and pinWorker.m in the path.
%   - N copies of the thermal model folder in <workerRoot>, named
%     <workerPrefix>1 .. <workerPrefix>N.

% ---- Machine / experiment parameters ----------------------------------
numWorkers  = 4;
gamsExe     = uc_config().gamsExe;
workerRoot  = uc_config().workerRoot;
workerPrefix= 'worker_th';                  % worker_th1 .. worker_th4
runFile     = 'MIP_GCI_SIMPL_sinPHSP.gms';
logName     = 'MIP_GCI_SIMPL_sinPHSP.log';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

% ---- Result containers -------------------------------------------------
Z = zeros(numDays,1);
U = zeros(numDays,16,24);  DELTA = zeros(numDays,16,24);  POWER = zeros(numDays,16,24);
WC = zeros(numDays,24);
COMP_TIME = zeros(numDays,1);  OPT_GAP = zeros(numDays,1);  T = zeros(numDays,1);
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
    wf = fullfile(workerRoot, sprintf('%s%d', workerPrefix, wid));
    try
        th_writeNetDemand(fullfile(wf,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));
        resFile = fullfile(wf, 'result.txt');
        if isfile(resFile), delete(resFile); end

        tstart = tic;
        system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runFile));
        T(d) = toc(tstart);

        res = th_readResults(resFile);
        Z(d) = res.z;
        U(d,:,:) = res.u;  DELTA(d,:,:) = res.delta;  POWER(d,:,:) = res.P;
        WC(d,:) = res.wc;
        [OPT_GAP(d), COMP_TIME(d)] = th_readLog(fullfile(wf, logName));
        OK(d) = true;
        fprintf('Day %4d done (worker %d): Z = %.2f, gap = %.5f, %.1f s\n', ...
            d, wid, Z(d), OPT_GAP(d), T(d));
    catch err
        fprintf(2, 'Day %4d FAILED (worker %d): %s\n', d, wid, err.message);
    end
end

% ---- Summary and save --------------------------------------------------
nFail = nnz(~OK);
fprintf('\nCompleted: %d/%d days. Mean gap %.5f, mean time %.1f s.\n', ...
    numDays-nFail, numDays, mean(OPT_GAP(OK)), mean(T(OK)));
fprintf('Days at the time limit: %d (%.1f %%)\n', ...
    nnz(T(OK) >= 0.98*tilim), 100*mean(T(OK) >= 0.98*tilim));
if nFail > 0
    fprintf(2, 'FAILED days: %s\n', mat2str(find(~OK)'));
end
save(outFile, 'ND','U','DELTA','POWER','WC','Z','COMP_TIME','OPT_GAP','T','OK');
fprintf('Saved results to %s\n', outFile);
