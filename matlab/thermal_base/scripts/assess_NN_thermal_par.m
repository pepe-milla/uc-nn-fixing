% NN-assisted assessment for the thermal case (no PSHP), STRATEGY 2:
% binary variables are fixed wherever the calibrated confidence bands
% certify the prediction. Parallel counterpart of assess_NN.m.
%
% To run STRATEGY 1 instead, set strategy1 = true below: the trustful-unit
% mask is then used and every variable of those units is fixed to its
% predicted value, as in the article.
%
% Requirements:
%   - Parallel Computing Toolbox and pinWorker.m in the path.
%   - N copies of the thermal model folder (<workerPrefix>1 .. N).
%   - NN.mat and Fixing_bands.mat (or Strategy1_units.mat) from the SAME run.

% ---- Machine / experiment parameters ----------------------------------
% strategy1 and outFile may be set by the caller (see run_campaigns_i9.m);
% the values below are only defaults for a standalone run.
if ~exist('strategy1','var'), strategy1 = false; end   % false = S2, true = S1

numWorkers  = 4;
gamsExe     = uc_config().gamsExe;
workerRoot  = uc_config().workerRoot;
workerPrefix= 'worker_th';
runFile     = 'MIP_GCI_SIMPL_sinPHSP.gms';
logName     = 'MIP_GCI_SIMPL_sinPHSP.log';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

% ---- Predictions and fixing masks (serial: inference is fast) ----------
FIXON  = false(numGen, 24, numDays);
FIXOFF = false(numGen, 24, numDays);
PT = zeros(numDays,1);
for d = 1:numDays
    t0 = tic;
    y = extractdata(predict(netTrained, dlarray(single(ND(d,:)'), 'TC')));
    if size(y,1) ~= numGen, y = y'; end
    u_values = double(y);                            % [16 x 24]

    if strategy1
        pred   = canonicalizeDay(u_values > 0.5, idGroups);
        fixON  =  pred & trustful(:);
        fixOFF = ~pred & trustful(:);
    else
        fixON  = u_values > hiBand(:);
        fixOFF = u_values < loBand(:);
        % Stacked-group consistency: ON at layer j forces ON at shallower
        % layers, OFF at layer j forces OFF at deeper ones.
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
    end
    PT(d) = toc(t0);
    FIXON(:,:,d)  = fixON;
    FIXOFF(:,:,d) = fixOFF;
end
fprintf('Fijadas: %.1f de 384 binarias por dia (%.1f %%)\n', ...
    mean(squeeze(sum(FIXON,[1 2])+sum(FIXOFF,[1 2]))), ...
    100*mean(squeeze(sum(FIXON,[1 2])+sum(FIXOFF,[1 2])))/384);

% ---- Result containers -------------------------------------------------
Z = zeros(numDays,1);
U = zeros(numDays,16,24);  DELTA = zeros(numDays,16,24);  POWER = zeros(numDays,16,24);
WC = zeros(numDays,24);
COMP_TIME = zeros(numDays,1);  OPT_GAP = zeros(numDays,1);  T = zeros(numDays,1);
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
        FIXED_VAR(d) = th_writeFixedVars(fullfile(wf,'data_gms'), ...
            FIXON(:,:,d), FIXOFF(:,:,d), d);

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
fprintf('Days at the time limit: %d (%.1f %%) | fixed %.1f/day (%.1f %%)\n', ...
    nnz(T(OK) >= 0.98*tilim), 100*mean(T(OK) >= 0.98*tilim), ...
    mean(FIXED_VAR(OK)), 100*mean(FIXED_VAR(OK))/384);
if nFail > 0
    fprintf(2, 'FAILED days (check for infeasibility): %s\n', mat2str(find(~OK)'));
end
save(outFile, 'ND','U','DELTA','POWER','WC','Z','COMP_TIME','OPT_GAP', ...
     'T','PT','FIXED_VAR','OK','strategy1');
fprintf('Saved results to %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function X = canonicalizeDay(X, groups)
% Stacked canonical encoding for a single day, X = [numGen x 24].
for k = 1:numel(groups)
    idx = groups{k};
    cnt = sum(X(idx,:), 1);
    for j = 1:numel(idx)
        X(idx(j),:) = cnt >= j;
    end
end
end
