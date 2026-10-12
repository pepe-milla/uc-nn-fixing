
% Strategy 2 on the hydrothermal model WITH intertemporal constraints: solve
% the 730 days of 2021-2023 with the binaries that the calibrated confidence
% bands certify, and leave the rest to the solver. Compared against
% generate_baseline_par.m over the same days, same settings, same machine.
%
% Requirements: Parallel Computing Toolbox and pinWorker.m on the path, plus
% a trained network and its bands - and they must come from the SAME run of
% define_NN.m / calibrate_fixing.m, because the bands certify the error of
% that network's probabilities and mean nothing for another one.
%
% Restartable: each day is cached as it completes.
%
% One step here has no counterpart in the simplified case. Before the fixings
% are written they go through phTR_filterFixings, which releases any thermal
% combination that the minimum up/down times make impossible. Without it GAMS
% returns INFEASIBLE on those days instead of a costlier solution - a failure
% mode that simply did not exist when any pattern of binaries was feasible.

% ---- Machine-specific ---------------------------------------------------
gamsExe      = uc_config().gamsExe;
workerRoot   = uc_config().workerRoot;
numWorkers   = 4;
threadsPerWorker = 2;
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
caseDir   = fileparts(scriptDir);
projDir   = fileparts(fileparts(caseDir));
modelDir  = fullfile(uc_config().gamsRoot, 'hydrothermal_extended');
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

todo = true(numDays,1);
for d = 1:numDays
    if isfile(fullfile(cacheDir, sprintf('day_%04d.mat', d))), todo(d) = false; end
end
fprintf('S2 PHES intertemporal: %d dias, %d pendientes, %d workers.\n', ...
    numDays, nnz(todo), numWorkers);

% ---- Parallel runs -----------------------------------------------------
if nnz(todo) > 0
    p = gcp('nocreate');
    if isempty(p) || p.NumWorkers ~= numWorkers
        delete(gcp('nocreate'));
        parpool('Processes', numWorkers);
    end

    tCampaign = tic;
    parfor d = 1:numDays
        if ~todo(d), continue; end
        pinWorker();
        tk = getCurrentTask();
        if isempty(tk), wid = 1; else, wid = tk.ID; end
        wf = fullfile(workerRoot, sprintf('%s%d', workerPrefix, wid));
        try
            phTR_writeNetDemand(fullfile(wf,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));

            % Network prediction. Row g is the probability of canonical
            % LAYER g, not of physical unit g; phTR_writeFixedVars does the
            % mapping.
            %
            % The input must be a FORMATTED dlarray ('TC' = 24 time steps x 1
            % channel): the i9 runs a release where dlnetwork/predict rejects
            % a plain numeric array, which newer releases accept. single() is
            % the precision the network works in, and the one used for the
            % campaign of the simplified case, so the probabilities match the
            % ones calibrate_fixing certified.
            tp = tic;
            y = extractdata(predict(netTrained, dlarray(single(ND(d,:)'), 'TC')));
            if size(y,1) ~= 28, y = y'; end
            prob = double(y);                           % [28 x 24]
            ptime = toc(tp);

            [fixON, fixOFF, nConf, nCap] = applyBands(prob, hiBand, loBand, ...
                idGroups, idxTurb, idxPump);
            [fixON, fixOFF, nDrop] = phTR_filterFixings(fixON, fixOFF, TU, TD, u0);
            [nTh, nHy] = phTR_writeFixedVars(fullfile(wf,'data_gms'), ...
                fixON, fixOFF, idxTurb, idxPump, d);

            resFile = fullfile(wf, 'result.txt');
            if isfile(resFile), delete(resFile); end

            tstart = tic;
            system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runFile));
            twall = toc(tstart);

            [gap, csec, infeas] = phTR_readLog(fullfile(wf, logName));
            if infeas || ~isfile(resFile)
                fprintf(2, 'Dia %4d (worker %d): INFACTIBLE o sin result.txt\n', d, wid);
                continue
            end
            res = phTR_readResults(resFile);
            [okDay, nViol, why] = phTR_checkSolution(res, ND(d,:));
            if ~okDay
                fprintf(2, 'Dia %4d (worker %d) DESCARTADO: %s\n', d, wid, why);
                continue
            end
            % The solution must honour every variable we fixed. This catches
            % a stale fixed_variables.gms and a GAMS run in the wrong folder.
            u28 = [round(res.u); round(res.ug); flipud(round(res.up))];
            if ~(all(u28(fixON) == 1) && all(u28(fixOFF) == 0))
                fprintf(2, 'Dia %4d (worker %d) DESCARTADO: la solucion no respeta las fijaciones\n', d, wid);
                continue
            end

            saveDay(fullfile(cacheDir, sprintf('day_%04d.mat', d)), res, gap, csec, ...
                    twall, ptime, nViol, nTh, nHy, nDrop, nConf, nCap);
            fprintf('Dia %4d (worker %d): Z = %9.0f | gap %.4f | %6.1f s | fijadas %d+%d | liberadas %d\n', ...
                d, wid, res.z, gap, csec, nTh, nHy, nDrop);
        catch err
            fprintf(2, 'Dia %4d (worker %d) FALLO: %s\n', d, wid, err.message);
        end
    end
    fprintf('\nFase de calculo terminada en %.2f h.\n', toc(tCampaign)/3600);
end

% ---- Assemble ----------------------------------------------------------
Z = nan(numDays,1);
U = nan(numDays,16,24);   Ug = nan(numDays,6,24);   Up = nan(numDays,6,24);
DELTA = nan(numDays,16,24); DELTAg = nan(numDays,6,24); DELTAp = nan(numDays,6,24);
POWER = nan(numDays,16,24); POWERg = nan(numDays,6,24); POWERp = nan(numDays,6,24);
V = nan(numDays,24);      WC = nan(numDays,24);
COMP_TIME = nan(numDays,1); OPT_GAP = nan(numDays,1); T = nan(numDays,1);
PT = nan(numDays,1);      OK = false(numDays,1);     VIOL = nan(numDays,1);
FIX_TH = nan(numDays,1);  FIX_HY = nan(numDays,1);
DROPPED = nan(numDays,1); NCONF = nan(numDays,1);    NCAP = nan(numDays,1);

for d = 1:numDays
    f = fullfile(cacheDir, sprintf('day_%04d.mat', d));
    if ~isfile(f), continue; end
    S = load(f);
    Z(d) = S.z;
    U(d,:,:) = S.u;      Ug(d,:,:) = S.ug;      Up(d,:,:) = S.up;
    DELTA(d,:,:) = S.delta; DELTAg(d,:,:) = S.deltag; DELTAp(d,:,:) = S.deltap;
    POWER(d,:,:) = S.P;  POWERg(d,:,:) = S.Pg;  POWERp(d,:,:) = S.Pp;
    V(d,:) = S.v;        WC(d,:) = S.wc;
    OPT_GAP(d) = S.gap;  COMP_TIME(d) = S.csec; T(d) = S.twall;  PT(d) = S.ptime;
    VIOL(d) = S.nViol;   FIX_TH(d) = S.nTh;     FIX_HY(d) = S.nHy;
    DROPPED(d) = S.nDrop; NCONF(d) = S.nConf;   NCAP(d) = S.nCap;
    OK(d) = true;
end

save(outFile, 'ND','U','Ug','Up','Z','DELTA','DELTAg','DELTAp', ...
     'POWER','POWERg','POWERp','V','WC','COMP_TIME','OPT_GAP','T','PT','OK', ...
     'VIOL','FIX_TH','FIX_HY','DROPPED','NCONF','NCAP');

fprintf('\n--- resumen ---\n');
fprintf('dias resueltos: %d de %d\n', nnz(OK), numDays);
if any(~OK)
    fprintf(2, 'dias SIN resolver: %s\n', mat2str(find(~OK)'));
    fprintf(2, ['si son mas de los que fallan en el baseline, las fijaciones estan ' ...
                'provocando infactibilidad y hay que mirarlo antes de seguir.\n']);
end
nTot = 672;   % 28 units x 24 hours
fprintf('binarias fijadas: %.1f %% (%.0f termicas + %.0f hidro de %d por dia)\n', ...
    100*mean(FIX_TH(OK)+FIX_HY(OK))/nTot, mean(FIX_TH(OK)), mean(FIX_HY(OK)), nTot);
fprintf('liberadas por el filtro de tiempos minimos: %d en total\n', sum(DROPPED(OK),'omitnan'));
fprintf('conflictos de grupo: %d | topes turbina+bomba: %d\n', ...
    sum(NCONF(OK),'omitnan'), sum(NCAP(OK),'omitnan'));
fprintf('coste medio %.0f EUR | gap medio %.4f\n', mean(Z(OK)), mean(OPT_GAP(OK)));
fprintf('tiempo por dia: media %.1f s | mediana %.1f s | max %.1f s | total %.1f h\n', ...
    mean(COMP_TIME(OK)), median(COMP_TIME(OK)), max(COMP_TIME(OK)), sum(COMP_TIME(OK))/3600);
fprintf('prediccion de la red: %.1f ms por dia\n', 1000*mean(PT(OK)));
fprintf('violaciones intertemporales: %d (debe ser 0)\n', sum(VIOL(OK),'omitnan'));
fprintf('Guardado en %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function [FixON, FixOFF, nConflict, nCapRel] = applyBands(P, hiBand, loBand, ...
    groups, grpT, grpP)
% Same policy as calibrate_fixing.m, on a single day: [28 x 24] here instead
% of [24 x 28 x n] there. A NaN band never fixes, because every comparison
% against NaN is false. ON at layer j forces ON at shallower layers, OFF at
% layer j forces OFF at deeper layers; contradictions are left free. Then at
% most 6 turbine+pump layers may be fixed ON per hour, releasing the least
% confident deep layers when that is violated.
FixON  = P > hiBand(:);
FixOFF = P < loBand(:);
for k = 1:numel(groups)
    idx = groups{k};
    for j = 2:numel(idx)
        FixON(idx(1:j-1),:)  = FixON(idx(1:j-1),:)  | FixON(idx(j),:);
        FixOFF(idx(j:end),:) = FixOFF(idx(j:end),:) | FixOFF(idx(j-1),:);
    end
end
conflict = FixON & FixOFF;
nConflict = nnz(conflict);
FixON(conflict)  = false;
FixOFF(conflict) = false;

nCapRel = 0;
for t = 1:size(P,2)
    nT = sum(FixON(grpT,t));
    nP = sum(FixON(grpP,t));
    while nT + nP > numel(grpT)
        if nP == 0 || (nT > 0 && P(grpT(nT),t) <= P(grpP(nP),t))
            FixON(grpT(nT),t) = false;  nT = nT - 1;
        else
            FixON(grpP(nP),t) = false;  nP = nP - 1;
        end
        nCapRel = nCapRel + 1;
    end
end
end

function saveDay(f, res, gap, csec, twall, ptime, nViol, nTh, nHy, nDrop, nConf, nCap)
z = res.z; u = res.u; ug = res.ug; up = res.up;
delta = res.delta; deltag = res.deltag; deltap = res.deltap;
P = res.P; Pg = res.Pg; Pp = res.Pp; v = res.v; wc = res.wc;
save(f, 'z','u','ug','up','delta','deltag','deltap','P','Pg','Pp','v','wc', ...
     'gap','csec','twall','ptime','nViol','nTh','nHy','nDrop','nConf','nCap');
end

function pinWorker()
% Pins this worker to logical processors 4(k-1) and 4(k-1)+2, the mapping
% of the previous campaigns - and it must be the same one as the baseline,
% or the speedup mixes two configurations. Local function so it travels
% with the parfor body. Never throws: it is called outside the loop's try
% block, and a failure here would abort the whole parfor.
persistent done
if ~isempty(done), return; end
done = true;
try
    tk = getCurrentTask();
    if isempty(tk), return; end
    k = tk.ID;
    a = 4*(k-1);  b = a + 2;
    proc = System.Diagnostics.Process.GetCurrentProcess();
    proc.ProcessorAffinity = System.IntPtr(int64(2^a + 2^b));
    fprintf('Worker %d anclado a logicos %d y %d\n', k, a, b);
catch err
    fprintf(2, 'Worker: no se pudo anclar (%s). Sigue sin anclar.\n', err.message);
end
end
