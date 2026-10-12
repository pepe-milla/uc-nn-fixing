
% Strategy 2 on the THERMAL model with intertemporal constraints: solve the
% 730 days of 2021-2023 with the binaries that the calibrated confidence
% bands certify, and leave the rest to the solver. Compared against
% generate_baseline_par.m over the same days, same settings, same machine.
%
% The serial version of this campaign has already been run on the laptop
% (Assess_NN_S2_TRamp.mat). This one exists so that the thermal and
% hydrothermal intertemporal cases are measured on the same machine and the
% same worker configuration as each other and as the simplified cases.
%
% Requirements: Parallel Computing Toolbox and pinWorker.m on the path, plus
% the trained network and its bands - which must come from the SAME run of
% define_NN.m / calibrate_fixing.m, because the bands certify the error of
% that network's probabilities and mean nothing for another one.
%
% Restartable: each day is cached as it completes.

% ---- Machine-specific ---------------------------------------------------
gamsExe      = uc_config().gamsExe;
workerRoot   = uc_config().workerRoot;
numWorkers   = 4;
threadsPerWorker = 2;
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
caseDir   = fileparts(scriptDir);
projDir   = fileparts(fileparts(caseDir));
modelDir  = fullfile(uc_config().gamsRoot, 'thermal_extended');
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

todo = true(numDays,1);
for d = 1:numDays
    if isfile(fullfile(cacheDir, sprintf('day_%04d.mat', d))), todo(d) = false; end
end
fprintf('S2 termico intertemporal: %d dias, %d pendientes, %d workers.\n', ...
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
            tr_writeNetDemand(fullfile(wf,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));

            % Network prediction. Row g is the probability of canonical LAYER
            % g, not of physical unit g; fixing layer j onto physical unit
            % idGroups{k}(j) is valid because they are interchangeable.
            %
            % The input must be a FORMATTED dlarray ('TC' = 24 time steps x 1
            % channel): the i9 runs a release where dlnetwork/predict rejects
            % a plain numeric array, which newer releases accept. single() is
            % the precision the network works in, and the one used for the
            % campaign of the simplified case, so the probabilities match the
            % ones calibrate_fixing certified.
            tp = tic;
            y = extractdata(predict(netTrained, dlarray(single(ND(d,:)'), 'TC')));
            if size(y,1) ~= 16, y = y'; end
            prob = double(y);                           % [16 x 24]
            ptime = toc(tp);

            [fixON, fixOFF, nConf] = applyBands(prob, hiBand, loBand, idGroups);
            [fixON, fixOFF, nDrop] = tr_filterFixings(fixON, fixOFF, TU, TD, u0);
            nFix = tr_writeFixedVars(fullfile(wf,'data_gms'), fixON, fixOFF, d);

            resFile = fullfile(wf, 'result.txt');
            if isfile(resFile), delete(resFile); end

            tstart = tic;
            system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runFile));
            twall = toc(tstart);

            [gap, csec, infeas] = tr_readLog(fullfile(wf, logName));
            if infeas || ~isfile(resFile)
                fprintf(2, 'Dia %4d (worker %d): INFACTIBLE o sin result.txt\n', d, wid);
                continue
            end
            res = tr_readResults(resFile);
            [okDay, nViol, why] = tr_checkSolution(res, ND(d,:));
            if ~okDay
                fprintf(2, 'Dia %4d (worker %d) DESCARTADO: %s\n', d, wid, why);
                continue
            end
            % The solution must honour every variable we fixed. This catches
            % a stale fixed_variables.gms and a GAMS run in the wrong folder.
            uu = round(res.u);
            if ~(all(uu(fixON) == 1) && all(uu(fixOFF) == 0))
                fprintf(2, 'Dia %4d (worker %d) DESCARTADO: la solucion no respeta las fijaciones\n', d, wid);
                continue
            end

            saveDay(fullfile(cacheDir, sprintf('day_%04d.mat', d)), res, gap, csec, ...
                    twall, ptime, nViol, nFix, nDrop, nConf);
            fprintf('Dia %4d (worker %d): Z = %9.0f | gap %.5f | %6.1f s | fijadas %3d | liberadas %d\n', ...
                d, wid, res.z, gap, csec, nFix, nDrop);
        catch err
            fprintf(2, 'Dia %4d (worker %d) FALLO: %s\n', d, wid, err.message);
        end
    end
    fprintf('\nFase de calculo terminada en %.2f h.\n', toc(tCampaign)/3600);
end

% ---- Assemble ----------------------------------------------------------
Z = nan(numDays,1);
U = nan(numDays,16,24); DELTA = nan(numDays,16,24); POWER = nan(numDays,16,24);
WC = nan(numDays,24);
COMP_TIME = nan(numDays,1); OPT_GAP = nan(numDays,1); T = nan(numDays,1);
PT = nan(numDays,1);  OK = false(numDays,1);  VIOL = nan(numDays,1);
FIXED_VAR = nan(numDays,1); DROPPED = nan(numDays,1); NCONF = nan(numDays,1);

for d = 1:numDays
    f = fullfile(cacheDir, sprintf('day_%04d.mat', d));
    if ~isfile(f), continue; end
    S = load(f);
    Z(d) = S.z;  U(d,:,:) = S.u;  DELTA(d,:,:) = S.delta;  POWER(d,:,:) = S.P;
    WC(d,:) = S.wc;
    OPT_GAP(d) = S.gap;  COMP_TIME(d) = S.csec;  T(d) = S.twall;  PT(d) = S.ptime;
    VIOL(d) = S.nViol;   FIXED_VAR(d) = S.nFix;  DROPPED(d) = S.nDrop;
    NCONF(d) = S.nConf;  OK(d) = true;
end

save(outFile, 'ND','U','Z','DELTA','POWER','WC','COMP_TIME','OPT_GAP','T','PT', ...
     'OK','VIOL','FIXED_VAR','DROPPED','NCONF');

fprintf('\n--- resumen ---\n');
fprintf('dias resueltos: %d de %d\n', nnz(OK), numDays);
if any(~OK)
    fprintf(2, 'dias SIN resolver: %s\n', mat2str(find(~OK)'));
    fprintf(2, ['los dias 660 y 661 son infactibles por construccion; si hay mas, las ' ...
                'fijaciones estan provocando infactibilidad y hay que mirarlo.\n']);
end
fprintf('binarias fijadas: %.1f %% (%.0f de 384 por dia)\n', ...
    100*mean(FIXED_VAR(OK))/384, mean(FIXED_VAR(OK)));
fprintf('liberadas por el filtro de tiempos minimos: %d en total\n', sum(DROPPED(OK),'omitnan'));
fprintf('conflictos de grupo: %d\n', sum(NCONF(OK),'omitnan'));
fprintf('coste medio %.0f EUR | gap medio %.5f\n', mean(Z(OK)), mean(OPT_GAP(OK)));
fprintf('tiempo por dia: media %.2f s | mediana %.2f s | max %.2f s | total %.2f h\n', ...
    mean(COMP_TIME(OK)), median(COMP_TIME(OK)), max(COMP_TIME(OK)), sum(COMP_TIME(OK))/3600);
fprintf('prediccion de la red: %.1f ms por dia\n', 1000*mean(PT(OK)));
fprintf('violaciones intertemporales: %d (debe ser 0)\n', sum(VIOL(OK),'omitnan'));
fprintf('Guardado en %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function [FixON, FixOFF, nConflict] = applyBands(P, hiBand, loBand, groups)
% Same policy as calibrate_fixing.m, on a single day: [16 x 24] here instead
% of [24 x 16 x n] there. A NaN band never fixes, because every comparison
% against NaN is false. ON at layer j forces ON at shallower layers, OFF at
% layer j forces OFF at deeper layers; contradictions are left free.
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
end

function saveDay(f, res, gap, csec, twall, ptime, nViol, nFix, nDrop, nConf)
z = res.z; u = res.u; delta = res.delta; P = res.P; wc = res.wc;
save(f, 'z','u','delta','P','wc','gap','csec','twall','ptime','nViol', ...
     'nFix','nDrop','nConf');
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
