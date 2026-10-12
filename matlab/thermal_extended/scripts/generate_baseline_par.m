
% Parallel baseline (pure MIP, no fixed variables) for the THERMAL model with
% intertemporal constraints, over the 730 days of 2021-2023.
%
% This period was already solved serially on the laptop. Re-running it on the
% i9 with 4 workers costs about an hour and is worth it: it puts the thermal
% and hydrothermal intertemporal cases, and the simplified cases they are
% compared against, on the same machine with the same worker configuration.
% Absolute times only mean something within one such setup, so this is what
% makes the comparison between formulations valid at all.
%
% Requirements: Parallel Computing Toolbox and pinWorker.m on the path.
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
fprintf('Baseline termico intertemporal: %d dias, %d pendientes, %d workers.\n', ...
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
            resFile = fullfile(wf, 'result.txt');
            if isfile(resFile), delete(resFile); end

            tstart = tic;
            system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runFile));
            twall = toc(tstart);

            [gap, csec, infeas] = tr_readLog(fullfile(wf, logName));
            if infeas || ~isfile(resFile)
                fprintf(2, 'Dia %4d (worker %d): infactible o sin result.txt\n', d, wid);
                continue
            end
            res = tr_readResults(resFile);
            [okDay, nViol, why] = tr_checkSolution(res, ND(d,:));
            if ~okDay
                fprintf(2, 'Dia %4d (worker %d) DESCARTADO: %s\n', d, wid, why);
                continue
            end

            saveDay(fullfile(cacheDir, sprintf('day_%04d.mat', d)), res, gap, csec, twall, nViol);
            fprintf('Dia %4d (worker %d): Z = %9.0f | gap %.5f | %6.1f s\n', d, wid, res.z, gap, csec);
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
OK = false(numDays,1);      VIOL = nan(numDays,1);

for d = 1:numDays
    f = fullfile(cacheDir, sprintf('day_%04d.mat', d));
    if ~isfile(f), continue; end
    S = load(f);
    Z(d) = S.z;  U(d,:,:) = S.u;  DELTA(d,:,:) = S.delta;  POWER(d,:,:) = S.P;
    WC(d,:) = S.wc;
    OPT_GAP(d) = S.gap;  COMP_TIME(d) = S.csec;  T(d) = S.twall;  VIOL(d) = S.nViol;
    OK(d) = true;
end

save(outFile, 'ND','U','Z','DELTA','POWER','WC','COMP_TIME','OPT_GAP','T','OK','VIOL');

fprintf('\n--- resumen ---\n');
fprintf('dias resueltos: %d de %d\n', nnz(OK), numDays);
if any(~OK)
    fprintf(2, 'dias SIN resolver: %s\n', mat2str(find(~OK)'));
    fprintf(2, ['los dias 660 y 661 son infactibles por construccion: su demanda en la ' ...
                'hora 1 supera lo que el parque puede alcanzar desde el estado inicial.\n']);
end
fprintf('coste medio %.0f EUR | gap medio %.5f | max %.5f\n', ...
    mean(Z(OK)), mean(OPT_GAP(OK)), max(OPT_GAP(OK)));
fprintf('dias sin converger al %.1f %%: %d\n', 100*epgap, nnz(OPT_GAP(OK) > epgap + 1e-9));
fprintf('tiempo por dia: media %.2f s | mediana %.2f s | max %.2f s | total %.2f h\n', ...
    mean(COMP_TIME(OK)), median(COMP_TIME(OK)), max(COMP_TIME(OK)), sum(COMP_TIME(OK))/3600);
fprintf('vertido eolico: %d dias, %.1f MWh\n', nnz(sum(WC(OK,:),2) > 1e-6), sum(WC(OK,:),'all'));
fprintf('violaciones intertemporales: %d (debe ser 0)\n', sum(VIOL(OK),'omitnan'));
fprintf('Guardado en %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function saveDay(f, res, gap, csec, twall, nViol)
z = res.z; u = res.u; delta = res.delta; P = res.P; wc = res.wc;
save(f, 'z','u','delta','P','wc','gap','csec','twall','nViol');
end

function pinWorker()
% Pins this worker to logical processors 4(k-1) and 4(k-1)+2, the mapping
% of the previous campaigns. Local function so it travels with the parfor
% body. Never throws: it is called outside the loop's try block, and a
% failure here would abort the whole parfor. The timings of this campaign
% are reported, so check that the four "anclado" lines appear at the start.
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
