
% Parallel generation of UC labels for the hydrothermal model WITH
% intertemporal constraints (Modelo GC simpl_T_Ramp), over the 731 days of
% 2019-2021 - the same period used by every other case in this project.
%
% Splits the days among N GAMS workers running concurrently, each in its own
% copy of the model folder, and assembles the same matrices as the PHES
% pipeline.
%
% Requirements:
%   - Parallel Computing Toolbox and pinWorker.m on the path.
%   - Nothing else: the worker folders are created from the model directory
%     on first use.
%
% Only two settings are machine-specific (gamsExe and workerRoot); the rest
% of the paths are derived from the location of this file, so the script
% works wherever the project tree lives.
%
% Each day is written to a small cache file as it completes. If the run dies
% halfway - 14 hours is a long time - relaunching skips everything already
% cached instead of starting over.

% ---- Machine-specific ---------------------------------------------------
gamsExe      = uc_config().gamsExe;
workerRoot   = uc_config().workerRoot;
numWorkers   = 4;
threadsPerWorker = 2;
% -------------------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
caseDir   = fileparts(scriptDir);                       % the case folder
projDir   = fileparts(fileparts(caseDir));              % ...\Trabajo
modelDir  = fullfile(uc_config().gamsRoot, 'hydrothermal_extended');
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2019-2021');
numDays = size(ND, 1);

% ---- Days still to solve -----------------------------------------------
todo = true(numDays,1);
for d = 1:numDays
    if isfile(fullfile(cacheDir, sprintf('day_%04d.mat', d))), todo(d) = false; end
end
fprintf('Etiquetas PHES intertemporal: %d dias, %d pendientes, %d workers.\n', ...
    numDays, nnz(todo), numWorkers);
if nnz(todo) == 0
    fprintf('Todo cacheado. Ejecuta la fase de ensamblado directamente.\n');
else
    fprintf('Estimacion: ~%.1f h de reloj.\n\n', nnz(todo)*277/numWorkers/0.97/3600);
end

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
            resFile = fullfile(wf, 'result.txt');
            if isfile(resFile), delete(resFile); end

            tstart = tic;
            system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runFile));
            twall = toc(tstart);

            [gap, csec, infeas] = phTR_readLog(fullfile(wf, logName));
            if infeas || ~isfile(resFile)
                fprintf(2, 'Dia %4d (worker %d): infactible o sin result.txt\n', d, wid);
                continue
            end
            res = phTR_readResults(resFile);
            [okDay, nViol, why] = phTR_checkSolution(res, ND(d,:));
            if ~okDay
                fprintf(2, 'Dia %4d (worker %d) DESCARTADO: %s\n', d, wid, why);
                continue
            end

            saveDay(fullfile(cacheDir, sprintf('day_%04d.mat', d)), res, gap, csec, twall, nViol);
            fprintf('Dia %4d (worker %d): Z = %9.0f | gap %.4f | %6.1f s\n', d, wid, res.z, gap, csec);
        catch err
            fprintf(2, 'Dia %4d (worker %d) FALLO: %s\n', d, wid, err.message);
        end
    end
    fprintf('\nFase de calculo terminada en %.2f h.\n', toc(tCampaign)/3600);
end

% ---- Assemble the cache into the campaign file -------------------------
Z = nan(numDays,1);
U = nan(numDays,16,24);   Ug = nan(numDays,6,24);   Up = nan(numDays,6,24);
DELTA = nan(numDays,16,24); DELTAg = nan(numDays,6,24); DELTAp = nan(numDays,6,24);
POWER = nan(numDays,16,24); POWERg = nan(numDays,6,24); POWERp = nan(numDays,6,24);
V = nan(numDays,24);      WC = nan(numDays,24);
COMP_TIME = nan(numDays,1); OPT_GAP = nan(numDays,1); T = nan(numDays,1);
OK = false(numDays,1);    VIOL = nan(numDays,1);

for d = 1:numDays
    f = fullfile(cacheDir, sprintf('day_%04d.mat', d));
    if ~isfile(f), continue; end
    S = load(f);
    Z(d) = S.z;
    U(d,:,:) = S.u;      Ug(d,:,:) = S.ug;      Up(d,:,:) = S.up;
    DELTA(d,:,:) = S.delta; DELTAg(d,:,:) = S.deltag; DELTAp(d,:,:) = S.deltap;
    POWER(d,:,:) = S.P;  POWERg(d,:,:) = S.Pg;  POWERp(d,:,:) = S.Pp;
    V(d,:) = S.v;        WC(d,:) = S.wc;
    OPT_GAP(d) = S.gap;  COMP_TIME(d) = S.csec; T(d) = S.twall; VIOL(d) = S.nViol;
    OK(d) = true;
end

save(outFile, 'ND','U','Ug','Up','Z','DELTA','DELTAg','DELTAp', ...
     'POWER','POWERg','POWERp','V','WC','COMP_TIME','OPT_GAP','T','OK','VIOL');

fprintf('\n--- resumen ---\n');
fprintf('dias resueltos: %d de %d\n', nnz(OK), numDays);
if any(~OK)
    fprintf(2, 'dias SIN resolver: %s\n', mat2str(find(~OK)'));
    fprintf(2, 'relanza este script: solo repetira esos.\n');
end
fprintf('coste medio %.0f EUR | gap medio %.4f | max %.4f\n', ...
    mean(Z(OK)), mean(OPT_GAP(OK)), max(OPT_GAP(OK)));
fprintf('dias sin converger al %.1f %%: %d (%.1f %%)\n', 100*epgap, ...
    nnz(OPT_GAP(OK) > epgap + 1e-9), 100*mean(OPT_GAP(OK) > epgap + 1e-9));
fprintf('tiempo por dia: media %.1f s | mediana %.1f s | max %.1f s | total %.1f h\n', ...
    mean(COMP_TIME(OK)), median(COMP_TIME(OK)), max(COMP_TIME(OK)), sum(COMP_TIME(OK))/3600);
fprintf('vertido eolico: %d dias, %.1f MWh\n', nnz(sum(WC(OK,:),2) > 1e-6), sum(WC(OK,:),'all'));
fprintf('violaciones intertemporales: %d (debe ser 0)\n', sum(VIOL(OK),'omitnan'));
fprintf('Guardado en %s\n', outFile);
fprintf(['\nEl cache de dias sueltos queda en %s. Borralo cuando hayas comprobado\n' ...
         'el fichero de campana, no antes.\n'], cacheDir);

% ---- Local functions ---------------------------------------------------
function saveDay(f, res, gap, csec, twall, nViol)
% Written from inside the parfor, so it has to be a function: saving a whole
% struct keeps the day self-contained and lets a broken run be resumed.
z = res.z; u = res.u; ug = res.ug; up = res.up;delta = res.delta; deltag = res.deltag; deltap = res.deltap;P = res.P; Pg = res.Pg; Pp = res.Pp; v = res.v; wc = res.wc;save(f, 'z','u','ug','up','delta','deltag','deltap','P','Pg','Pp','v','wc', ...
     'gap','csec','twall','nViol');
end

function pinWorker()
% Pins this parallel worker to two logical processors, so that each GAMS
% process gets cores of its own: worker k -> logical processors 4(k-1) and
% 4(k-1)+2, the mapping of the previous campaigns ("Worker 2 anclado a
% logicos 4 y 6"). Child processes launched with system() - cmd, GAMS,
% CPLEX - inherit the affinity. Runs once per worker process.
%
% It is a local function, as in the earlier campaigns, because local
% functions of a script travel with the parfor body; a pinWorker.m file
% would have to be on every worker's path.
%
% Never throws. It is called outside the try block of the loop, so an error
% here would abort the whole parfor. If pinning fails the worker runs
% unpinned, which only adds noise to the timings - and for a label
% campaign the timings are not even reported.
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
