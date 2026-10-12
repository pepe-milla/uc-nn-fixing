
% Generation of UC labels for the hydrothermal model WITH intertemporal
% constraints (Modelo GC simpl_T_Ramp), over the same 731 days of 2019-2021
% used by every other case in this project.
%
% This laptop has no Parallel Computing Toolbox, so there is no parfor here.
% Parallelism is obtained instead by launching several MATLAB instances at
% once, each with a different `instance` number below. Each instance owns
% its own model folder, its own output file and a disjoint, interleaved
% subset of the days; merge_training_data.m then assembles them.
%
%   Instance 1:  set instance = 1, run.
%   Instance 2:  open a second MATLAB, set instance = 2, run.
%   When both finish:  merge_training_data
%
% Two instances of two threads each beat one instance of four: the parallel
% speedup CPLEX gets from 2 to 4 threads on a branch-and-bound tree is well
% short of 2x, so splitting the days wins. Do not go beyond two on a 4-core
% machine.
%
% Parallelism is free here in a way it is not elsewhere: the solve times of
% a label campaign are never reported, only the labels matter. The timing
% campaigns must still run serially so their measurements stay clean.
%
% It is restartable. Results are saved every `chunkSize` days and, if the
% output file already exists for the same period, the days already solved
% are skipped. A run of two and a half days should not have to start over
% because the laptop rebooted.

% ---- Instance ----------------------------------------------------------
instance     = 1;    % <<-- 1 for the first MATLAB, 2 for the second
numInstances = 2;
% ------------------------------------------------------------------------

% ---- Machine / experiment parameters -----------------------------------
threadsPerRun = 2;   % with numInstances = 2 on 4 physical cores
gamsExe    = uc_config().gamsExe;
modelDir   = fullfile(uc_config().gamsRoot, 'hydrothermal_extended');
workerRoot = uc_config().workerRoot;
runFile    = 'MIP_GCI_SIMPL.gms';
logName    = 'MIP_GCI_SIMPL.log';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2019-2021');
numDays = size(ND, 1);

% Days of this instance, interleaved so both cover the whole year evenly
myDays = instance:numInstances:numDays;

% Containers, sized for the full period so the parts merge trivially
Z = nan(numDays,1);
U = nan(numDays,16,24);   Ug = nan(numDays,6,24);   Up = nan(numDays,6,24);
DELTA = nan(numDays,16,24); DELTAg = nan(numDays,6,24); DELTAp = nan(numDays,6,24);
POWER = nan(numDays,16,24); POWERg = nan(numDays,6,24); POWERp = nan(numDays,6,24);
V = nan(numDays,24);      WC = nan(numDays,24);
COMP_TIME = nan(numDays,1); OPT_GAP = nan(numDays,1); T = nan(numDays,1);
OK = false(numDays,1);    VIOL = nan(numDays,1);
MINE = false(numDays,1);  MINE(myDays) = true;

% Resume, but only if the stored file really is the same period
if isfile(outFile)
    S = load(outFile);
    if isfield(S,'OK') && isfield(S,'ND') && isequal(size(S.ND), size(ND)) ...
            && max(abs(S.ND - ND), [], 'all') < 1e-6
        Z = S.Z; U = S.U; Ug = S.Ug; Up = S.Up;
        DELTA = S.DELTA; DELTAg = S.DELTAg; DELTAp = S.DELTAp;
        POWER = S.POWER; POWERg = S.POWERg; POWERp = S.POWERp;
        V = S.V; WC = S.WC; COMP_TIME = S.COMP_TIME; OPT_GAP = S.OPT_GAP;
        T = S.T; OK = S.OK;
        if isfield(S,'VIOL'), VIOL = S.VIOL; end
        fprintf('Reanudando la instancia %d: %d de %d dias ya resueltos.\n', ...
            instance, nnz(OK), numel(myDays));
    else
        warning('%s existe pero no corresponde a este periodo. No se reanuda.', outFile);
    end
end

pending = myDays(~OK(myDays));
if isempty(pending)
    fprintf('La instancia %d no tiene nada pendiente.\n', instance);
    return
end
fprintf(['Instancia %d de %d: %d dias pendientes de %d.\n' ...
         'Estimacion inicial ~%.1f h (a 333 s por dia, %d hilos).\n\n'], ...
    instance, numInstances, numel(pending), numel(myDays), ...
    numel(pending)*333/3600, threadsPerRun);

tStart = tic; doneHere = 0;
for c = 1:numel(pending)
    d = pending(c);
    try
        writeNetDemand(fullfile(wf,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));
        if isfile(resFile), delete(resFile); end

        t0 = tic;
        [~,~] = system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, runFile));
        T(d) = toc(t0);
        [OPT_GAP(d), COMP_TIME(d)] = readLog(logFile);

        if isInfeasible(logFile) || ~isfile(resFile)
            fprintf(2, 'Dia %4d: infactible o sin result.txt\n', d);
        else
            res = readResults(resFile);
            gen = sum(res.P,1) + sum(res.Pg,1) - sum(res.Pp,1);
            if ~(res.z > 0) || max(abs(gen - (ND(d,:) + res.wc))) > 1
                fprintf(2, 'Dia %4d: incoherente (Z = %.2f, balance roto)\n', d, res.z);
            else
                VIOL(d) = countViolations(round(res.u), round(res.delta), res.P);
                if VIOL(d) > 0
                    fprintf(2, 'Dia %4d: viola %d restricciones intertemporales\n', d, VIOL(d));
                else
                    Z(d) = res.z;
                    U(d,:,:) = res.u;      Ug(d,:,:) = res.ug;     Up(d,:,:) = res.up;
                    DELTA(d,:,:) = res.delta; DELTAg(d,:,:) = res.deltag; DELTAp(d,:,:) = res.deltap;
                    POWER(d,:,:) = res.P;  POWERg(d,:,:) = res.Pg; POWERp(d,:,:) = res.Pp;
                    V(d,:) = res.v;        WC(d,:) = res.wc;
                    OK(d) = true;
                    fprintf('Dia %4d: Z = %9.0f | gap %.4f | %6.1f s\n', ...
                        d, res.z, OPT_GAP(d), COMP_TIME(d));
                end
            end
        end
    catch err
        fprintf(2, 'Dia %4d FALLO: %s\n', d, err.message);
    end

    doneHere = doneHere + 1;
    if mod(doneHere, chunkSize) == 0 || c == numel(pending)
        save(outFile, 'ND','U','Ug','Up','Z','DELTA','DELTAg','DELTAp', ...
             'POWER','POWERg','POWERp','V','WC','COMP_TIME','OPT_GAP','T','OK', ...
             'VIOL','MINE','instance','numInstances');
        el = toc(tStart);
        fprintf('--- guardado: %d/%d de esta instancia | %.1f h | ~%.1f h restantes ---\n\n', ...
            doneHere, numel(pending), el/3600, el/3600*(numel(pending)-doneHere)/doneHere);
    end
end

fprintf('\nInstancia %d completada: %d de %d dias.\n', instance, nnz(OK(myDays)), numel(myDays));
bad = myDays(~OK(myDays));
if ~isempty(bad), fprintf(2, 'Dias sin resolver: %s\n', mat2str(bad)); end
fprintf('gap medio %.4f | tiempo medio %.1f s\n', ...
    mean(OPT_GAP(OK),'omitnan'), mean(COMP_TIME(OK),'omitnan'));
fprintf('Cuando terminen las %d instancias, ejecuta merge_training_data.\n', numInstances);

% ---- Local functions ---------------------------------------------------
function writeNetDemand(dataGmsFolder, DayDemand, DayWind)
% wcpen = 0: the penalty never activates in the hydrothermal model, since
% the pumps absorb the surplus, and a value of 1e6 would turn any day where
% curtailment does become unavoidable into a multi-million-euro outlier.
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
fprintf(fid, '* Data reading from Matlab with generate_training_data.m\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /0/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
fclose(fid);
end

function tf = isInfeasible(logFile)
% Only the solver's verdict counts. The bare word "infeasible" also appears
% in the status column of any branch-and-bound node whose relaxation is
% infeasible, which is routine; matching it discarded a good day earlier in
% this project.
tf = false;
if ~isfile(logFile), return; end
txt = lower(fileread(logFile));
tf = contains(txt, 'integer infeasible') || contains(txt, 'problem is infeasible');
end

function res = readResults(resFile)
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result.txt no encontrado');
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
res.v      = readMat(fid, 'Solution: variable v(t=col)',              1);
res.wc     = readMat(fid, 'Solution: variable wc(t=col)',             1);
end

function M = readMat(fid, header, nRows)
while true
    line = fgetl(fid);
    assert(ischar(line), 'Seccion no encontrada: %s', header);
    if contains(line, header), break; end
end
M = fscanf(fid, '%f', [24, nRows])';
end

function [gap, csec] = readLog(logFile)
% GAMS writes the job time as H:MM:SS.sss and the fields must be combined:
% taking only the last number of the line drops the minutes.
gap = NaN; csec = NaN;
if ~isfile(logFile), return; end
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

function v = countViolations(u, delta, P)
% Ramps and minimum times of the thermal units. The PSH units carry none of
% these constraints by design and are not checked. A non-zero count means
% GAMS did not solve this model.
Pmin = [6.79 6.79 27.84 27.84 75.5 75.5 6.79 6.79 14.09 14.09 4.58 4.58 4.58 6.79 17.7 17.7]';
RU   = [24.26 24.26 29.70 29.70 103.05 113.50 24.26 24.26 20.51 20.51 8.51 8.51 8.51 13.23 22.22 22.22]';
TU   = [1 1 7 7 4 4 1 1 1 1 1 1 1 1 6 6];
TD   = [1 1 7 7 3 3 1 1 1 1 1 1 1 1 5 5];
P0   = [32.3 32.3 37 0 0 0 32.3 32.3 20.5 20.5 8.51 8.51 8.51 17.64 37 0]';
u0   = [1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 0]';
uprev = [u0, u(:,1:end-1)];  Pprev = [P0, P(:,1:end-1)];
w  = delta - (u - uprev);
vU = (P - Pprev) - (RU.*uprev + Pmin.*delta);
vD = (Pprev - P) - (RU.*u     + Pmin.*w);
v  = nnz(vU > 0.05) + nnz(vD > 0.05);
v  = v + nnz(u == 1 & P < Pmin - 1e-4) + nnz(u == 0 & P > 1e-4);
for i = 1:16
    for t = 1:24
        if t >= TU(i) && sum(delta(i, max(1,t-TU(i)+1):t)) > u(i,t) + 1e-6, v = v + 1; end
        if t >= TD(i) && sum(w(i,     max(1,t-TD(i)+1):t)) > 1 - u(i,t) + 1e-6, v = v + 1; end
    end
end
end
