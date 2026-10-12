function run_solver_metrics(cases, k, K, onlyDays)
%RUN_SOLVER_METRICS  Re-solves baseline and Strategy 2 keeping the CPLEX log,
%   to extract solver metrics (nodes, iterations, cuts, root vs branching
%   time, incumbents) that the i9 campaigns did not record.
%
%   run_solver_metrics(cases, k, K)
%     cases  cell array with any of 'thermal_simpl', 'thermal_ext', 'hydro_simpl',
%            'hydro_ext'
%            (processed in that order), or a single char
%     k, K   this is instance k of K: it takes the days with mod(d-1,K) == k-1,
%            so K MATLAB processes can share the work without overlapping
%
%   Every day is solved twice, baseline first and Strategy 2 right after, in
%   a private copy of the model (work\<case>_<k>). The settings are those of
%   the i9 campaigns, including threads=2 and deterministic parallel mode: on
%   the days that close before the time limit CPLEX should then follow the
%   same search as on the i9, which sm_summary checks by comparing the
%   objective values. Days that stop at the time limit are machine-dependent
%   and must not be compared with the i9 campaigns.
%
%   The fixings of Strategy 2 are recomputed from the stored network and
%   bands exactly as in the campaigns (group stacking, turbine+pump cap,
%   minimum-time filter in the extended formulation).
%
%     onlyDays  optional: restrict every case to these days (for testing)
%
%   Restartable: each (case, mode, day) is cached as it completes, together
%   with a copy of its CPLEX log. Nothing in the original model folders or
%   campaign results is written to.

if ischar(cases), cases = {cases}; end
scriptDir = fileparts(mfilename('fullpath'));
caseRoot  = fileparts(scriptDir);                         % ...\6_Solver metrics
projDir   = uc_config().repoRoot;
M         = fullfile(projDir, 'matlab');
gamsExe   = uc_config().gamsExe;
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

TU = [1 1 7 7 4 4 1 1 1 1 1 1 1 1 6 6];
TD = [1 1 7 7 3 3 1 1 1 1 1 1 1 1 5 5];
u0 = [1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 0];
gT = {[1 2],[3 4],[7 8],[9 10],[11 12 13],[15 16]};

for ci = 1:numel(cases)
    cfg = caseConfig(cases{ci}, M, projDir);
    addpath(cfg.helperDir);
    G = gT;  if cfg.nG == 28, G = [gT, {17:22}, {23:28}]; end
    N  = load(cfg.netFile, 'netTrained');
    Bd = load(cfg.bandsFile, 'hiBand', 'loBand');

    wf = fullfile(caseRoot, 'work', sprintf('%s_%d', cfg.name, k));
    prepareWork(cfg.modelDir, wf, cfg.runFile, cfg.cplexOpt);
    cacheDir = fullfile(caseRoot, 'results', 'cache', cfg.name);
    rawDir   = fullfile(caseRoot, 'results', 'raw_logs', cfg.name);
    if ~isfolder(cacheDir), mkdir(cacheDir); end
    if ~isfolder(rawDir),   mkdir(rawDir);   end

    if nargin >= 4 && ~isempty(onlyDays), cfg.days = intersect(cfg.days, onlyDays); end
    myDays = cfg.days(mod(cfg.days-1, K) == k-1);
    say(sprintf('=== %s: %d dias para la instancia %d de %d ===', cfg.name, numel(myDays), k, K));
    tCase = tic;
    for d = myDays(:)'
        for mc = {'base', 'S2'}
            md = mc{1};
            f = fullfile(cacheDir, sprintf('%s_%04d.mat', md, d));
            if isfile(f), continue; end
            try
                % Fixings: none for the baseline, the campaign's for Strategy 2
                % formatted dlarray: works on the older release of the i9 too
                y = extractdata(predict(N.netTrained, dlarray(single(ND(d,:)'), 'TC')));
                if size(y,1) ~= cfg.nG, y = y'; end
                P = double(y);                                     % [nG x 24]
                if strcmp(md, 'base')
                    on = false(cfg.nG, 24);  off = false(cfg.nG, 24);
                else
                    [on, off] = applyBands(P, Bd.hiBand, Bd.loBand, G, cfg.nG);
                end
                cfg.writeND(fullfile(wf,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));
                nFix = cfg.writeFix(fullfile(wf,'data_gms'), on, off, d, TU, TD, u0);

                logFile = fullfile(wf, cfg.logName);
                if isfile(logFile), delete(logFile); end
                resFile = fullfile(wf, 'result.txt');
                if isfile(resFile), delete(resFile); end
                t0 = tic;
                [~, ~] = system(sprintf('cd /d "%s" && "%s" %s lo=2', wf, gamsExe, cfg.runFile));
                twall = toc(t0);

                assert(isfile(logFile), 'GAMS no escribio el log');
                copyfile(logFile, fullfile(rawDir, sprintf('%s_%04d.log', md, d)));
                [gap, csec, infeas] = cfg.readLog(logFile);
                sm = sm_parseCplexLog(logFile);
                z = NaN;
                if ~infeas && isfile(resFile)
                    res = cfg.readRes(resFile);  z = res.z;
                end
                save(f, 'd', 'md', 'z', 'gap', 'csec', 'twall', 'infeas', 'nFix', 'sm');
                say(sprintf('%s dia %4d %-4s: Z %10.0f | gap %.4f | %6.1f s | nodos %7g | cortes %4d | fijadas %3d', ...
                    cfg.name, d, md, z, gap, csec, sm.nodes, sm.cuts, nFix));
            catch err
                say(sprintf('%s dia %4d %s FALLO: %s', cfg.name, d, md, err.message));
            end
        end
    end
    say(sprintf('=== %s terminado en %.2f h ===', cfg.name, toc(tCase)/3600));
end
end

% ========================================================================
function cfg = caseConfig(name, M, projDir)
G = fullfile(projDir, 'Gams');
opt = @(epgap, tilim) sprintf(['epgap=%g\nitlim=100000000\nmipemphasis=2\nnodefileind=2\n' ...
    'parallelmode=1\nthreads=2\ntilim=%g\nworkmem=6000\n'], epgap, tilim);
cfg.name = name;
switch name
    case 'thermal_simpl'
        d = fullfile(M, 'thermal_base');
        cfg.modelDir = fullfile(G, 'thermal_base');
        cfg.runFile = 'MIP_GCI_SIMPL_sinPHSP.gms';  cfg.logName = 'MIP_GCI_SIMPL_sinPHSP.log';
        cfg.netFile = fullfile(d, 'trained models', 'NN.mat');
        cfg.bandsFile = fullfile(d, 'data', 'Fixing_bands_eps005.mat');
        cfg.helperDir = fullfile(d, 'scripts');
        cfg.nG = 16;  cfg.cplexOpt = opt(0.001, 120);  cfg.days = 1:730;
        cfg.writeND  = @th_writeNetDemand;
        cfg.writeFix = @(g, on, off, dd, ~, ~, ~) th_writeFixedVars(g, on, off, dd);
        cfg.readLog  = @(f) readLog3(@th_readLog, f);
        cfg.readRes  = @th_readResults;
    case 'thermal_ext'
        d = fullfile(M, 'thermal_extended');
        cfg.modelDir = fullfile(G, 'Modelo GC simpl_T_Ramp_sinPSHP');
        cfg.runFile = 'MIP_GCI_SIMPL_sinPHSP.gms';  cfg.logName = 'MIP_GCI_SIMPL_sinPHSP.log';
        cfg.netFile = fullfile(d, 'trained models', 'NN.mat');
        cfg.bandsFile = fullfile(d, 'data', 'Fixing_bands.mat');
        cfg.helperDir = fullfile(d, 'scripts');
        cfg.nG = 16;  cfg.cplexOpt = opt(0.001, 120);  cfg.days = 1:730;
        cfg.writeND  = @tr_writeNetDemand;
        cfg.writeFix = @(g, on, off, dd, TU, TD, u0) tr_writeFixedVars(g, on, off, dd, TU, TD, u0);
        cfg.readLog  = @tr_readLog;
        cfg.readRes  = @tr_readResults;
    case 'hydro_simpl'
        % The campaign's own writers are local functions of assess_NN_par.m in
        % "hydrothermal_base"; their file formats are identical to the phTR_*
        % ones (same pump mapping, same result.txt sections), which are reused
        % here with the two differences of the base formulation: curtailment
        % penalty of 1e6 and no minimum-time filter.
        % NOTE: the copy of assess_NN_par.m on the laptop still carries the
        % old settings (0.01, 600 s, 4 threads, opportunistic mode); the
        % stored results show 0.015 and 1200 s, and threads=2 in
        % deterministic mode is assumed as in the other campaigns. The
        % comparison of objective values in sm_summary tests that assumption.
        d = fullfile(M, 'hydrothermal_base');
        cfg.modelDir = fullfile(G, 'Modelo GC simplificado');
        cfg.runFile = 'MIP_GCI_SIMPL.gms';  cfg.logName = 'MIP_GCI_SIMPL.log';
        cfg.netFile = fullfile(d, 'trained models', 'NN.mat');
        cfg.bandsFile = fullfile(d, 'data', 'Fixing_bands_S2.mat');
        cfg.helperDir = fullfile(M, 'hydrothermal_extended', 'scripts');
        cfg.nG = 28;  cfg.cplexOpt = opt(0.015, 1200);
        cfg.days = 6:12:730;                                  % same sample as hydro_ext
        cfg.writeND  = @writeNetDemandBase;
        cfg.writeFix = @writeFixHydroSimpl;
        cfg.readLog  = @phTR_readLog;
        cfg.readRes  = @phTR_readResults;
    case 'hydro_ext'
        d = fullfile(M, 'hydrothermal_extended');
        cfg.modelDir = fullfile(G, 'Modelo GC simpl_T_Ramp');
        cfg.runFile = 'MIP_GCI_SIMPL.gms';  cfg.logName = 'MIP_GCI_SIMPL.log';
        cfg.netFile = fullfile(d, 'trained models', 'NN.mat');
        cfg.bandsFile = fullfile(d, 'data', 'Fixing_bands.mat');
        cfg.helperDir = fullfile(d, 'scripts');
        cfg.nG = 28;  cfg.cplexOpt = opt(0.015, 1200);
        cfg.days = 6:12:730;                                  % 61 days spread over the period
        cfg.writeND  = @phTR_writeNetDemand;
        cfg.writeFix = @writeFixHydroExt;
        cfg.readLog  = @phTR_readLog;
        cfg.readRes  = @phTR_readResults;
    otherwise
        error('Caso desconocido: %s', name);
end
end

function n = writeFixHydroSimpl(g, on, off, d, ~, ~, ~)
[nT, nH] = phTR_writeFixedVars(g, on, off, 17:22, 23:28, d);
n = nT + nH;
end

function writeNetDemandBase(dataGmsFolder, DayDemand, DayWind)
% Same file as writeNetDemand in "hydrothermal_base"\assess_NN_par.m: base
% formulation, curtailment penalty of 1e6 EUR/MWh.
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
fprintf(fid, '* Data reading from Matlab with run_solver_metrics.m\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /1000000/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
fclose(fid);
end

function n = writeFixHydroExt(g, on, off, d, TU, TD, u0)
[on, off] = phTR_filterFixings(on, off, TU, TD, u0);
[nT, nH] = phTR_writeFixedVars(g, on, off, 17:22, 23:28, d);
n = nT + nH;
end

function [gap, csec, infeas] = readLog3(fun, f)
[gap, csec] = fun(f);
txt = fileread(f);
infeas = contains(txt, 'integer infeasible') || contains(txt, 'problem is infeasible');
end

function prepareWork(modelDir, wf, runFile, cplexOpt)
% Private copy of the model: the run file and data_gms, nothing else.
if ~isfolder(fullfile(wf,'data_gms')), mkdir(fullfile(wf,'data_gms')); end
copyfile(fullfile(modelDir, runFile), fullfile(wf, runFile));
g = dir(fullfile(modelDir, 'data_gms', '*.gms'));
for i = 1:numel(g)
    copyfile(fullfile(g(i).folder, g(i).name), fullfile(wf, 'data_gms', g(i).name));
end
fid = fopen(fullfile(wf, 'cplex.op9'), 'w');  fprintf(fid, '%s', cplexOpt);  fclose(fid);
end

function [ON, OFF] = applyBands(P, hi, lo, G, nG)
% Same policy as the campaigns: bands, group stacking, contradictions left
% free, and at most six turbine+pump positions fixed ON per hour.
ON = P > hi(:);  OFF = P < lo(:);
for kk = 1:numel(G)
    idx = G{kk};
    for j = 2:numel(idx)
        ON(idx(1:j-1),:) = ON(idx(1:j-1),:) | ON(idx(j),:);
        OFF(idx(j:end),:) = OFF(idx(j:end),:) | OFF(idx(j-1),:);
    end
end
cf = ON & OFF;  ON(cf) = false;  OFF(cf) = false;
if nG == 28
    gt = 17:22;  gp = 23:28;
    for t = 1:size(P,2)
        nT = sum(ON(gt,t));  nP = sum(ON(gp,t));
        while nT + nP > 6
            if nP == 0 || (nT > 0 && P(gt(nT),t) <= P(gp(nP),t)), ON(gt(nT),t) = false; nT = nT-1;
            else, ON(gp(nP),t) = false; nP = nP-1; end
        end
    end
end
end

function f = firstExisting(c)
f = '';
for i = 1:numel(c), if isfile(c{i}), f = c{i}; return; end, end
error('No se encuentra ninguno de: %s', strjoin(c, ' | '));
end

function say(msg)
fprintf('[%s] %s\n', string(datetime('now'),'dd-MMM HH:mm:ss'), msg);
end
