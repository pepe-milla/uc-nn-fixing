% Parallel generation of UC training labels (PHES case).
% Splits the days among N GAMS workers running concurrently, each in its
% own copy of the model folder (worker1..workerN), and assembles the same
% output matrices as the serial generate_training_data.m.
%
% Requirements:
%   - Parallel Computing Toolbox (parpool/parfor).
%   - N copies of the GAMS model folder: <workerRoot>\worker1 .. workerN
%     (full copies, including data_gms). Their cplex.op9 is overwritten
%     by this script, and their fixed_variables.gms is emptied (baseline
%     guard), so stale content cannot contaminate the labels.
%
% All machine-specific settings are parameters below.
%
% This script was used twice in this case, with the parameters edited each
% time:
%   - to generate the LABELS: period 2019-2021, time limit 600 s, output
%     ../data/Training_Data.mat. The stored labels show that time limit.
%   - to run the unassisted BASELINE of the evaluation period: period
%     2021-2023, time limit 1200 s, output ../results/MIP_Baseline_2021-2023.mat.
% The values below are those of the baseline, which is the run reported in
% the paper; the assisted campaigns of this case use the same 1200 s.

% ---- Machine / experiment parameters ----------------------------------
numWorkers = 4;                    % 
gamsExe    = uc_config().gamsExe;
workerRoot = uc_config().workerRoot;
runFile    = 'MIP_GCI_SIMPL.gms';
logName    = 'MIP_GCI_SIMPL.log';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2021-2023');
numDays = size(ND, 1);

% Result containers (same layout as the serial script)
Z = zeros(numDays,1);
U = zeros(numDays,16,24);   Ug = zeros(numDays,6,24);   Up = zeros(numDays,6,24);
DELTA = zeros(numDays,16,24); DELTAg = zeros(numDays,6,24); DELTAp = zeros(numDays,6,24);
POWER = zeros(numDays,16,24); POWERg = zeros(numDays,6,24); POWERp = zeros(numDays,6,24);
V = zeros(numDays,24);      WC = zeros(numDays,24);
COMP_TIME = zeros(numDays,1); OPT_GAP = zeros(numDays,1); T = zeros(numDays,1);
OK = false(numDays,1);

% Parallel pool
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
        fprintf('Day %4d done (worker %d): Z = %.2f, gap = %.4f, %.1f s\n', ...
            d, wid, Z(d), OPT_GAP(d), T(d));
    catch err
        fprintf(2, 'Day %4d FAILED (worker %d): %s\n', d, wid, err.message);
    end
end

% Summary and save
nFail = nnz(~OK);
fprintf('\nCompleted: %d/%d days. Mean gap %.4f, mean time %.1f s.\n', ...
    numDays-nFail, numDays, mean(OPT_GAP(OK)), mean(T(OK)));
if nFail > 0
    fprintf(2, 'FAILED days (rerun them): %s\n', mat2str(find(~OK)'));
end
save(outFile, 'ND','U','Ug','Up','Z','DELTA','DELTAg','DELTAp', ...
     'POWER','POWERg','POWERp','V','WC','COMP_TIME','OPT_GAP','T','OK');
fprintf('Saved results to %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function writeNetDemand(dataGmsFolder, DayDemand, DayWind)
% Writes net_demand.gms with the same format as write_net_demand.m
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
fprintf(fid, '* Data reading from Matlab with generate_training_data_par.m\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /1000000/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
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
