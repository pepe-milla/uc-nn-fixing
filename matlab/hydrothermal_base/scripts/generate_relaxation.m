% Generates the LP-relaxation solutions used as extra input by Strategy 3
% (PHES case). Serial version, intended for the laptop: each relaxation is
% an LP that takes a few seconds, so the whole training period runs in
% about an hour and no parallelism is needed.
%
% IMPORTANT: use this only for the TRAINING period. For the test period the
% relaxation must be solved inside the assessment campaign, on the same
% machine and configuration as the MIP runs, because its computation time
% counts as overhead of Strategy 3.
%
% Output: ../data/Training_Data_S3.mat = the labels of Training_Data.mat
% plus UR, URg, URp (relaxed commitments) and ZR (relaxed objective).

% ---- Parameters --------------------------------------------------------
gamsExe    = uc_config().gamsExe;
modelDir   = fullfile(uc_config().gamsRoot, 'hydrothermal_base');
runFile    = 'RMIP_GCI_SIMPL.gms';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2019-2021');
numDays = size(ND, 1);
NDrel = DAYDEM - DAYWIND;

% Labels of the same period, to be merged into the output file
S = load(labelFile, 'ND', 'U', 'Ug', 'Up');
assert(size(S.ND,1) == numDays, ...
    ['El periodo (%d dias) no coincide con las etiquetas de %s (%d dias). ' ...
     'Revisa init_time / final_time.'], numDays, labelFile, size(S.ND,1));
assert(max(abs(S.ND(:) - NDrel(:))) < 1e-6, ...
    'La demanda neta calculada no coincide con la de las etiquetas.');

% ---- Solve the relaxation for every day --------------------------------
UR  = zeros(numDays,16,24);
URg = zeros(numDays, 6,24);
URp = zeros(numDays, 6,24);
ZR  = zeros(numDays,1);
TR  = zeros(numDays,1);
OKR = false(numDays,1);

for d = 1:numDays
    try
        writeNetDemand(fullfile(modelDir,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));
        resFile = fullfile(modelDir, 'result_rmip.txt');
        if isfile(resFile), delete(resFile); end

        cmd = sprintf('cd /d "%s" && "%s" %s lo=2', modelDir, gamsExe, runFile);
        t0 = tic;
        system(cmd);
        TR(d) = toc(t0);

        res = readRelaxation(resFile);
        ZR(d)      = res.z;
        UR(d,:,:)  = res.u;
        URg(d,:,:) = res.ug;
        URp(d,:,:) = res.up;
        OKR(d) = true;
        if mod(d,25)==0 || d==1
            fprintf('Dia %4d/%d: ZR = %.2f, %.1f s\n', d, numDays, ZR(d), TR(d));
        end
    catch err
        fprintf(2, 'Dia %4d FALLO: %s\n', d, err.message);
    end
end
cd(scriptDir);

% ---- Report and save ---------------------------------------------------
fprintf('\nCompletados %d/%d dias. Tiempo medio %.1f s (total %.1f min).\n', ...
    nnz(OKR), numDays, mean(TR(OKR)), sum(TR)/60);
frac = mean([UR(:); URg(:); URp(:)] > 1e-6 & [UR(:); URg(:); URp(:)] < 1-1e-6);
fprintf('Fraccion de variables relajadas con valor fraccionario: %.1f %%\n', 100*frac);
if any(~OKR)
    fprintf(2, 'Dias fallidos: %s\n', mat2str(find(~OKR)'));
end

ND = S.ND;  U = S.U;  Ug = S.Ug;  Up = S.Up;
save(outFile, 'ND','U','Ug','Up','UR','URg','URp','ZR','TR','OKR');
fprintf('Guardado en %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function writeNetDemand(dataGmsFolder, DayDemand, DayWind)
fid = fopen(fullfile(dataGmsFolder, 'net_demand.gms'), 'w');
fprintf(fid, '* Data reading from Matlab with generate_relaxation.m\n');
fprintf(fid, 'Scalar\nwcpen wind curtailment penalty (EUR:MWh) /1000000/;\n');
fprintf(fid, 'Parameter\nDemand(t) Hourly REAL power demands\n/\n');
for h = 1:24, fprintf(fid, '%d %s\n', h, num2str(DayDemand(h))); end
fprintf(fid, '/\nWind(t) Scenarios of wind power production\n/\n');
for h = 1:24, fprintf(fid, '%d  %s\n', h, num2str(DayWind(h))); end
fprintf(fid, '/\n;\nParameter\nD(t) Net power demand (MW);\nD(t)=Demand(t)-Wind(t)\n;\n');
fclose(fid);
end

function res = readRelaxation(resFile)
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result_rmip.txt no encontrado');
cleaner = onCleanup(@() fclose(fid));
fgetl(fid);                                   % cabecera
res.z  = fscanf(fid, '%f', 1);
res.u  = readMat(fid, 'Solution: variable u(i=row, t=col)',  16);
res.ug = readMat(fid, 'Solution: variable ug(c=row, t=col)',  6);
res.up = readMat(fid, 'Solution: variable up(c=row, t=col)',  6);
end

function M = readMat(fid, header, nRows)
while true
    line = fgetl(fid);
    assert(ischar(line), 'Seccion no encontrada: %s', header);
    if contains(line, header), break; end
end
M = fscanf(fid, '%f', [24, nRows])';
end
