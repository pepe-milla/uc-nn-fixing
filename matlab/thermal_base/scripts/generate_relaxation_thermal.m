% Generates the LP-relaxation solutions used as extra input by Strategy 3
% (thermal case). Serial version for the laptop: each relaxation is an LP of
% a couple of seconds, so the whole training period runs in minutes.
%
% Use this only for the TRAINING period. For the test period the relaxation
% must be solved inside the assessment campaign, on the same machine and
% configuration as the MIP runs, since its time is overhead of Strategy 3.
%
% Output: ../data/Training_Data_S3.mat = the labels of the training period
% plus UR (relaxed commitments) and ZR (relaxed objective).

% ---- Parameters --------------------------------------------------------
gamsExe    = uc_config().gamsExe;
modelDir   = fullfile(uc_config().gamsRoot, 'thermal_base');
runFile    = 'RMIP_GCI_SIMPL_sinPHSP.gms';
% Hourly demand and wind of the period, from the data of this repository
[DAYDEM, DAYWIND, ND] = uc_load_profiles('2019-2021');
numDays = size(ND, 1);
NDrel = DAYDEM - DAYWIND;

% Labels of the same period, merged into the output file
S = load(labelFile, 'ND', 'U');
assert(size(S.ND,1) == numDays, ...
    ['El periodo (%d dias) no coincide con las etiquetas de %s (%d dias). ' ...
     'Revisa init_time / final_time.'], numDays, labelFile, size(S.ND,1));
assert(max(abs(S.ND(:) - NDrel(:))) < 1e-6, ...
    'La demanda neta calculada no coincide con la de las etiquetas.');

% ---- Solve the relaxation for every day --------------------------------
UR  = zeros(numDays,16,24);
ZR  = zeros(numDays,1);
TR  = zeros(numDays,1);
OKR = false(numDays,1);

for d = 1:numDays
    try
        th_writeNetDemand(fullfile(modelDir,'data_gms'), DAYDEM(d,:), DAYWIND(d,:));
        resFile = fullfile(modelDir, 'result_rmip.txt');
        if isfile(resFile), delete(resFile); end

        cmd = sprintf('cd /d "%s" && "%s" %s lo=2', modelDir, gamsExe, runFile);
        t0 = tic;
        system(cmd);
        TR(d) = toc(t0);

        res = readRelaxation(resFile);
        ZR(d)     = res.z;
        UR(d,:,:) = res.u;
        OKR(d) = true;
        if mod(d,50)==0 || d==1
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
frac = mean(UR(:) > 1e-6 & UR(:) < 1-1e-6);
fprintf('Fraccion de variables relajadas con valor fraccionario: %.1f %%\n', 100*frac);

% Sanity check: if the relaxed values came back as exact integers only, the
% results file is rounding them and the Strategy 3 input would be useless.
assert(frac > 0, ...
    ['Ninguna variable relajada es fraccionaria: revisa que results_rmip.gms ' ...
     'escribe con decimales y que el modelo se resolvio como RMIP.']);

% How often the integral relaxed values agree with the MIP solution: this is
% the reference against which the network has to prove its worth.
integral = (UR < 1e-6) | (UR > 1-1e-6);
match = round(UR) == S.U;
fprintf('Relajada entera en %.1f %% de las variables; de esas coinciden con el MIP %.2f %%\n', ...
    100*mean(integral(:)), 100*mean(match(integral)));

if any(~OKR)
    fprintf(2, 'Dias fallidos: %s\n', mat2str(find(~OKR)'));
end

ND = S.ND;  U = S.U;
save(outFile, 'ND','U','UR','ZR','TR','OKR');
fprintf('Guardado en %s\n', outFile);

% ---- Local functions ---------------------------------------------------
function res = readRelaxation(resFile)
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result_rmip.txt no encontrado');
cleaner = onCleanup(@() fclose(fid));
fgetl(fid);                                   % cabecera
res.z = fscanf(fid, '%f', 1);
while true
    line = fgetl(fid);
    assert(ischar(line), 'Seccion de u no encontrada');
    if contains(line, 'Solution: variable u(i=row, t=col)'), break; end
end
res.u = fscanf(fid, '%f', [24, 16])';
end
