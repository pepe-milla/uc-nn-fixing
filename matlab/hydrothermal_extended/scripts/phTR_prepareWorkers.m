function phTR_prepareWorkers(workerRoot, workerPrefix, numWorkers, modelDir, cplexOpt, fixedVarsText)
%PHTR_PREPAREWORKERS  Creates and configures the per-worker model folders.
%
%   Each GAMS process needs its own copy of the model directory, because
%   they all write net_demand.gms, fixed_variables.gms and result.txt with
%   the same names. Folders are created from modelDir on first use, so there
%   is no manual set-up step to forget.
%
%   Emptying (or writing) fixed_variables.gms on every launch is not
%   optional. A file left there by a previous NN campaign silently fixed
%   units 3, 4, 15 and 16 for all 24 hours of all 731 days of a label
%   campaign in this project, and nothing in the output said so.
for k = 1:numWorkers
    wf = fullfile(workerRoot, sprintf('%s%d', workerPrefix, k));
    if ~isfolder(wf)
        fprintf('Creando la carpeta del worker %d desde %s ...\n', k, modelDir);
        assert(isfolder(modelDir), 'No existe el directorio del modelo: %s', modelDir);
        copyfile(modelDir, wf);
    end

    fid = fopen(fullfile(wf, 'cplex.op9'), 'w');
    assert(fid ~= -1, 'No se pudo escribir cplex.op9 en %s', wf);
    fprintf(fid, '%s', cplexOpt);
    fclose(fid);

    fid = fopen(fullfile(wf, 'data_gms', 'fixed_variables.gms'), 'w');
    assert(fid ~= -1, 'No se pudo escribir fixed_variables.gms en %s', wf);
    fprintf(fid, '%s', fixedVarsText);
    fclose(fid);

    if isfile(fullfile(wf, 'result.txt')), delete(fullfile(wf, 'result.txt')); end

    % The model must be the intertemporal one. Checking a parameter that
    % only exists there costs nothing and rules out a folder copied from the
    % simplified model, which would run happily and produce wrong labels.
    dt = fileread(fullfile(wf, 'data_gms', 'data_TGUs_PSHP.gms'));
    assert(contains(dt, 'RU(i)') && contains(dt, 'TU(i)'), ...
        ['La carpeta %s no contiene el modelo con restricciones intertemporales ' ...
         '(faltan RU o TU en data_TGUs_PSHP.gms).'], wf);
end
end
