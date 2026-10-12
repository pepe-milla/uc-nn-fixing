
% Assembles the parts produced by the instances of generate_training_data.m
% into a single MIP_Results_2019-2021.mat, and reports what the campaign
% actually produced.
%
% Every part is sized for the whole period and only fills the days it owns
% (MINE), so merging is a matter of taking each day from the part that owns
% it. The script refuses to merge parts that do not agree on the period, and
% warns about days no part solved.

scriptDir = fileparts(mfilename('fullpath'));
resDir    = fullfile(scriptDir, '..', 'results');
parts     = dir(fullfile(resDir, 'MIP_Results_2019-2021_part*.mat'));
assert(~isempty(parts), 'No hay ficheros part*.mat en %s', resDir);
outFile   = fullfile(resDir, 'MIP_Results_2019-2021.mat');

fprintf('Partes encontradas: %d\n', numel(parts));
S = load(fullfile(resDir, parts(1).name));
numDays = size(S.ND,1);
ND = S.ND;

fields3 = {'U','Ug','Up','DELTA','DELTAg','DELTAp','POWER','POWERg','POWERp'};
fields2 = {'V','WC'};
fields1 = {'Z','COMP_TIME','OPT_GAP','T','VIOL'};

acc = struct();
for f = [fields3 fields2 fields1]
    acc.(f{1}) = nan(size(S.(f{1})));
end
OK = false(numDays,1);
owner = zeros(numDays,1);

for k = 1:numel(parts)
    P = load(fullfile(resDir, parts(k).name));
    assert(isequal(size(P.ND), size(ND)) && max(abs(P.ND - ND), [], 'all') < 1e-6, ...
        'La parte %s no corresponde al mismo periodo.', parts(k).name);
    mine = P.MINE(:) & P.OK(:);
    dup  = mine & OK;
    if any(dup)
        warning('Dias solapados entre partes (se conserva el primero): %s', mat2str(find(dup)'));
        mine = mine & ~OK;
    end
    for f = fields3
        acc.(f{1})(mine,:,:) = P.(f{1})(mine,:,:);
    end
    for f = fields2
        acc.(f{1})(mine,:) = P.(f{1})(mine,:);
    end
    for f = fields1
        if isfield(P, f{1}), acc.(f{1})(mine) = P.(f{1})(mine); end
    end
    OK(mine) = true;
    owner(mine) = k;
    fprintf('  %-40s aporta %3d dias\n', parts(k).name, nnz(mine));
end

U = acc.U; Ug = acc.Ug; Up = acc.Up;
DELTA = acc.DELTA; DELTAg = acc.DELTAg; DELTAp = acc.DELTAp;
POWER = acc.POWER; POWERg = acc.POWERg; POWERp = acc.POWERp;
V = acc.V; WC = acc.WC;
Z = acc.Z; COMP_TIME = acc.COMP_TIME; OPT_GAP = acc.OPT_GAP; T = acc.T; VIOL = acc.VIOL;

save(outFile, 'ND','U','Ug','Up','Z','DELTA','DELTAg','DELTAp', ...
     'POWER','POWERg','POWERp','V','WC','COMP_TIME','OPT_GAP','T','OK','VIOL');

fprintf('\n--- resumen de la campana ---\n');
fprintf('dias resueltos: %d de %d\n', nnz(OK), numDays);
if any(~OK)
    fprintf(2, 'dias SIN resolver: %s\n', mat2str(find(~OK)'));
    fprintf(2, 'relanzalos antes de entrenar: una etiqueta que falta es un dia menos.\n');
end
fprintf('coste medio %.0f EUR | gap medio %.4f | max %.4f\n', ...
    mean(Z(OK)), mean(OPT_GAP(OK)), max(OPT_GAP(OK)));
fprintf('dias sin converger al 1.5 %%: %d (%.1f %%)\n', ...
    nnz(OPT_GAP(OK) > 0.015 + 1e-9), 100*mean(OPT_GAP(OK) > 0.015 + 1e-9));
fprintf('tiempo: media %.1f s | mediana %.1f s | max %.1f s | total %.1f h\n', ...
    mean(COMP_TIME(OK)), median(COMP_TIME(OK)), max(COMP_TIME(OK)), sum(COMP_TIME(OK))/3600);
fprintf('vertido eolico: %d dias, %.1f MWh\n', ...
    nnz(sum(WC(OK,:),2) > 1e-6), sum(WC(OK,:), 'all'));
fprintf('violaciones intertemporales: %d (debe ser 0)\n', sum(VIOL(OK), 'omitnan'));
fprintf('Guardado en %s\n', outFile);
