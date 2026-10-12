function res = tr_readResults(resFile)
%TR_READRESULTS  Reads result.txt of the thermal model with intertemporal
%   constraints: objective value, u, delta and P (16 units x 24 hours) and
%   the wind curtailment vector, in the order written by
%   data_gms/results.gms.
%
%   Each block is one row per unit and 24 columns, so fscanf with size
%   [24 nRows] fills columns with file rows and the transpose gives
%   [units x hours]. wc is returned as a 1 x 24 row.
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result.txt no encontrado: %s', resFile);
cleaner = onCleanup(@() fclose(fid));

fgetl(fid);                                          % header line
res.z     = fscanf(fid, '%f', 1);
res.u     = readMat(fid, 'Solution: variable u(i=row, t=col)',     16);
res.delta = readMat(fid, 'Solution: variable delta(i=row, t=col)', 16);
res.P     = readMat(fid, 'Solution: variable P(i=row, t=col)',     16);
res.wc    = readMat(fid, 'Solution: variable wc(t=col)',            1);
end

function M = readMat(fid, header, nRows)
while true
    line = fgetl(fid);
    assert(ischar(line), 'Seccion no encontrada en result.txt: %s', header);
    if contains(line, header), break; end
end
M = fscanf(fid, '%f', [24, nRows])';
end
