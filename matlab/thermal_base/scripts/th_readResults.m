function res = th_readResults(resFile)
% Reads result.txt of the thermal model (no PSHP): objective value and the
% matrices u, delta and P (16 units x 24 hours) plus the wind curtailment
% vector, in the order written by data_gms/results.gms.
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result.txt no encontrado');
cleaner = onCleanup(@() fclose(fid));
fgetl(fid);                                          % header line
res.z     = fscanf(fid, '%f', 1);
res.u     = readMat(fid, 'Solution: variable u(i=row, t=col)',     16);
res.delta = readMat(fid, 'Solution: variable delta(i=row, t=col)', 16);
res.P     = readMat(fid, 'Solution: variable P(i=row, t=col)',     16);
res.wc    = readMat(fid, 'Solution: variable wc(t=col)',            1)';
end

function M = readMat(fid, header, nRows)
while true
    line = fgetl(fid);
    assert(ischar(line), 'Seccion no encontrada: %s', header);
    if contains(line, header), break; end
end
M = fscanf(fid, '%f', [24, nRows])';
end
