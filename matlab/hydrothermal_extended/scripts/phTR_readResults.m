function res = phTR_readResults(resFile)
%PHTR_READRESULTS  Reads result.txt of the hydrothermal model.
%   Sections come in the order written by data_gms/results.gms. Each block
%   is one row per unit and 24 columns, so fscanf with size [24 nRows] fills
%   columns with file rows and the transpose gives [units x hours].
fid = fopen(resFile, 'r');
assert(fid ~= -1, 'result.txt no encontrado: %s', resFile);
cleaner = onCleanup(@() fclose(fid));

fgetl(fid);                                  % header line
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
    assert(ischar(line), 'Seccion no encontrada en result.txt: %s', header);
    if contains(line, header), break; end
end
M = fscanf(fid, '%f', [24, nRows])';
end
