function [nThermal, nHydro] = phTR_writeFixedVars(dataGmsFolder, fixON, fixOFF, idxTurb, idxPump, d)
%PHTR_WRITEFIXEDVARS  Writes fixed_variables.gms for the hydrothermal model
%   with intertemporal constraints.
%
%   Masks are [28 x 24] logicals in the canonical (stacked) space. Thermal
%   layers map one-to-one onto u(i,t). Turbine layers map to physical units
%   upwards (layer j -> unit j) and pump layers DOWNWARDS (layer j ->
%   unit 7-j): together with the nT + nP <= 6 cap enforced during
%   calibration, that guarantees no reversible unit is fixed to turbine and
%   pump in the same hour, which LogicPSHP forbids.
%
%   The masks must already have gone through phTR_filterFixings: this
%   function only writes what it is given.
fid = fopen(fullfile(dataGmsFolder, 'fixed_variables.gms'), 'w');
assert(fid ~= -1, 'No se pudo escribir fixed_variables.gms en %s', dataGmsFolder);
fprintf(fid, '* Variables for d = %d\n', d);

nThermal = 0;
for i = 1:16
    for t = 1:24
        if fixON(i,t)
            fprintf(fid, ' u.fx(''%d'',''%d'') = 1;\n', i, t); nThermal = nThermal + 1;
        elseif fixOFF(i,t)
            fprintf(fid, ' u.fx(''%d'',''%d'') = 0;\n', i, t); nThermal = nThermal + 1;
        end
    end
end

nHydro = 0;
for j = 1:numel(idxTurb)
    i = idxTurb(j);  c = j;
    for t = 1:24
        if fixON(i,t)
            fprintf(fid, ' ug.fx(''%d'',''%d'') = 1;\n', c, t); nHydro = nHydro + 1;
        elseif fixOFF(i,t)
            fprintf(fid, ' ug.fx(''%d'',''%d'') = 0;\n', c, t); nHydro = nHydro + 1;
        end
    end
end
for j = 1:numel(idxPump)
    i = idxPump(j);  c = 7 - j;
    for t = 1:24
        if fixON(i,t)
            fprintf(fid, ' up.fx(''%d'',''%d'') = 1;\n', c, t); nHydro = nHydro + 1;
        elseif fixOFF(i,t)
            fprintf(fid, ' up.fx(''%d'',''%d'') = 0;\n', c, t); nHydro = nHydro + 1;
        end
    end
end

fclose(fid);
end
