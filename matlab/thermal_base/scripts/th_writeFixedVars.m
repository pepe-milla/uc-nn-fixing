function nFixed = th_writeFixedVars(dataGmsFolder, fixON, fixOFF, d)
% Writes fixed_variables.gms for the thermal model. fixON/fixOFF are
% [16 x 24] logical masks in the canonical (stacked) space; since the units
% of a group are interchangeable, layer j of a group maps directly onto the
% j-th physical unit of that group.
fid = fopen(fullfile(dataGmsFolder, 'fixed_variables.gms'), 'w');
fprintf(fid, '* Variables for d = %d\n', d);
nFixed = 0;
for i = 1:16
    for t = 1:24
        if fixON(i,t)
            fprintf(fid, ' u.fx(''%d'',''%d'') = 1;\n', i, t); nFixed = nFixed + 1;
        elseif fixOFF(i,t)
            fprintf(fid, ' u.fx(''%d'',''%d'') = 0;\n', i, t); nFixed = nFixed + 1;
        end
    end
end
fclose(fid);
end
