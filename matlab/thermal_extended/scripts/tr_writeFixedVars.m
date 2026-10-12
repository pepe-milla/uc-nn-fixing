function [nFixed, nDropped] = tr_writeFixedVars(dataGmsFolder, fixON, fixOFF, d, TU, TD, u0)
%TR_WRITEFIXEDVARS  Writes fixed_variables.gms for the TRamp thermal model.
%
%   [nFixed, nDropped] = TR_WRITEFIXEDVARS(dataGmsFolder, fixON, fixOFF, d)
%
%   Same interface and same canonical (stacked) convention as
%   th_writeFixedVars in "thermal_base", with one addition: the fixings
%   are passed through tr_filterFixings first, so that no unit is handed a
%   pattern that the minimum up/down times make impossible. Without that
%   filter GAMS returns INFEASIBLE on those days rather than a costlier
%   solution - a failure mode that did not exist in the simplified
%   formulation. See tr_filterFixings for the reasoning and its scope.
%
%   nDropped is the number of fixings released by the filter. Accumulate it
%   over the campaign: it is the cost of the filter in coverage, and worth
%   reporting alongside the fixing rate.

if nargin < 5, TU = []; end
if nargin < 6, TD = []; end
if nargin < 7, u0 = []; end

[fixON, fixOFF, nDropped] = tr_filterFixings(fixON, fixOFF, TU, TD, u0);

fid = fopen(fullfile(dataGmsFolder, 'fixed_variables.gms'), 'w');
assert(fid > 0, 'tr_writeFixedVars: no se pudo abrir fixed_variables.gms en %s', dataGmsFolder);
fprintf(fid, '* Variables for d = %d\n', d);
nFixed = 0;
for i = 1:size(fixON, 1)
    for t = 1:size(fixON, 2)
        if fixON(i,t)
            fprintf(fid, ' u.fx(''%d'',''%d'') = 1;\n', i, t); nFixed = nFixed + 1;
        elseif fixOFF(i,t)
            fprintf(fid, ' u.fx(''%d'',''%d'') = 0;\n', i, t); nFixed = nFixed + 1;
        end
    end
end
fclose(fid);
end
