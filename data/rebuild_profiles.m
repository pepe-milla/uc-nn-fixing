
% Rebuilds the demand and wind series of this folder from the 5-minute
% records of Red Eléctrica de España. See README.md for the source.
%
% The .mat must hold a struct Data with fields ts (timestamps), dem (demand,
% MW) and eol (wind generation, MW), sampled every 5 minutes.
%
% The values are written at full double precision (%.17g), not rounded: the
% campaigns of the paper fed GAMS with num2str of the unrounded hourly means
% (5 significant digits below 1 MW) and the network with the unrounded net
% demand, so any rounding here makes some days solve differently.

rawFile = 'Data_REE_GC_2014_2024.mat';    % <- edit: path to the 5-minute records

here = fileparts(mfilename('fullpath'));
S = load(rawFile, 'Data');  Data = S.Data;

% Same start and end instants as the campaigns (731 and 730 days)
periods = { '2019-2021', '2019-12-19 20:00', '2021-12-19 20:00'
            '2021-2023', '2021-12-19 20:00', '2023-12-19 20:00' };
for p = 1:size(periods,1)
    t0 = datetime(periods{p,2}, 'InputFormat','yyyy-MM-dd HH:mm');
    t1 = datetime(periods{p,3}, 'InputFormat','yyyy-MM-dd HH:mm');
    n  = days(duration(t1 - t0));
    r  = 1 + 12*hours(t0 - Data.ts(1));
    DEM = zeros(n,24);  WIN = zeros(n,24);
    for d = 1:n
        for h = 1:24
            DEM(d,h) = mean(Data.dem(r:(r+11)));   % 12 records of 5 min = 1 hour
            WIN(d,h) = mean(Data.eol(r:(r+11)));
            r = r + 12;
        end
    end
    writeFull(DEM, fullfile(here, sprintf('demand_%s.csv', periods{p,1})));
    writeFull(WIN, fullfile(here, sprintf('wind_%s.csv',   periods{p,1})));
    fprintf('%s: %d dias escritos\n', periods{p,1}, n);
end

function writeFull(M, file)
% One row per day, comma-separated, every value with 17 significant digits
% so that readmatrix returns exactly the same doubles.
fid = fopen(file, 'w');
assert(fid ~= -1, 'No se pudo escribir %s', file);
fmt = [repmat('%.17g,', 1, size(M,2)-1) '%.17g\n'];
fprintf(fid, fmt, M.');
fclose(fid);
end
