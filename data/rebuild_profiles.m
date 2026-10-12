
% Rebuilds the demand and wind series of this folder from the 5-minute
% records of Red Eléctrica de España. See README.md for the source.
%
% The .mat must hold a struct Data with fields ts (timestamps), dem (demand,
% MW) and eol (wind generation, MW), sampled every 5 minutes.

rawFile = 'Data_REE_GC_2014_2024.mat';    % <- edit: path to the 5-minute records

here = fileparts(mfilename('fullpath'));
S = load(rawFile, 'Data');  Data = S.Data;

periods = { '2019-2021', datetime('2019-12-19 20:00','InputFormat','yyyy-MM-dd HH:mm')
            '2021-2023', datetime('2021-12-19 20:00','InputFormat','yyyy-MM-dd HH:mm') };
for p = 1:size(periods,1)
    t0 = periods{p,2};
    t1 = t0 + years(2);
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
    writematrix(round(DEM,4), fullfile(here, sprintf('demand_%s.csv', periods{p,1})));
    writematrix(round(WIN,4), fullfile(here, sprintf('wind_%s.csv',   periods{p,1})));
    fprintf('%s: %d dias escritos\n', periods{p,1}, n);
end
