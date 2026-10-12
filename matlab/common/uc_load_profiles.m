function [DAYDEM, DAYWIND, ND] = uc_load_profiles(period)
%UC_LOAD_PROFILES  Hourly demand and wind profiles of one of the two periods.
%
%   [DAYDEM, DAYWIND, ND] = UC_LOAD_PROFILES(period), with period equal to
%   '2019-2021' (731 days, training) or '2021-2023' (730 days, evaluation).
%   Each output is [days x 24] in MW; ND is the net demand, DAYDEM - DAYWIND,
%   which is the only input of the neural network.
%
%   The series come from the 5-minute records published by Red Eléctrica de
%   España for the Gran Canaria system, averaged to hourly values. See
%   data/README.md for the source and for how to rebuild them.

if nargin < 1, period = '2021-2023'; end
assert(ismember(period, {'2019-2021', '2021-2023'}), ...
    'uc_load_profiles: period must be ''2019-2021'' or ''2021-2023''.');

c = uc_config();
DAYDEM  = readmatrix(fullfile(c.dataRoot, ['demand_' period '.csv']));
DAYWIND = readmatrix(fullfile(c.dataRoot, ['wind_'   period '.csv']));
assert(isequal(size(DAYDEM), size(DAYWIND)), 'Demand and wind series differ in size.');
ND = DAYDEM - DAYWIND;
end
