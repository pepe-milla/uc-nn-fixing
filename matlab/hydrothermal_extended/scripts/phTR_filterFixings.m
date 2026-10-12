function [fixON, fixOFF, nDropped] = phTR_filterFixings(fixON, fixOFF, TU, TD, u0)
%PHTR_FILTERFIXINGS  Releases thermal fixings that no schedule can satisfy.
%
%   [fixON, fixOFF, nDropped] = PHTR_FILTERFIXINGS(fixON, fixOFF, TU, TD, u0)
%
%   Masks are [28 x 24] logicals in the canonical (stacked) space: rows 1-16
%   thermal, 17-22 turbine, 23-28 pump. Only the thermal rows are examined -
%   the PSH units carry no minimum times, so no pattern of theirs can be
%   infeasible on its own.
%
%   In the simplified formulation every pattern of binaries was feasible, so
%   fixed variables could only make the MILP more expensive. With minimum
%   up/down times that is no longer true: fixing u(i,t) = 1, u(i,t+1) = 0 and
%   u(i,t+2) = 1 is infeasible on its own as soon as TD(i) >= 2, and GAMS
%   then returns INFEASIBLE instead of an expensive solution. The conflict is
%   a property of the combination, not of any single fixing: all three can
%   come from a perfectly calibrated confidence band and the set still be
%   impossible, which is why the bands cannot catch this and this filter has
%   to.
%
%   The function checks each thermal unit against MinTimeOnTGU and
%   MinTimeOffTGU as they are written in data_gms/equations.gms - not the
%   textbook form, which differs at the ends of the horizon - and releases
%   the fewest fixings needed when the check fails. The output is always a
%   subset of the input: an impossible pattern degrades into lower coverage,
%   never into a wrong answer.
%
%   Canonicalization is sound here because the units of each identical
%   thermal group share TU and TD, so a canonical layer and the physical unit
%   it maps onto obey the same constraints.
%
%   Scope: this guarantees that each thermal unit admits some feasible
%   schedule. It says nothing about the power balance, the reserves, the
%   ramps or the reservoir, which couple units and cannot be checked without
%   solving. Those were equally uncheckable in the simplified case; what is
%   new, and what this removes, is the purely combinatorial failure mode.

if nargin < 3 || isempty(TU), TU = [1 1 7 7 4 4 1 1 1 1 1 1 1 1 6 6]; end
if nargin < 4 || isempty(TD), TD = [1 1 7 7 3 3 1 1 1 1 1 1 1 1 5 5]; end
if nargin < 5 || isempty(u0), u0 = [1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 0]; end

assert(~any(fixON & fixOFF, "all"), ...
    'phTR_filterFixings: hay horas marcadas a la vez en fixON y fixOFF.');

numHours = size(fixON, 2);
nDropped = 0;

for i = 1:numel(TU)                       % thermal rows only
    if ~any(fixON(i,:)) && ~any(fixOFF(i,:)), continue; end

    % NaN = free hour, 1 = fixed on, 0 = fixed off
    fix = nan(1, numHours);
    fix(fixON(i,:))  = 1;
    fix(fixOFF(i,:)) = 0;

    while true
        tBad = blockTime(fix, TU(i), TD(i), u0(i));
        if tBad == 0, break; end
        % A free hour can never block on its own (holding the unit in its
        % current state is always admissible), so some earlier fixing is
        % responsible. Releasing the last one at or before the conflict
        % guarantees progress, and with every hour free the constant
        % schedule u = u0 is feasible, so the loop terminates.
        k = find(~isnan(fix(1:tBad)), 1, "last");
        fix(k)   = NaN;
        nDropped = nDropped + 1;
    end

    fixON(i,:)  = (fix == 1);
    fixOFF(i,:) = (fix == 0);
end
end

% ---- Local functions ----
function tBad = blockTime(fix, TU, TD, u0)
% First hour at which no schedule can continue, or 0 if the whole horizon is
% feasible. Forward reachability over the last K+1 values of u packed into a
% bitmask: u(t) in bit 0, u(t-1) in bit 1, ... u(t-K) in bit K. Hours before
% the horizon carry u0, which is the right boundary condition because
% TUR = TDR = 0 in initial_conditions.gms, i.e. no residual obligation is
% inherited from the previous day. K+1 values are enough because delta and w
% at the oldest hour of a window, t-TU+1, need u at t-TU, and TU, TD <= K.
numHours = numel(fix);
K  = max(TU, TD);
nS = 2^(K+1);

reach = false(nS, 1);
if u0, reach(nS) = true; else, reach(1) = true; end   % all K+1 bits = u0
tBad = 0;

for t = 1:numHours
    nxt = false(nS, 1);
    for s = find(reach).'
        st = s - 1;
        for v = 0:1
            if ~isnan(fix(t)) && fix(t) ~= v, continue; end
            ns = bitand(bitshift(st, 1) + v, nS - 1);

            % delta(t-r) = max(0, u(t-r) - u(t-r-1)) = bits r and r+1 of ns
            if t >= TU
                dsum = 0;
                for r = 0:TU-1
                    dsum = dsum + max(0, bitget(ns, r+1) - bitget(ns, r+2));
                end
                if dsum > v, continue; end       % MinTimeOnTGU
            end
            if t >= TD
                wsum = 0;
                for r = 0:TD-1
                    wsum = wsum + max(0, bitget(ns, r+2) - bitget(ns, r+1));
                end
                if wsum > 1 - v, continue; end   % MinTimeOffTGU
            end

            nxt(ns+1) = true;
        end
    end
    reach = nxt;
    if ~any(reach), tBad = t; return; end
end
end
