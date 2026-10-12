
% Self-test for tr_filterFixings. Run it once before wiring the filter into
% the campaign; it should print PASS three times.
%
%  1. The reachability test is checked against exhaustive enumeration of
%     every 8-hour schedule, for all combinations of TU, TD in 1..4 and both
%     initial states, using the constraints exactly as GAMS writes them.
%  2. The 731 real training schedules are fed in fully fixed. They come from
%     feasible MILP solutions, so the filter must not touch a single one.
%  3. Deliberately impossible patterns must come back feasible and as a
%     subset of the input.

rng(0);

%% 1. Reachability vs brute force
H = 8; nBad = 0; nCase = 0;
allSeq = dec2bin(0:2^H-1) - '0';                 % [256 x 8]
for TU = 1:4
    for TD = 1:4
        for u0 = [0 1]
            ok = false(size(allSeq,1),1);
            for k = 1:size(allSeq,1)
                ok(k) = gamsOK(allSeq(k,:), u0, TU, TD);
            end
            feas = allSeq(ok,:);
            for rep = 1:200
                pat = randi([0 3], 1, H);        % 0,1 = free; 2 = off; 3 = on
                fixON  = (pat == 3);
                fixOFF = (pat == 2);
                fix = nan(1,H); fix(fixON) = 1; fix(fixOFF) = 0;
                free = isnan(fix);
                brute = any(all(feas(:,~free) == fix(~free), 2));
                [~,~,dropped] = tr_filterFixings(fixON, fixOFF, TU, TD, u0);
                nCase = nCase + 1;
                if brute ~= (dropped == 0), nBad = nBad + 1; end
            end
        end
    end
end
fprintf('1. Reachability: %d casos, %d discrepancias  -> %s\n', ...
    nCase, nBad, ternary(nBad == 0, 'PASS', 'FALLO'));

%% 2. The real labels must survive untouched
S = load('../results/MIP_Results_2019-2021.mat', 'U');
U = double(S.U);
totalDropped = 0;
for d = 1:size(U,1)
    Ud = squeeze(U(d,:,:)) == 1;                 % [16 x 24]
    [~,~,dropped] = tr_filterFixings(Ud, ~Ud);
    totalDropped = totalDropped + dropped;
end
fprintf('2. Etiquetas reales (%d dias x 16 unidades): %d fijaciones liberadas -> %s\n', ...
    size(U,1), totalDropped, ternary(totalDropped == 0, 'PASS', 'FALLO'));

%% 3. Impossible patterns come back feasible and as a subset
TU = [1 1 7 7 4 4 1 1 1 1 1 1 1 1 6 6];
TD = [1 1 7 7 3 3 1 1 1 1 1 1 1 1 5 5];
nBad = 0; nRep = 0;
for rep = 1:500
    pat = randi([0 3], 16, 24);
    inON = (pat == 3); inOFF = (pat == 2);
    [oON, oOFF, dropped] = tr_filterFixings(inON, inOFF, TU, TD);
    if dropped == 0, continue; end
    nRep = nRep + 1;
    if any(oON & ~inON, "all") || any(oOFF & ~inOFF, "all"), nBad = nBad + 1; continue; end
    [~,~,again] = tr_filterFixings(oON, oOFF, TU, TD);   % must be a fixed point
    if again ~= 0, nBad = nBad + 1; end
end
fprintf('3. Patrones imposibles: %d reparados, %d fallos -> %s\n', ...
    nRep, nBad, ternary(nBad == 0, 'PASS', 'FALLO'));

% ---- Local functions ----
function ok = gamsOK(u, u0, TU, TD)
% MinTimeOnTGU / MinTimeOffTGU as written in data_gms/equations.gms.
uu = [u0 u];
delta = max(0, uu(2:end) - uu(1:end-1));
w     = max(0, uu(1:end-1) - uu(2:end));
ok = true;
for t = 1:numel(u)
    if t >= TU && sum(delta(t-TU+1:t)) > u(t),     ok = false; return; end
    if t >= TD && sum(w(t-TD+1:t))     > 1 - u(t), ok = false; return; end
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
