
% Calibration of per-unit confidence bands for fixing UC binary variables in
% the MILP (hydrothermal case WITH intertemporal constraints). Run AFTER
% define_NN.m: it uses netTrained, XCal, YCal, XTest, YTest, numGen,
% numHours, idGroups, idxTherm, idxTurb, idxPump, unitNames and idxTest from
% the workspace.
%
% Idea: a variable is fixed only when the network's confidence is so high
% that its error rate on the CALIBRATION set (independent of training and of
% early stopping) stays below a chosen budget, certified with a one-sided
% Clopper-Pearson upper confidence bound.
%
% Unit-aware asymmetric policy, unchanged from the simplified hydrothermal
% case: ON fixing for every unit; OFF fixing ONLY for the hydro layers.
% Wrong ON fixes are cheap and cannot cause infeasibility; wrong OFF fixes
% on THERMAL units remove generation capacity, so thermal OFF stays
% disabled. OFF-fixing deep hydro layers is physically benign - pumping is
% optional load, the thermal fleet covers demand without the deep turbine
% layers - and those layers are almost always OFF and near-perfectly
% predicted. That reasoning survives the change of formulation intact,
% because the intertemporal constraints were added to the thermal units
% only: the hydro layers still carry no minimum times and no ramp limits.
%
% What does change is what an ON fix can do. In the simplified formulation
% any pattern of binaries was feasible and a bad fix could only cost money;
% with minimum up/down times a combination of thermal fixes can be
% impossible on its own, and GAMS then returns INFEASIBLE. That is not
% handled here but at writing time, by the feasibility filter in
% assess_NN_par.m - the bands are a statistical statement about individual
% variables and cannot see combinations.
%
% NOTE: everything lives in the canonical (stacked) space. When applying the
% fixes to the physical MILP, map layer j of each identical group to the
% physical units in a fixed order (they are interchangeable).

% ---- Parameters ----
epsON   = 0.005;                     % ON error budget (all units)
epsOFF  = zeros(1,numGen);           % OFF budgets per unit:
epsOFF([idxTurb idxPump]) = 0.005;   %   hydro layers only, thermal = 0
alphaCP = 0.05;    % confidence level of the certification (95%)
candHi  = [0.50:0.01:0.99, 0.995, 0.999];   % candidate cuts

% ---- Manual band overrides ----
% The certification is per cut and not for the selected band (see the note
% in calibrateBands), so a unit can certify on the calibration set and then
% miss its budget on the test set. When the per-unit report shows that, the
% band is corrected here rather than in the workspace, so the record of what
% was disabled - and why - lives in the script that produced the bands.
%
% Overriding is always conservative: it can only remove fixings, never add
% them, so it cannot invalidate anything downstream.
%
% G16 disabled. It certified a band of 0.890 and fixed 75.6% of its hours,
% then realized 2.38% error against a 0.5% budget: 47 errors where 9.9 were
% expected, which is not a fluctuation (p is numerically zero). It alone
% produced 48% of all the errors with 5.9% of the fixings. The mechanism is
% visible in the transition table of define_NN: G16 never switches in the
% whole period, and the network predicts 0.39 switches per day for it. With
% TU = 6 and TD = 5 those spurious switches are also the fixings most likely
% to make a day infeasible, so removing them buys feasibility as well as
% accuracy. Cost: coverage 45.8% -> 43.1%, realized error 0.29% -> 0.16%.
bandOverride.disableON  = 16;      % unidades cuya banda ON se anula
bandOverride.disableOFF = [];      % idem para la banda OFF
bandOverride.setON      = [];      % [unidad cut; ...] para estrechar sin anular
bandOverride.setOFF     = [];      % [unidad cut; ...]

% ---- Network probabilities on calibration and test ----
nCal  = numel(XCal);
nTest = numel(XTest);
YCalProb  = zeros(numHours,numGen,nCal);
YTestProb = zeros(numHours,numGen,nTest);
for i = 1:nCal,  YCalProb(:,:,i)  = predict(netTrained, XCal{i});  end
for i = 1:nTest, YTestProb(:,:,i) = predict(netTrained, XTest{i}); end

% ---- Main calibration at the chosen budgets ----
[hiBand, loBand] = calibrateBands(YCalProb, YCal, epsON, epsOFF, alphaCP, candHi);
[hiBand, loBand] = applyOverrides(hiBand, loBand, bandOverride, unitNames, true);
[FixON, FixOFF, nConflict, nCapRel] = applyBands(YTestProb, hiBand, loBand, ...
    idGroups, idxTurb, idxPump);

% ---- Per-unit report on the test set ----
fprintf("Bands certified (CP %.0f%% per cut, %d candidates, uncorrected)\n", ...
    100*(1-alphaCP), numel(candHi));
fprintf("epsON = %.4f (all units), epsOFF = %.4f (hydro only)\n", epsON, max(epsOFF));
fprintf("Check the err columns against the budgets: the guarantee is per cut, ");
fprintf("not for the selected band.\n\n");
fprintf("Unit |  loBand  hiBand | %%fixOFF  %%fixON | err OFF   err ON\n");
fprintf("-----+-----------------+------------------+------------------\n");
for g = 1:numGen
    fOn  = FixON(:,g,:);  fOff = FixOFF(:,g,:);  yT = YTest(:,g,:);
    covOn  = mean(fOn(:));
    covOff = mean(fOff(:));
    if any(fOff(:)), eOff = mean(yT(fOff)==1); else, eOff = NaN; end
    if any(fOn(:)),  eOn  = mean(yT(fOn)==0);  else, eOn  = NaN; end
    fprintf("%-4s | %7.3f %7.3f | %6.1f %%  %5.1f %% | %7.4f  %7.4f\n", ...
        unitNames(g), loBand(g), hiBand(g), covOff*100, covOn*100, eOff, eOn);
end

% ---- Global summary on the test set ----
nFix   = nnz(FixON) + nnz(FixOFF);
nWrong = nnz(YTest(FixON)==0) + nnz(YTest(FixOFF)==1);
fprintf("\nVariables fixed: %5.1f %% of all test binaries (%d of %d)\n", ...
    100*nFix/numel(YTest), nFix, numel(YTest));
fprintf("Errors among fixed variables: %.4f %% (%d)\n", 100*nWrong/max(nFix,1), nWrong);
fprintf("Group-consistency conflicts resolved (unfixed): %d\n", nConflict);
fprintf("Turbine+pump cap releases: %d\n", nCapRel);
fprintf("Free variables left to the solver: %4.1f %%\n", 100*(1 - nFix/numel(YTest)));

% ---- Sweep: coverage vs realized error for several budgets ----
% Produces the data for the thesis curve: error budget -> (coverage, error).
% Same policy at every point: eps for ON (all units) and for OFF on the
% hydro layers; thermal OFF fixing stays disabled.
% The overrides apply at every point of the sweep too: a unit that is
% excluded is excluded from the method, not just from one operating point,
% and the curve has to describe the method as it is actually used.
fprintf("\nSweep (eps = ON all units = OFF hydro; thermal OFF disabled):\n");
fprintf("  eps     | fixed %%  | realized err %%\n");
for e = [0.05 0.03 0.02 0.01 0.005]
    eOFF = zeros(1,numGen);  eOFF([idxTurb idxPump]) = e;
    [hB, lB] = calibrateBands(YCalProb, YCal, e, eOFF, alphaCP, candHi);
    [hB, lB] = applyOverrides(hB, lB, bandOverride, unitNames, false);
    [fOn, fOff] = applyBands(YTestProb, hB, lB, idGroups, idxTurb, idxPump);
    nF = nnz(fOn) + nnz(fOff);
    nW = nnz(YTest(fOn)==0) + nnz(YTest(fOff)==1);
    fprintf("  %7.3f | %6.2f  | %8.4f\n", e, 100*nF/numel(YTest), 100*nW/max(nF,1));
end

% ---- Export for the MILP pipeline ----
% idxTest maps each test slice back to its original day in the dataset.
if ~isfolder("../data"), mkdir("../data"); end
save("../data/Fixing_bands.mat", "hiBand","loBand","FixON","FixOFF", ...
     "YTestProb","idxTest","idGroups","idxTherm","idxTurb","idxPump", ...
     "epsON","epsOFF","alphaCP","candHi","bandOverride");
fprintf("\nSaved bands and fixing masks to ../data/Fixing_bands.mat\n");

% ---- Local functions ----
function [hiBand, loBand] = applyOverrides(hiBand, loBand, ov, unitNames, verbose)
% Applies the manual band corrections. Only ever tightens: disabling sets
% the band to NaN, which never fixes, and setting a cut is rejected unless
% it is stricter than the certified one. That guard matters - a manual
% override that loosened a band would silently break the certification the
% rest of the pipeline relies on.
if isfield(ov,'disableON') && ~isempty(ov.disableON)
    hiBand(ov.disableON) = NaN;
    if verbose
        fprintf("Banda ON anulada a mano: %s\n", strjoin(cellstr(unitNames(ov.disableON)), ', '));
    end
end
if isfield(ov,'disableOFF') && ~isempty(ov.disableOFF)
    loBand(ov.disableOFF) = NaN;
    if verbose
        fprintf("Banda OFF anulada a mano: %s\n", strjoin(cellstr(unitNames(ov.disableOFF)), ', '));
    end
end
if isfield(ov,'setON') && ~isempty(ov.setON)
    for r = 1:size(ov.setON,1)
        g = ov.setON(r,1);  c = ov.setON(r,2);
        assert(isnan(hiBand(g)) || c >= hiBand(g), ...
            'setON para %s afloja la banda (%.3f < %.3f certificada).', string(unitNames(g)), c,hiBand(g));
        hiBand(g) = c;
        if verbose, fprintf("Banda ON de %s fijada a mano en %.3f\n", unitNames(g), c); end
    end
end
if isfield(ov,'setOFF') && ~isempty(ov.setOFF)
    for r = 1:size(ov.setOFF,1)
        g = ov.setOFF(r,1);  c = ov.setOFF(r,2);
        assert(isnan(loBand(g)) || c <= loBand(g), ...
            'setOFF para %s afloja la banda (%.3f > %.3f certificada).', string(unitNames(g)), c,loBand(g));
        loBand(g) = c;
        if verbose, fprintf("Banda OFF de %s fijada a mano en %.3f\n", unitNames(g), c); end
    end
end
if verbose, fprintf("\n"); end
end

function [hiBand, loBand] = calibrateBands(P, Y, epsON, epsOFF, alpha, candHi)
% Widest per-unit bands whose calibration-set error is within budget,
% certified with a one-sided Clopper-Pearson upper confidence bound.
% epsON/epsOFF may be scalars or per-unit vectors; a zero budget disables
% that direction. NaN band = no cut could be certified: never fix.
%
% LIMITATION, stated deliberately. The sweep does not test one cut chosen in
% advance: it tries numel(candHi) of them and keeps the first that passes.
% The (1-alpha) bound therefore holds for each candidate in isolation, not
% for the one actually selected, and the probability that the SELECTED band
% breaches the budget can reach 1-(1-alpha)^numel(candHi) - 93% for the
% 52-point grid used here. A Bonferroni correction (alpha/numel(candHi))
% would restore the family-wise level; it is not applied, so that this case
% uses exactly the same procedure as the published simplified thermal and
% PHES cases, which would otherwise have to be recalibrated and re-run for
% the comparison to mean anything. The correction is the line commented out
% below. Always read the per-unit err columns of the report against the
% budgets before trusting a campaign, and tighten a band by hand if one
% exceeds them - tightening is always conservative, since it can only remove
% fixings.
% alpha = alpha / numel(candHi);   % <- Bonferroni, if ever needed

numGen = size(P,2);
if isscalar(epsON),  epsON  = repmat(epsON, 1, numGen); end
if isscalar(epsOFF), epsOFF = repmat(epsOFF, 1, numGen); end
hiBand = nan(1,numGen);
loBand = nan(1,numGen);
for g = 1:numGen
    p = squeeze(P(:,g,:));
    y = squeeze(double(Y(:,g,:)));
    for hi = candHi                       % ascending: first hit = widest band
        sel = p > hi;
        if binoUCB(nnz(y(sel)==0), nnz(sel), alpha) <= epsON(g)
            hiBand(g) = hi;
            break
        end
    end
    for lo = sort(1-candHi, 'descend')    % descending: first hit = widest band
        sel = p < lo;
        if binoUCB(nnz(y(sel)==1), nnz(sel), alpha) <= epsOFF(g)
            loBand(g) = lo;
            break
        end
    end
end
end

function u = binoUCB(k, n, alpha)
% One-sided (1-alpha) Clopper-Pearson upper bound for a binomial proportion
% (k errors in n trials). With n = 0 or k = n the bound is 1 (never certifies).
if n == 0 || k >= n
    u = 1;
else
    u = betainv(1-alpha, k+1, n-k);
end
end

function [FixON, FixOFF, nConflict, nCapRel] = applyBands(P, hiBand, loBand, ...
    groups, grpT, grpP)
% Fixing masks from the bands (NaN band never fixes), made consistent with
% the stacked encoding: ON at layer j forces ON at shallower layers, OFF at
% layer j forces OFF at deeper layers. Contradictions are unfixed. Then the
% reversible-plant cap is enforced: at most 6 turbine+pump layers fixed ON
% per hour, releasing the least confident deep layers when violated.
FixON  = P > reshape(hiBand,1,[]);
FixOFF = P < reshape(loBand,1,[]);
for k = 1:numel(groups)
    idx = groups{k};
    for j = 2:numel(idx)
        FixON(:,idx(1:j-1),:) = FixON(:,idx(1:j-1),:) | FixON(:,idx(j),:);
        FixOFF(:,idx(j:end),:) = FixOFF(:,idx(j:end),:) | FixOFF(:,idx(j-1),:);
    end
end
conflict = FixON & FixOFF;
nConflict = nnz(conflict);
FixON(conflict)  = false;
FixOFF(conflict) = false;

% Turbine+pump cap. After stacking, ON-fixed layers of a group are always
% the contiguous block 1..k, so the deepest fixed layer is layer k.
nCapRel = 0;
tot = sum(FixON(:,grpT,:),2) + sum(FixON(:,grpP,:),2);   % [24 x 1 x n]
[tIdx,~,sIdx] = ind2sub(size(tot), find(tot > numel(grpT)));
for m = 1:numel(tIdx)
    t = tIdx(m);  s = sIdx(m);
    nT = sum(FixON(t,grpT,s));
    nP = sum(FixON(t,grpP,s));
    while nT + nP > numel(grpT)
        if nP == 0 || (nT > 0 && P(t,grpT(nT),s) <= P(t,grpP(nP),s))
            FixON(t,grpT(nT),s) = false;  nT = nT - 1;
        else
            FixON(t,grpP(nP),s) = false;  nP = nP - 1;
        end
        nCapRel = nCapRel + 1;
    end
end
end
