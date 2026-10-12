% Calibration of per-unit confidence bands for fixing UC binary variables
% in the MILP (PHES case). Run AFTER define_NN.m: it uses netTrained,
% XCal, YCal, XTest, YTest, numGen, numHours, idGroups, idxTurb, idxPump,
% unitNames and idxTest from the workspace.
%
% Idea: a variable is fixed only when the network's confidence is so high
% that its error rate on the CALIBRATION set (independent of training and
% of early stopping) stays below a chosen budget, certified with a
% one-sided Clopper-Pearson upper confidence bound.
%
% Unit-aware asymmetric policy: ON fixing for every unit; OFF fixing ONLY
% for the hydro layers. Rationale: wrong ON fixes are cheap and can never
% cause infeasibility; wrong OFF fixes on THERMAL units remove generation
% capacity (the infeasibility mechanism seen in the no-PHES case), so
% thermal OFF stays disabled — but OFF-fixing deep hydro layers is
% physically benign (pumping is optional load, and the thermal fleet
% covers demand without the deep turbine layers), and those layers are
% almost always OFF and near-perfectly predicted.
%
% PHES-specific: a physical reversible unit cannot turbine and pump in the
% same hour, so at most 6 reversible layers in total may be fixed ON per
% hour. When violated, the least confident deep layers are released.
% (Mapping layers to physical units: assign turbine layers from unit 1
% upwards and pump layers from unit 6 downwards — always compatible while
% nT + nP <= 6.)
%
% NOTE: everything lives in the canonical (stacked) space. When applying
% the fixes to the physical MILP, map layer j of each identical group to
% the physical units in a fixed order (they are interchangeable).

% ---- Parameters ----
epsON   = 0.005;                     % ON error budget (all units)
epsOFF  = zeros(1,numGen);           % OFF budgets per unit:
epsOFF([idxTurb idxPump]) = 0.005;   %   hydro layers only, thermal = 0
alphaCP = 0.05;    % confidence level of the certification (95%)
candHi  = [0.50:0.01:0.99, 0.995, 0.999];   % candidate cuts for ON

% ---- Network probabilities on calibration and test ----
nCal  = numel(XCal);
nTest = numel(XTest);
YCalProb  = zeros(numHours,numGen,nCal);
YTestProb = zeros(numHours,numGen,nTest);
for i = 1:nCal,  YCalProb(:,:,i)  = predict(netTrained, XCal{i});  end
for i = 1:nTest, YTestProb(:,:,i) = predict(netTrained, XTest{i}); end

% ---- Main calibration at the chosen budgets ----
[hiBand, loBand] = calibrateBands(YCalProb, YCal, epsON, epsOFF, alphaCP, candHi);
[FixON, FixOFF, nConflict, nCapRel] = applyBands(YTestProb, hiBand, loBand, ...
    idGroups, idxTurb, idxPump);

% ---- Per-unit report on the test set ----
fprintf("Bands certified (CP %.0f%%): epsON = %.3f (all units), epsOFF = %.3f (hydro only)\n\n", ...
    100*(1-alphaCP), epsON, max(epsOFF));
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
fprintf("\nSweep (eps = ON all units = OFF hydro; thermal OFF disabled):\n");
fprintf("  eps     | fixed %%  | realized err %%\n");
for e = [0.05 0.03 0.02 0.01 0.005]
    eOFF = zeros(1,numGen);  eOFF([idxTurb idxPump]) = e;
    [hB, lB] = calibrateBands(YCalProb, YCal, e, eOFF, alphaCP, candHi);
    [fOn, fOff] = applyBands(YTestProb, hB, lB, idGroups, idxTurb, idxPump);
    nF = nnz(fOn) + nnz(fOff);
    nW = nnz(YTest(fOn)==0) + nnz(YTest(fOff)==1);
    fprintf("  %7.3f | %6.2f  | %8.4f\n", e, 100*nF/numel(YTest), 100*nW/max(nF,1));
end

% ---- Export for the MILP pipeline ----
% idxTest maps each test slice back to its original day in the dataset.
save("../data/Fixing_bands.mat", "hiBand","loBand","FixON","FixOFF", ...
     "YTestProb","idxTest","idGroups","idxTurb","idxPump","epsON","epsOFF");
fprintf("\nSaved bands and fixing masks to ../data/Fixing_bands.mat\n");

% ---- Local functions ----
function [hiBand, loBand] = calibrateBands(P, Y, epsON, epsOFF, alpha, candHi)
% Widest per-unit bands whose calibration-set error is within budget,
% certified with a one-sided Clopper-Pearson upper confidence bound.
% epsON/epsOFF may be scalars or per-unit vectors; a zero budget disables
% that direction. NaN band = no cut could be certified: never fix.
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
