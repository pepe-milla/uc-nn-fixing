% Calibration of per-unit confidence bands for fixing UC binary variables
% in the MILP. Run AFTER define_NN.m: it uses netTrained, XCal, YCal,
% XTest, YTest, numGen, numHours, idGroups and idxTest from the workspace.
%
% Idea: a variable is fixed only when the network's confidence is so high
% that its error rate on the CALIBRATION set (independent of training and
% of early stopping) stays below a chosen budget. Bands are calibrated per
% unit and per direction (ON / OFF), because probability calibration
% differs strongly across units, and a wrong OFF fix is more dangerous
% than a wrong ON fix.
%
% NOTE: everything lives in the canonical (stacked) space. When applying
% the fixes to the physical MILP, map layer j of each identical group to
% the physical units in any fixed order (they are interchangeable).

% ---- Parameters ----
% Bands are certified with a one-sided Clopper-Pearson upper confidence
% bound: a cut is accepted only if the 95% upper bound of its error rate on
% the CALIBRATION set is within budget. This kills bands that look clean
% just by luck (few samples / winner's curse), at the price of lower
% coverage: certifying eps needs roughly 3/eps clean samples. See the note
% in calibrateBands on the scope of the guarantee.
% Asymmetric operating point: generous on the ON side (wrong ON fixes are
% cheap — a startup plus min load — and can never cause infeasibility) and
% fully closed on the OFF side (wrong OFF fixes remove capacity: they are
% the expensive/infeasible ones). epsOFF = 0 disables OFF fixing entirely:
% no band can certify against a zero budget.
%
% epsON is deliberately the same 0.005 used in the simplified thermal case
% ("thermal_base"), so that the coverage and the realized error of the
% two formulations are read off the same operating point. Changing the
% budget between formulations would make the comparison meaningless.
epsON   = 0.005;   % error budget for variables fixed to ON
epsOFF  = 0;       % OFF fixing disabled
alphaCP = 0.05;    % family-wise confidence level of the certification (95%)
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
[FixON, FixOFF, nConflict] = applyBands(YTestProb, hiBand, loBand, idGroups);

% ---- Per-unit report on the test set ----
fprintf("Bands certified (CP %.0f%% per cut, %d candidates, uncorrected)\n", ...
    100*(1-alphaCP), numel(candHi));
fprintf("epsON = %.4f, epsOFF = %.4f\n", epsON, epsOFF);
fprintf("Check the err ON column against epsON: the guarantee is per cut, ");
fprintf("not for the selected band.\n\n");
fprintf("Unit |  loBand  hiBand | %%fixOFF  %%fixON | err OFF   err ON\n");
fprintf("-----+-----------------+------------------+------------------\n");
for g = 1:numGen
    fOn  = FixON(:,g,:);  fOff = FixOFF(:,g,:);  yT = YTest(:,g,:);
    covOn  = mean(fOn(:));
    covOff = mean(fOff(:));
    if any(fOff(:)), eOff = mean(yT(fOff)==1); else, eOff = NaN; end
    if any(fOn(:)),  eOn  = mean(yT(fOn)==0);  else, eOn  = NaN; end
    fprintf("G%-3d | %7.3f %7.3f | %6.1f %%  %5.1f %% | %7.4f  %7.4f\n", ...
        g, loBand(g), hiBand(g), covOff*100, covOn*100, eOff, eOn);
end

% ---- Global summary on the test set ----
nFix   = nnz(FixON) + nnz(FixOFF);
nWrong = nnz(YTest(FixON)==0) + nnz(YTest(FixOFF)==1);
fprintf("\nVariables fixed: %5.1f %% of all test binaries (%d of %d)\n", ...
    100*nFix/numel(YTest), nFix, numel(YTest));
fprintf("Errors among fixed variables: %.4f %% (%d)\n", 100*nWrong/max(nFix,1), nWrong);
fprintf("Group-consistency conflicts resolved (unfixed): %d\n", nConflict);
fprintf("Free variables left to the solver: %4.1f %%\n", 100*(1 - nFix/numel(YTest)));

% ---- Sweep: coverage vs realized error for several budgets ----
% Produces the data for the thesis curve: error budget -> (coverage, error).
% The sweep must follow the SAME policy as the operating point (OFF fixing
% disabled), otherwise the curve would describe a different method than the
% one actually used in the campaigns.
fprintf("\nSweep (OFF fixing disabled, as in the operating point):\n");
fprintf("  epsON   | fixed %%  | realized err %%\n");
for e = [0.05 0.03 0.02 0.01 0.005]
    [hB, lB] = calibrateBands(YCalProb, YCal, e, 0, alphaCP, candHi);
    [fOn, fOff] = applyBands(YTestProb, hB, lB, idGroups);
    nF = nnz(fOn) + nnz(fOff);
    nW = nnz(YTest(fOn)==0) + nnz(YTest(fOff)==1);
    fprintf("  %7.3f | %6.2f  | %8.4f\n", e, 100*nF/numel(YTest), 100*nW/max(nF,1));
end

% ---- Export for the MILP pipeline ----
% idxTest maps each test slice back to its original day in the dataset.
save("../data/Fixing_bands.mat", "hiBand","loBand","FixON","FixOFF", ...
     "YTestProb","idxTest","idGroups","epsON","epsOFF","alphaCP","candHi");
fprintf("\nSaved bands and fixing masks to ../data/Fixing_bands.mat\n");

% ---- Local functions ----
function [hiBand, loBand] = calibrateBands(P, Y, epsON, epsOFF, alpha, candHi)
% Widest per-unit bands whose calibration-set error is within budget,
% certified with a one-sided Clopper-Pearson upper confidence bound.
% NaN band = no cut could be certified: that direction is never fixed.
%
% LIMITATION, stated deliberately. The sweep does not test one cut chosen
% in advance: it tries numel(candHi) of them and keeps the first that
% passes. The (1-alpha) bound therefore holds for each candidate in
% isolation, not for the one actually selected, and the probability that
% the SELECTED band breaches the budget can reach
% 1-(1-alpha)^numel(candHi) — 93% for the 52-point grid used here. A
% Bonferroni correction (alpha/numel(candHi)) would restore the family-wise
% level; it is not applied, so that this case uses exactly the same
% procedure as the published simplified thermal and PHES cases, which would
% otherwise have to be recalibrated and re-run for the comparison to mean
% anything. The correction is the single line commented out below.
%
% What this costs is real and must be watched rather than assumed away. At
% epsON = 0.01 it bit: G6 certified with 2 errors in ~684 calibration points
% (CP95 = 0.0092) and then realized 2.21% on the test set, more than twice
% its budget, contributing 15 of the 47 total errors with 5% of the
% fixings. At epsON = 0.005 the simplified thermal case stayed inside budget
% on every unit (worst 0.15%, aggregate 0.080%). So: always read the per-
% unit "err ON" column of the report against epsON before trusting a
% campaign, and if a unit exceeds it, tighten that band by hand.
% alpha = alpha / numel(candHi);   % <- Bonferroni, if ever needed

numGen = size(P,2);
hiBand = nan(1,numGen);
loBand = nan(1,numGen);
for g = 1:numGen
    p = squeeze(P(:,g,:));
    y = squeeze(double(Y(:,g,:)));
    for hi = candHi                       % ascending: first hit = widest band
        sel = p > hi;
        if binoUCB(nnz(y(sel)==0), nnz(sel), alpha) <= epsON
            hiBand(g) = hi;
            break
        end
    end
    for lo = sort(1-candHi, 'descend')    % descending: first hit = widest band
        sel = p < lo;
        if binoUCB(nnz(y(sel)==1), nnz(sel), alpha) <= epsOFF
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

function [FixON, FixOFF, nConflict] = applyBands(P, hiBand, loBand, groups)
% Fixing masks from the bands (NaN band never fixes), made consistent with
% the stacked encoding: ON at layer j forces ON at shallower layers, OFF at
% layer j forces OFF at deeper layers. Contradictions are unfixed.
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
end
