
% Reproduces every table and figure of the paper from the stored results,
% without solving a single MILP. Run it from anywhere:
%
%     cd <repo>/reproduce ; reproduce_all
%
% It needs MATLAB with the Deep Learning Toolbox (to evaluate the stored
% networks) and the Statistics Toolbox (for the confidence bounds and the
% rank correlations). It takes a couple of minutes, most of it spent
% evaluating the networks over the 730 days of the evaluation period.
%
% What it prints, in order:
%   Table 5  classification metrics of the four networks
%   Table 6  the same metrics without the canonical encoding
%   Table 7  the network against two trivial rules
%   Table 8  acceleration of the solver (Strategies 1, 2 and 3)
%   Table 9  the three budgets of the thermal simplified case
%   Table 10 realized error of the fixed variables, out of sample
%   Table 11 wrong fixings against solver performance
%   Table 12 solver statistics
%   Tables 13 and 14 (Annex A) per-unit statistics
% and it writes Figures 3 and 4 into reproduce/figures.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'matlab', 'common'));
c = uc_config();
M = c.matlabRoot;
figDir = fullfile(here, 'figures');
if ~isfolder(figDir), mkdir(figDir); end
addpath(fullfile(M, 'thermal_extended', 'scripts'), fullfile(M, 'hydrothermal_extended', 'scripts'));

gT = {[1 2],[3 4],[7 8],[9 10],[11 12 13],[15 16]};
TU = [1 1 7 7 4 4 1 1 1 1 1 1 1 1 6 6];
TD = [1 1 7 7 3 3 1 1 1 1 1 1 1 1 5 5];
u0 = [1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 0];

% name | folder | net | bands | labels | baseline | S2 | epgap | tilim | nG | extended
C = { 'Thermal, simplified',      'thermal_base',          'Fixing_bands_eps005.mat', 'data/Training_Data_2019-2021.mat',    'MIP_Baseline_thermal_2021-2023.mat', 'Assess_NN_S2_thermal_eps005.mat', 0.001, 120,  16, false
      'Thermal, extended',        'thermal_extended',      'Fixing_bands.mat',        'data/Training_Data_2019-2021.mat',    'MIP_Baseline_2021-2023_i9.mat',      'Assess_NN_S2_TRamp_i9.mat',       0.001, 120,  16, true
      'Hydrothermal, simplified', 'hydrothermal_base',     'Fixing_bands_S2.mat',     'data/Training_Data.mat',              'MIP_Baseline_2021-2023.mat',         'Assess_NN_par.mat',               0.015, 1200, 28, false
      'Hydrothermal, extended',   'hydrothermal_extended', 'Fixing_bands.mat',        'results/MIP_Results_2019-2021.mat',   'MIP_Baseline_2021-2023.mat',         'Assess_NN_S2_PHES_TRamp.mat',     0.015, 1200, 28, true };
nC = size(C,1);
S = struct();

fprintf('\nEvaluando las redes sobre el periodo de evaluacion...\n');
for k = 1:nC
    d = fullfile(M, C{k,2});
    S(k).name = C{k,1};  S(k).nG = C{k,9};   S(k).ext = C{k,10};
    S(k).epgap = C{k,7}; S(k).tilim = C{k,8};
    S(k).G = gT;  if S(k).nG == 28, S(k).G = [gT, {17:22}, {23:28}]; end
    N = load(fullfile(d, 'trained models', 'NN.mat'), 'netTrained');
    S(k).net = N.netTrained;
    B = load(fullfile(d, 'data', C{k,3}), 'hiBand', 'loBand');
    S(k).hiBand = B.hiBand;  S(k).loBand = B.loBand;
    S(k).base = load(fullfile(d, 'results', C{k,5}));
    S(k).s2   = load(fullfile(d, 'results', C{k,6}));
    S(k).lab  = load(fullfile(d, strrep(C{k,4}, '/', filesep)));
    nd = size(S(k).base.ND, 1);
    P = zeros(nd, S(k).nG, 24);
    for dday = 1:nd
        y = extractdata(predict(S(k).net, dlarray(single(S(k).base.ND(dday,:)'), 'TC')));
        if size(y,1) ~= S(k).nG, y = y'; end
        P(dday,:,:) = double(y);
    end
    S(k).P = P;
    Y = round(S(k).base.U);
    if S(k).nG == 28, Y = cat(2, Y, round(S(k).base.Ug), round(S(k).base.Up)); end
    S(k).Yraw = Y;  S(k).Y = canon3(Y, S(k).G);
    S(k).ok = logical(S(k).base.OK(:));
    S(k).closed = S(k).ok & S(k).base.OPT_GAP(:) <= S(k).epgap*1.001;
    fprintf('  %-26s %d dias, %d resueltos, %d cerrados\n', S(k).name, nd, nnz(S(k).ok), nnz(S(k).closed));
end

%% ---------------- Table 5: classification ------------------------------
fprintf('\n=== Table 5: classification metrics (threshold 0.5, canonical labels)\n');
fprintf('%-28s %6s %6s %6s %6s\n', 'Case / group', 'Acc', 'Prec', 'Rec', 'F1');
for k = 1:nC
    Yh = canon3(S(k).P > 0.5, S(k).G);
    S(k).Yhat = Yh;
    m = metrics(Yh(S(k).ok,:,:), S(k).Y(S(k).ok,:,:));
    fprintf('%-28s %6.3f %6.3f %6.3f %6.3f\n', S(k).name, m.acc, m.prec, m.rec, m.f1);
    if S(k).nG == 28
        grp = {1:16, 'thermal units'; 17:22, 'turbine mode'; 23:28, 'pump mode'};
        for q = 1:3
            m = metrics(Yh(S(k).ok,grp{q,1},:), S(k).Y(S(k).ok,grp{q,1},:));
            fprintf('   %-25s %6.3f %6.3f %6.3f %6.3f\n', grp{q,2}, m.acc, m.prec, m.rec, m.f1);
        end
    end
end
%% ---------------- Table 6: the same metrics, no canonical encoding -----
% Neither the predictions nor the labels are canonicalized here, i.e. the
% labels are taken exactly as the solver returned them. The gap against
% Table 5 is the effect of the label degeneracy of the identical units.
fprintf('\n=== Table 6: classification metrics without the canonical encoding\n');
fprintf('%-28s %6s %6s %6s %6s\n', 'Case / group', 'Acc', 'Prec', 'Rec', 'F1');
for k = 1:nC
    Yr = S(k).P > 0.5;
    m = metrics(Yr(S(k).ok,:,:), S(k).Yraw(S(k).ok,:,:));
    fprintf('%-28s %6.3f %6.3f %6.3f %6.3f\n', S(k).name, m.acc, m.prec, m.rec, m.f1);
    if S(k).nG == 28
        grp = {1:16, 'thermal units'; 17:22, 'turbine mode'; 23:28, 'pump mode'};
        for q = 1:3
            m = metrics(Yr(S(k).ok,grp{q,1},:), S(k).Yraw(S(k).ok,grp{q,1},:));
            fprintf('   %-25s %6.3f %6.3f %6.3f %6.3f\n', grp{q,2}, m.acc, m.prec, m.rec, m.f1);
        end
    end
end

%% ---------------- Table 7 and Annex A: trivial rules --------------------
fprintf('\n=== Table 7: the network against two trivial rules\n');
fprintf('%-28s %18s %18s\n', '', 'Accuracy', 'F1-score');
fprintf('%-28s %6s %6s %6s %6s %6s %6s\n', 'Case', 'Const', 'Thr', 'NN', 'Const', 'Thr', 'NN');
for k = 1:nC
    [U, cons, thr] = trivialRules(S(k), gT);
    S(k).unitTab = U;
    Y = S(k).Y(S(k).ok,:,:);  Yh = S(k).Yhat(S(k).ok,:,:);
    f1 = @(a,b) 2*nnz(a & b)/max(1, 2*nnz(a & b) + nnz(a & ~b) + nnz(~a & b));
    fprintf('%-28s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n', S(k).name, ...
        mean(cons(:)==Y(:)), mean(thr(:)==Y(:)), mean(Yh(:)==Y(:)), f1(cons,Y), f1(thr,Y), f1(Yh,Y));
end

%% ---------------- Table 8: acceleration --------------------------------
fprintf('\n=== Table 8: acceleration of the MIP solver\n');
fprintf('%-26s %-12s %6s %8s %8s %8s %8s %8s %7s %6s %8s\n', 'Case', 'Method', 'Fixed', 'mean s', 'med s', 'total h', 'aggr', 'med x', 'faster', 'open', 'cost %');
extra = { 'thermal_base', {'Assess_NN_S1_thermal.mat','Strategy 1'; 'Assess_NN_S3_thermal.mat','Strategy 3'}
          'hydrothermal_base', {'Assess_NN_S1_par.mat','Strategy 1'; 'Assess_NN_S3_par.mat','Strategy 3'} };
for k = 1:nC
    b = S(k).base;  both = S(k).ok;
    fprintf('%-26s %-12s %6s %8.2f %8.2f %8.1f %8s %8s %7s %6d %8s\n', S(k).name, 'Baseline', '--', ...
        mean(b.COMP_TIME(both)), median(b.COMP_TIME(both)), sum(b.COMP_TIME(both))/3600, ...
        '--', '--', '--', nnz(b.OPT_GAP(both) > S(k).epgap*1.001), '--');
    rows = {};
    idx = find(strcmp(extra(:,1), C{k,2}));
    if ~isempty(idx)
        for r = 1:size(extra{idx,2},1)
            f = fullfile(M, C{k,2}, 'results', extra{idx,2}{r,1});
            if isfile(f), rows(end+1,:) = {load(f), extra{idx,2}{r,2}}; end %#ok<SAGROW>
        end
    end
    rows(end+1,:) = {S(k).s2, 'Strategy 2'}; %#ok<SAGROW>
    [~, o] = sort(cellfun(@(x) string(x), rows(:,2)));
    rows = rows(o,:);
    for r = 1:size(rows,1)
        printSpeed(rows{r,2}, b, rows{r,1}, S(k), M, C{k,2});
    end
end

%% ---------------- Table 9: budgets (thermal simplified) ----------------
fprintf('\n=== Table 9: the three budgets of the thermal simplified case\n');
fprintf('%6s %8s %8s %10s %10s %9s %7s\n', 'eps %', 'fixed %', 'error %', 'aggregate', 'median', 'cost %', 'open');
camp = { 0.005, 'Assess_NN_S2_thermal_eps005.mat', 'Fixing_bands_eps005.mat'
         0.010, 'Assess_NN_S2_thermal_eps01.mat',  'Fixing_bands_eps01.mat'
         0.030, 'Assess_NN_S2_thermal_eps03.mat',  'Fixing_bands_eps03.mat' };
k = 1;  % thermal simplified
for r = 1:3
    Sc = load(fullfile(M,'thermal_base','results',camp{r,2}));
    Bd = load(fullfile(M,'thermal_base','data',camp{r,3}), 'hiBand','loBand');
    [cov, err] = bandsOnPeriod(S(k), Bd.hiBand, Bd.loBand, TU, TD, u0);
    a = speedStats(S(k).base, Sc, S(k).ok, S(k).epgap);
    fprintf('%6.1f %8.1f %8.3f %9.2fx %9.2fx %+9.3f %7d\n', 100*camp{r,1}, 100*cov, 100*err, ...
        a.agg, a.medRatio, a.costClosed, a.truncS);
end

%% ---------------- Table 10 and Table 11 --------------------------------
fprintf('\n=== Table 10: realized error of the fixed variables (eps = 0.5%%)\n');
fprintf('%-28s %6s %12s %8s %8s %8s\n', 'Case', 'days', 'fixings', 'err %', 'ON %', 'OFF %');
for k = 1:nC
    [~, ~, U] = bandsOnPeriod(S(k), S(k).hiBand, S(k).loBand, TU, TD, u0);
    S(k).errUnits = U;
    fprintf('%-28s %6d %12d %8.3f %8.3f %8s\n', S(k).name, nnz(S(k).closed), ...
        sum(U.nOn)+sum(U.nOff), 100*(sum(U.eOn)+sum(U.eOff))/max(1,sum(U.nOn)+sum(U.nOff)), ...
        100*sum(U.eOn)/max(1,sum(U.nOn)), ternary(sum(U.nOff)==0, '--', sprintf('%.3f', 100*sum(U.eOff)/max(1,sum(U.nOff)))));
end
fprintf('  units above the 0.5%% budget (one-sided binomial p):\n');
for k = 1:nC
    U = S(k).errUnits;  nT = nnz(U.nOn>0) + nnz(U.nOff>0);
    for g = 1:S(k).nG
        for s = 1:2
            n = [U.nOn(g) U.nOff(g)];  e = [U.eOn(g) U.eOff(g)];  d = {'ON','OFF'};
            if n(s) > 0 && e(s)/n(s) > 0.005
                p = 1 - binocdf(e(s)-1, n(s), 0.005);
                fprintf('    %-26s unit %2d %-3s %5.3f %% (p = %.4f, Bonferroni %.4f)\n', ...
                    S(k).name, g, d{s}, 100*e(s)/n(s), p, 0.05/nT);
            end
        end
    end
end

fprintf('\n=== Table 11: wrong fixings against solver performance\n');
fprintf('%-26s %5s %8s %6s %6s %9s %9s %15s\n', 'Case', 'eps %', 'days w/e', ...
    'rho_t', 'rho_c', 'cost/err', 'cost w/o', 'slow w/o err');
% One row per campaign: case index | eps [%] | bands | results
% (an empty pair means the eps = 0.5 % campaign already loaded in S)
rowsW = { 1, 0.5, '', ''
          1, 1.0, 'Fixing_bands_eps01.mat', 'Assess_NN_S2_thermal_eps01.mat'
          1, 3.0, 'Fixing_bands_eps03.mat', 'Assess_NN_S2_thermal_eps03.mat'
          3, 0.5, '', ''
          2, 0.5, '', ''
          4, 0.5, '', '' };
for r = 1:size(rowsW,1)
    k = rowsW{r,1};  s = S(k);
    if isempty(rowsW{r,3})
        sc = s.s2;
    else
        Bd = load(fullfile(M, C{k,2}, 'data', rowsW{r,3}), 'hiBand', 'loBand');
        s.hiBand = Bd.hiBand;  s.loBand = Bd.loBand;
        sc = load(fullfile(M, C{k,2}, 'results', rowsW{r,4}));
    end
    w = wrongPerDay(s, TU, TD, u0);
    ratio = s.base.COMP_TIME ./ sc.COMP_TIME;
    cost = 100*(sc.Z ./ s.base.Z - 1);
    use = s.closed & okOf(sc, numel(s.ok));
    cl = use & sc.OPT_GAP(:) <= s.epgap*1.001;
    c0 = cl & w == 0;                       % days closed by both and free of errors
    lm = [ones(nnz(cl),1), w(cl)] \ cost(cl);
    slow = use & ratio < 1;
    fprintf('%-26s %5.1f %7.1f%% %6.2f %6.2f %9.3f %9.3f %7d/%-4d (%3.0f%%)\n', s.name, rowsW{r,2}, ...
        100*mean(w(use)>0), corr(w(use), log(ratio(use)), 'Type','Spearman'), ...
        corr(w(cl), cost(cl), 'Type','Spearman'), lm(2), ...
        100*(sum(sc.Z(c0))/sum(s.base.Z(c0)) - 1), ...
        nnz(slow & w==0), nnz(slow), 100*nnz(slow & w==0)/max(1,nnz(slow)));
end

%% ---------------- Table 12: solver statistics --------------------------
f = fullfile(M, 'solver_metrics', 'results', 'solver_metrics.mat');
if isfile(f)
    fprintf('\n=== Table 12: solver statistics (medians over the days closed by both)\n');
    SM = load(f, 'SM');  SM = SM.SM;
    key = {'thermal_simpl','Thermal, simplified'; 'thermal_ext','Thermal, extended'
           'hydro_simpl','Hydrothermal, simplified'; 'hydro_ext','Hydrothermal, extended'};
    % Upper block: size of the problem and of the search
    fprintf('%-26s %6s %13s %15s %13s %11s %13s\n', 'Case', 'days', 'nodes', ...
        'iter [10^3]', 'binaries', 'cuts', 'root gap %');
    for r = 1:size(key,1)
        if ~isfield(SM, key{r,1}), continue; end
        R2 = SM.(key{r,1});  ib = R2.pairBase;  is = R2.pairS2;
        fprintf('%-26s %6d %6.0f %6.0f %7.1f %7.1f %6.0f %6.0f %5.0f %5.0f %6.1f %6.1f\n', key{r,2}, numel(ib), ...
            median(R2.base.nodes(ib)), median(R2.S2.nodes(is)), ...
            median(R2.base.iterations(ib))/1000, median(R2.S2.iterations(is))/1000, ...
            median(R2.base.binsPresolve(ib)), median(R2.S2.binsPresolve(is)), ...
            median(R2.base.cuts(ib)), median(R2.S2.cuts(is)), ...
            median(R2.base.rootGap(ib),'omitnan'), median(R2.S2.rootGap(is),'omitnan'));
    end
    % Lower block: where the time goes. Set-up is the part of the elapsed time
    % that the solver's own accounting does not cover (model generation and
    % writing of the solution), computed per day before taking the median.
    fprintf('\n%-26s %12s %12s %14s %14s %13s %15s\n', 'Case', 'set-up [s]', 'root [s]', ...
        'b&c [s]', 'total [s]', 'incumbents', 'solved at root %');
    for r = 1:size(key,1)
        if ~isfield(SM, key{r,1}), continue; end
        R2 = SM.(key{r,1});  ib = R2.pairBase;  is = R2.pairS2;
        sb = R2.base.csec(ib) - R2.base.rootTime(ib) - R2.base.bcTime(ib);
        ss = R2.S2.csec(is)   - R2.S2.rootTime(is)   - R2.S2.bcTime(is);
        fprintf('%-26s %6.1f %6.1f %6.1f %6.1f %7.1f %7.1f %7.1f %7.1f %5.0f %5.0f %7.1f %7.1f\n', key{r,2}, ...
            median(sb), median(ss), ...
            median(R2.base.rootTime(ib)), median(R2.S2.rootTime(is)), ...
            median(R2.base.bcTime(ib)), median(R2.S2.bcTime(is)), ...
            median(R2.base.csec(ib)), median(R2.S2.csec(is)), ...
            median(R2.base.incumbents(ib)), median(R2.S2.incumbents(is)), ...
            100*mean(R2.base.nodes(ib)==0), 100*mean(R2.S2.nodes(is)==0));
    end
else
    fprintf('\n(Table 12 needs matlab/solver_metrics/results/solver_metrics.mat)\n');
end

%% ---------------- Annex A: per-unit tables -----------------------------
fprintf('\n=== Tables 13 and 14 (Annex A): per-unit statistics\n');
pairs = {[1 2], 'Table 13, thermal cases'; [3 4], 'Table 14, hydrothermal cases'};
for q = 1:2
    kk = pairs{q,1};
    fprintf('\n%s (simplified | extended)\n', pairs{q,2});
    fprintf('%5s %7s %7s %6s %6s %6s   %7s %7s %6s %6s %6s\n', 'unit', 'ON %','trans','const','thr','NN', 'ON %','trans','const','thr','NN');
    nG = S(kk(1)).nG;
    nm = compose('%d',(1:16)');
    if nG == 28, nm = [nm; compose('T%d',(1:6)'); compose('P%d',(1:6)')]; end
    for g = 1:nG
        a = S(kk(1)).unitTab(g,:);  b = S(kk(2)).unitTab(g,:);
        fprintf('%5s %7.1f %7.2f %6.3f %6.3f %6.3f   %7.1f %7.2f %6.3f %6.3f %6.3f\n', nm{g}, ...
            100*a(1), a(2), a(4), a(5), a(3), 100*b(1), b(2), b(4), b(5), b(3));
    end
end

%% ---------------- Figures ---------------------------------------------
fprintf('\nEscribiendo las figuras en %s ...\n', figDir);
makeFigures(S, M, TU, TD, u0, figDir);
fprintf('Hecho.\n');

% ======================================================================
function X = canon3(X, G)
for k = 1:numel(G)
    idx = G{k};  cnt = sum(X(:,idx,:), 2);
    for j = 1:numel(idx), X(:,idx(j),:) = cnt >= j; end
end
end

function X = canon2(X, G)
for k = 1:numel(G)
    idx = G{k};  cnt = sum(X(idx,:), 1);
    for j = 1:numel(idx), X(idx(j),:) = cnt >= j; end
end
end

function m = metrics(Yh, Y)
tp = nnz(Yh & Y); tn = nnz(~Yh & ~Y); fp = nnz(Yh & ~Y); fn = nnz(~Yh & Y);
m.acc = (tp+tn)/numel(Y);  m.prec = tp/max(1,tp+fp);  m.rec = tp/max(1,tp+fn);
m.f1 = 2*m.prec*m.rec/max(eps, m.prec+m.rec);
end

function [ON, OFF] = applyBands(P, hi, lo, G, nG)
ON = P > hi(:);  OFF = P < lo(:);
for k = 1:numel(G)
    idx = G{k};
    for j = 2:numel(idx)
        ON(idx(1:j-1),:) = ON(idx(1:j-1),:) | ON(idx(j),:);
        OFF(idx(j:end),:) = OFF(idx(j:end),:) | OFF(idx(j-1),:);
    end
end
cf = ON & OFF;  ON(cf) = false;  OFF(cf) = false;
if nG == 28
    gt = 17:22;  gp = 23:28;
    for t = 1:size(P,2)
        nT = sum(ON(gt,t));  nP = sum(ON(gp,t));
        while nT + nP > 6
            if nP == 0 || (nT > 0 && P(gt(nT),t) <= P(gp(nP),t)), ON(gt(nT),t) = false; nT = nT-1;
            else, ON(gp(nP),t) = false; nP = nP-1; end
        end
    end
end
end

function [cov, err, U] = bandsOnPeriod(s, hi, lo, TU, TD, u0)
U = struct('nOn',zeros(s.nG,1),'eOn',zeros(s.nG,1),'nOff',zeros(s.nG,1),'eOff',zeros(s.nG,1));
nFix = 0; nTot = 0; nE = 0; nC = 0;
for d = find(s.ok)'
    P = squeeze(s.P(d,:,:));
    [on, off] = applyBands(P, hi, lo, s.G, s.nG);
    if s.ext
        if s.nG == 16, [on, off] = tr_filterFixings(on, off, TU, TD, u0);
        else,          [on, off] = phTR_filterFixings(on, off, TU, TD, u0); end
    end
    nFix = nFix + nnz(on) + nnz(off);  nTot = nTot + s.nG*24;
    if s.closed(d)
        tr = squeeze(s.Y(d,:,:));
        U.nOn = U.nOn + sum(on,2);    U.eOn = U.eOn + sum(on & tr==0, 2);
        U.nOff = U.nOff + sum(off,2); U.eOff = U.eOff + sum(off & tr==1, 2);
        nC = nC + nnz(on) + nnz(off);  nE = nE + nnz(on & tr==0) + nnz(off & tr==1);
    end
end
cov = nFix/nTot;  err = nE/max(1,nC);
end

function w = wrongPerDay(s, TU, TD, u0)
w = nan(size(s.ok));
for d = find(s.ok)'
    P = squeeze(s.P(d,:,:));
    [on, off] = applyBands(P, s.hiBand, s.loBand, s.G, s.nG);
    if s.ext
        if s.nG == 16, [on, off] = tr_filterFixings(on, off, TU, TD, u0);
        else,          [on, off] = phTR_filterFixings(on, off, TU, TD, u0); end
    end
    tr = squeeze(s.Y(d,:,:));
    w(d) = nnz(on & tr==0) + nnz(off & tr==1);
end
end

function a = speedStats(b, sc, sel, epg)
tb = b.COMP_TIME(sel);  ts = sc.COMP_TIME(sel);
gb = b.OPT_GAP(sel) > epg*1.001;  gs = sc.OPT_GAP(sel) > epg*1.001;
cl = ~gb & ~gs;  zb = b.Z(sel);  zs = sc.Z(sel);
a.meanS = mean(ts);  a.medS = median(ts);  a.agg = sum(tb)/sum(ts);
a.medRatio = median(tb./ts);  a.faster = 100*mean(ts < tb);
a.truncS = nnz(gs);  a.costClosed = 100*(sum(zs(cl))/sum(zb(cl)) - 1);
a.totB = sum(tb)/3600;  a.totS = sum(ts)/3600;   % total solver time [h]
end

function printSpeed(label, b, sc, s, M, caseDir)
ok = s.ok;
if isfield(sc, 'OK'), ok = ok & logical(sc.OK(:)); end
if isfield(sc, 'TREL') && contains(label, '3')
    sc.COMP_TIME = sc.COMP_TIME(:) + sc.TREL(:);
end
a = speedStats(b, sc, ok, s.epgap);
fx = NaN;  nTot = s.nG*24;
if isfield(sc,'FIXED_VAR'), fx = 100*mean(sc.FIXED_VAR(ok))/nTot;
elseif isfield(sc,'FIX_TH'), fx = 100*mean(sc.FIX_TH(ok)+sc.FIX_HY(ok))/nTot;
elseif isfield(sc,'FIXED_GEN'), fx = 100*mean(sc.FIXED_GEN(ok)+sc.FIXED_HGEN(ok))/nTot;
end
fs = '--';  if ~isnan(fx), fs = sprintf('%.1f', fx); end
fprintf('%-26s %-12s %6s %8.2f %8.2f %8.1f %7.2fx %7.2fx %6.1f%% %6d %+8.3f\n', '', label, fs, ...
    a.meanS, a.medS, a.totS, a.agg, a.medRatio, a.faster, a.truncS, a.costClosed);
end

function [U, cons, thr] = trivialRules(s, gT)
G = s.G;  nG = s.nG;
Ytr = round(s.lab.U);
if nG == 28 && isfield(s.lab,'Ug'), Ytr = cat(2, Ytr, round(s.lab.Ug), round(s.lab.Up)); end
keep = true(size(Ytr,1),1);
if isfield(s.lab,'OK'), keep = logical(s.lab.OK(:)); end
Ytr = canon3(Ytr(keep,:,:), G);  NDtr = s.lab.ND(keep,:);
Yev = s.Y(s.ok,:,:);  NDev = s.base.ND(s.ok,:);  Yh = s.Yhat(s.ok,:,:);
U = zeros(nG, 5);
cons = false(size(Yev));  thr = false(size(Yev));
for g = 1:nG
    ytr = squeeze(Ytr(:,g,:));  yev = squeeze(Yev(:,g,:));  yh = squeeze(Yh(:,g,:));
    isOn = mean(ytr(:)) >= 0.5;
    cons(:,g,:) = isOn;
    x = NDtr(:);  y = ytr(:);  cand = unique(quantile(x, 0:0.005:1));
    best = -1;  th = NaN;  sgn = 1;
    for sg = [1 -1]
        for q = cand(:)'
            a = mean((sg*(x - q) > 0) == y);
            if a > best, best = a; th = q; sgn = sg; end
        end
    end
    thr(:,g,:) = sgn*(NDev - th) > 0;
    U(g,:) = [mean(yev(:)), mean(sum(abs(diff(yev,1,2)),2)), mean(yh(:)==yev(:)), ...
              mean(yev(:)==isOn), mean((sgn*(NDev(:)-th) > 0) == yev(:))];
end
thr = canon3(thr, G);
end

function ok = okOf(s, n)
ok = true(n,1);
if isfield(s,'OK'), ok = logical(s.OK(:)); end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function makeFigures(S, M, TU, TD, u0, figDir)
cTh = [42 120 214]/255;  cHy = [235 104 52]/255;
ink = [11 11 11]/255;  ink2 = [82 81 78]/255;  gr = [0.86 0.86 0.84];
col = {cTh, cTh, cHy, cHy};  fs = 7;
order = [1 3 2 4];   % simplified row, extended row
f = figure('Units','centimeters','Position',[2 2 8.4 8.8],'Color','w','Visible','off');
tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for q = 1:4
    k = order(q);  ax = nexttile(tl);  hold(ax,'on');
    both = S(k).ok;  if isfield(S(k).s2,'OK'), both = both & logical(S(k).s2.OK(:)); end
    xy = [S(k).base.COMP_TIME(both), S(k).s2.COMP_TIME(both)];
    lim = S(k).tilim;  lo = 0.5;  hi = lim*1.6;
    plot(ax,[lo hi],[lo hi],'-','Color',ink2,'LineWidth',0.75);
    plot(ax,[lim lim],[lo hi],':','Color',ink2,'LineWidth',0.75);
    plot(ax,[lo hi],[lim lim],':','Color',ink2,'LineWidth',0.75);
    scatter(ax, xy(:,1), xy(:,2), 5, col{k}, 'filled', 'MarkerFaceAlpha',0.45, 'MarkerEdgeColor','none');
    tk = 10.^(0:ceil(log10(lim)));  tk = tk(tk <= hi);
    set(ax,'XScale','log','YScale','log','XLim',[lo hi],'YLim',[lo hi],'FontSize',fs, ...
        'XColor',ink2,'YColor',ink2,'GridColor',gr,'GridAlpha',1,'Box','off','LineWidth',0.5, ...
        'XTick',tk,'YTick',tk,'XTickLabel',string(tk),'YTickLabel',string(tk),'TickLength',[0.015 0.015]);
    grid(ax,'on'); ax.XMinorGrid='off'; ax.YMinorGrid='off'; axis(ax,'square');
    title(ax, S(k).name, 'FontSize', fs+0.5, 'FontWeight','normal', 'Color', ink);
    a = speedStats(S(k).base, S(k).s2, both, S(k).epgap);
    text(ax, lo*1.4, lim/2.2, sprintf('%.2f\\times aggregate\n%.2f\\times median', a.agg, a.medRatio), ...
        'FontSize', fs-0.5, 'Color', ink, 'VerticalAlignment','top');
    text(ax, lim/1.12, lo*1.5, sprintf('limit %d s', lim), 'FontSize', fs-1, 'Color', ink2, ...
        'HorizontalAlignment','right','VerticalAlignment','bottom');
end
xlabel(tl,'Baseline MIP time per day [s]','FontSize',fs,'Color',ink);
ylabel(tl,'Strategy 2 time per day [s]','FontSize',fs,'Color',ink);
exportgraphics(f, fullfile(figDir,'S2_time_per_day.pdf'), 'ContentType','vector');
close(f);

% budget sweep
eps_ = [0.002 0.005 0.01 0.02 0.03 0.05];
cov = zeros(4,numel(eps_));  err = zeros(4,numel(eps_));
for k = 1:4
    for e = 1:numel(eps_)
        [hb, lb] = calibrate(S(k), eps_(e));
        [cov(k,e), err(k,e)] = bandsOnPeriod(S(k), hb, lb, TU, TD, u0);
    end
end
f = figure('Units','centimeters','Position',[2 2 8.4 9.6],'Color','w','Visible','off');
tl = tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');
style = {'--','-','--','-'};  mk = {'o','s','o','s'};  fill_ = [false true false true];
short = {'Thermal, simpl.','Thermal, ext.','Hydrothermal, simpl.','Hydrothermal, ext.'};
x = 100*eps_;
ax1 = nexttile(tl); hold(ax1,'on');
ax2 = nexttile(tl); hold(ax2,'on');
plot(ax2,[0.15 7],[0.15 7],'-','Color',ink2,'LineWidth',0.75);
text(ax2, 3.2, 4.6, 'budget','FontSize',fs-1,'Color',ink2,'Rotation',33,'HorizontalAlignment','center');
h2 = gobjects(4,1);
for k = 1:4
    mf = 'w';  if fill_(k), mf = col{k}; end
    plot(ax1, x, 100*cov(k,:), style{k}, 'Color', col{k}, 'LineWidth',1.25, 'Marker', mk{k}, ...
        'MarkerSize',4,'MarkerFaceColor',mf,'MarkerEdgeColor',col{k});
    plot(ax2, x, 100*err(k,:), style{k}, 'Color', col{k}, 'LineWidth',1.25, 'Marker', mk{k}, ...
        'MarkerSize',4,'MarkerFaceColor',mf,'MarkerEdgeColor',col{k});
    h2(k) = plot(ax2, NaN, NaN, style{k}, 'Color', col{k}, 'LineWidth',1.25, 'Marker', mk{k}, ...
        'MarkerSize',4,'MarkerFaceColor',mf,'MarkerEdgeColor',col{k});
end
for ax = [ax1 ax2]
    set(ax,'XScale','log','XLim',[0.15 7],'XTick',[0.2 0.5 1 2 5],'XTickLabel',{'0.2','0.5','1','2','5'}, ...
        'FontSize',fs,'XColor',ink2,'YColor',ink2,'GridColor',gr,'GridAlpha',1,'Box','off', ...
        'LineWidth',0.5,'TickLength',[0.01 0.01]);
    grid(ax,'on'); ax.XMinorGrid='off'; ax.YMinorGrid='off';
    xline(ax, 0.5, ':', 'Color', ink2, 'LineWidth', 0.75);
end
set(ax1,'YLim',[15 70]);  ylabel(ax1,'Fixed binary variables [%]','FontSize',fs,'Color',ink);
text(ax1, 0.53, 67, 'operating point','FontSize',fs-1,'Color',ink2,'VerticalAlignment','top');
set(ax2,'YScale','log','YLim',[0.005 7],'YTick',[0.01 0.1 1],'YTickLabel',{'0.01','0.1','1'});
ylabel(ax2,'Realized error [%]','FontSize',fs,'Color',ink);
xlabel(ax2,'Error budget \epsilon [%]','FontSize',fs,'Color',ink);
lg = legend(ax2, h2, short, 'Location','southeast','FontSize',fs-0.5,'Box','off','TextColor',ink);
lg.ItemTokenSize = [14 8];
exportgraphics(f, fullfile(figDir,'Budget_sweep.pdf'), 'ContentType','vector');
close(f);
end

function [hi, lo] = calibrate(s, e)
% Recalibrates the bands at a given budget, exactly as calibrate_fixing does.
cand = [0.50:0.01:0.99, 0.995, 0.999];
nG = s.nG;  hyd = [];  if nG == 28, hyd = 17:28; end
Y = round(s.lab.U);
if nG == 28 && isfield(s.lab,'Ug'), Y = cat(2, Y, round(s.lab.Ug), round(s.lab.Up)); end
keep = true(size(Y,1),1);  if isfield(s.lab,'OK'), keep = logical(s.lab.OK(:)); end
Y = canon3(Y(keep,:,:), s.G);  ND = s.lab.ND(keep,:);
rng(0);
n = size(ND,1);
cv1 = cvpartition(n,'HoldOut',0.40);  idxHold = find(test(cv1));
cv2 = cvpartition(numel(idxHold),'HoldOut',0.75);  idxRest = idxHold(test(cv2));
cv3 = cvpartition(numel(idxRest),'HoldOut',0.5);
idxCal = idxRest(training(cv3));
P = zeros(24, nG, numel(idxCal));
for i = 1:numel(idxCal), P(:,:,i) = predict(s.net, ND(idxCal(i),:)'); end
Yc = permute(Y(idxCal,:,:), [3 2 1]);
eON = e*ones(1,nG);  eOFF = zeros(1,nG);  eOFF(hyd) = e;
hi = nan(1,nG);  lo = nan(1,nG);
for g = 1:nG
    p = squeeze(P(:,g,:));  y = squeeze(double(Yc(:,g,:)));
    for h = cand
        sel = p > h;  k = nnz(y(sel)==0);  m = nnz(sel);
        u = 1;  if m > 0 && k < m, u = betainv(0.95, k+1, m-k); end
        if u <= eON(g), hi(g) = h; break; end
    end
    for l = sort(1-cand,'descend')
        sel = p < l;  k = nnz(y(sel)==1);  m = nnz(sel);
        u = 1;  if m > 0 && k < m, u = betainv(0.95, k+1, m-k); end
        if u <= eOFF(g), lo(g) = l; break; end
    end
end
if isfield(s,'disableON'), hi(s.disableON) = NaN; end
end
