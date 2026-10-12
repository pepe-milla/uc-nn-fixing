
% Summary of the solver-metrics campaigns (run_solver_metrics), and check
% against the i9 campaigns of the article.
%
% The check: on a day that both machines close within the tolerance, with
% the same CPLEX version, threads=2 and deterministic parallel mode, the
% search should be the same, and so should the objective value. Where it
% is, the metrics measured here describe the runs reported in the article.

scriptDir = fileparts(mfilename('fullpath'));
caseRoot  = fileparts(scriptDir);
M         = fileparts(caseRoot);
cases = { 'thermal_simpl', 0.001, fullfile(M,'thermal_base','results','MIP_Baseline_thermal_2021-2023.mat'), fullfile(M,'thermal_base','results','Assess_NN_S2_thermal_eps005.mat')
          'thermal_ext',   0.001, fullfile(M,'thermal_extended','results','MIP_Baseline_2021-2023_i9.mat'),  fullfile(M,'thermal_extended','results','Assess_NN_S2_TRamp_i9.mat')
          'hydro_simpl',   0.015, fullfile(M,'hydrothermal_base','results','MIP_Baseline_2021-2023.mat'),                fullfile(M,'hydrothermal_base','results','Assess_NN_par.mat')
          'hydro_ext',     0.015, fullfile(M,'hydrothermal_extended','results','MIP_Baseline_2021-2023.mat'),        fullfile(M,'hydrothermal_extended','results','Assess_NN_S2_PHES_TRamp.mat') };
fld = {'nodes','iterations','cuts','binsPresolve','rootTime','bcTime','incumbents','rootGap','maxLeft'};
SM = struct();

for c = 1:size(cases,1)
    name = cases{c,1};  epg = cases{c,2};
    cdir = fullfile(caseRoot, 'results', 'cache', name);
    if ~isfolder(cdir), continue; end
    I9 = {load(cases{c,3}), load(cases{c,4})};
    R = struct();
    for m = 1:2
        md = {'base','S2'};  md = md{m};
        f = dir(fullfile(cdir, [md '_*.mat']));
        n = numel(f);
        R.(md).day = zeros(n,1);  R.(md).z = nan(n,1);  R.(md).gap = nan(n,1);  R.(md).csec = nan(n,1);
        R.(md).nFix = nan(n,1);   R.(md).infeas = false(n,1);
        for q = fld, R.(md).(q{1}) = nan(n,1); end
        for i = 1:n
            S = load(fullfile(f(i).folder, f(i).name));
            R.(md).day(i) = S.d;  R.(md).z(i) = S.z;  R.(md).gap(i) = S.gap;  R.(md).csec(i) = S.csec;
            R.(md).nFix(i) = S.nFix;  R.(md).infeas(i) = S.infeas;
            for q = fld, R.(md).(q{1})(i) = S.sm.(q{1}); end
        end
        % check against the i9
        zi = I9{m}.Z(R.(md).day);  gi = I9{m}.OPT_GAP(R.(md).day);
        closed = R.(md).gap <= epg*1.001 & gi(:) <= epg*1.001 & ~R.(md).infeas;
        same = closed & abs(R.(md).z - zi(:)) <= 1e-6*abs(zi(:));
        R.(md).sameAsI9 = same;  R.(md).bothClosed = closed;
        fprintf('%-13s %-4s: %4d dias | cerrados en ambas maquinas %4d | mismo objetivo que el i9 %4d (%.1f %%)\n', ...
            name, md, n, nnz(closed), nnz(same), 100*nnz(same)/max(1,nnz(closed)));
    end

    % paired days, both modes closed on this machine
    [dd, ib, is] = intersect(R.base.day, R.S2.day);
    ok = R.base.gap(ib) <= epg*1.001 & R.S2.gap(is) <= epg*1.001 & ~R.base.infeas(ib) & ~R.S2.infeas(is);
    ib = ib(ok);  is = is(ok);
    fprintf('  dias emparejados y cerrados en ambos modos: %d\n', numel(ib));
    fprintf('  %-14s %12s %12s %8s\n', 'metrica (mediana)', 'baseline', 'S2', 'cociente');
    for q = fld
        b = R.base.(q{1})(ib);  s = R.S2.(q{1})(is);
        fprintf('  %-14s %12.1f %12.1f %8.2f\n', q{1}, median(b,'omitnan'), median(s,'omitnan'), median(b,'omitnan')/median(s,'omitnan'));
    end
    fprintf('  resueltos en la raiz (0 nodos): baseline %.1f %% | S2 %.1f %%\n', ...
        100*mean(R.base.nodes(ib) == 0), 100*mean(R.S2.nodes(is) == 0));
    fprintf('  cuota del tiempo en la raiz: baseline %.1f %% | S2 %.1f %%\n', ...
        100*median(R.base.rootTime(ib)./(R.base.rootTime(ib)+R.base.bcTime(ib)),'omitnan'), ...
        100*median(R.S2.rootTime(is)./(R.S2.rootTime(is)+R.S2.bcTime(is)),'omitnan'));
    tr = R.base.csec(ib)./R.S2.csec(is);  nr = (R.base.nodes(ib)+1)./(R.S2.nodes(is)+1);
    fprintf('  Spearman(cociente de tiempos, cociente de nodos) = %.2f\n', corr(log(tr), log(nr), 'Type','Spearman'));
    R.pairBase = ib;  R.pairS2 = is;  R.days = dd(ok);
    SM.(name) = R;
end
save(fullfile(caseRoot, 'results', 'solver_metrics.mat'), 'SM');
fprintf('\nGuardado en %s\n', fullfile(caseRoot, 'results', 'solver_metrics.mat'));
