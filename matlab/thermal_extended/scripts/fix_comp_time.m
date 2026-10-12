
% Repairs the COMP_TIME stored by the campaigns that were run before the log
% parser was fixed.
%
% The old parser took the last number of the GAMS line "Job ... Stop ...
% elapsed H:MM:SS.sss", which is the seconds field alone: a solve of 2 min
% 45 s was recorded as 45 s. Anything under a minute was recorded correctly,
% which is why it went unnoticed.
%
% The loss is always a whole number of minutes, and every file also stores T,
% the wall time measured with tic/toc around the GAMS call. T exceeds the
% true elapsed time by a constant process-launch overhead - measured at
% 0.239 s with a standard deviation of 0.032 s, uncorrelated with the length
% of the solve - so the true value can be recovered exactly:
%
%     COMP_TIME_true = COMP_TIME_raw + 60 * round((T - overhead - COMP_TIME_raw)/60)
%
% This is a reconstruction, not an estimate: on the three thermal campaigns
% the residual |T - overhead - reconstructed| never exceeds 0.087 s, and the
% distribution of T - COMP_TIME_raw is cleanly trimodal at 0.24, 60.24 and
% 120.24 s with nothing in between. The script verifies both properties on
% each file and refuses to touch anything that does not satisfy them.
%
% The original column is kept as COMP_TIME_raw and a copy of the untouched
% file is written alongside, so nothing is lost. Running it twice is safe:
% a file that already carries COMP_TIME_raw is skipped.

files = {
    fullfile('..','results','MIP_Results_2019-2021.mat')
    fullfile('..','results','MIP_Results_2021-2023.mat')
    fullfile('..','results','Assess_NN_S2_TRamp.mat')
    % (the pilot file of the hydrothermal extended case is not distributed)
    };

scriptDir = fileparts(mfilename('fullpath'));
for k = 1:numel(files)
    f = fullfile(scriptDir, files{k});
    if ~isfile(f)
        fprintf('  (no existe, se omite)  %s\n', files{k});
        continue
    end
    repairFile(f);
end

% ---- Local functions ----
function repairFile(f)
S = load(f);
[~, name] = fileparts(f);

if isfield(S, 'COMP_TIME_raw')
    fprintf('  ya corregido, se omite   %s\n', name);
    return
end
if ~isfield(S, 'COMP_TIME') || ~isfield(S, 'T')
    fprintf('  sin COMP_TIME o sin T    %s\n', name);
    return
end

CT = S.COMP_TIME(:);  T = S.T(:);
m  = isfinite(CT) & isfinite(T);
if nnz(m) == 0
    fprintf('  sin datos utilizables    %s\n', name);
    return
end

d = T(m) - CT(m);
low = d < 30;                       % days the old parser got right
if nnz(low) < 0.5*nnz(m)
    fprintf(2, '  ABORTA: menos de la mitad de los dias tienen desfase normal  %s\n', name);
    return
end
ov = median(d(low));

% The gap must sit at multiples of 60 s plus the constant overhead, with
% nothing in between; otherwise the wrap hypothesis does not hold here.
if max(abs(d - ov - 60*round((d - ov)/60))) > 1
    fprintf(2, '  ABORTA: el desfase no cae en multiplos de 60 s  %s\n', name);
    return
end

kmin = round((T(m) - ov - CT(m))/60);
if all(kmin == 0)
    fprintf('  sin minutos perdidos     %s (nada que corregir)\n', name);
    return
end

CTfix = CT;
CTfix(m) = CT(m) + 60*kmin;
res = max(abs(T(m) - ov - CTfix(m)));
if res > 1
    fprintf(2, '  ABORTA: residuo de %.3f s tras reconstruir  %s\n', res, name);
    return
end

% Back up the untouched file, then write
bak = strrep(f, '.mat', '_rawtime_backup.mat');
if ~isfile(bak), copyfile(f, bak); end

S.COMP_TIME_raw = S.COMP_TIME;
S.COMP_TIME = reshape(CTfix, size(S.COMP_TIME));
S.COMP_TIME_overhead = ov;
S.COMP_TIME_note = ['COMP_TIME reconstruido: el parser antiguo perdia los minutos. ' ...
                    'COMP_TIME_raw conserva el valor original. Ver fix_comp_time.m'];
save(f, '-struct', 'S');

fprintf('  CORREGIDO  %-28s  %d dias con minutos perdidos (k=1: %d, k=2: %d, k>=3: %d)\n', ...
    name, nnz(kmin > 0), nnz(kmin == 1), nnz(kmin == 2), nnz(kmin >= 3));
fprintf('             media %6.2f s -> %6.2f s   max %7.2f s -> %7.2f s   residuo max %.3f s\n', ...
    mean(CT(m)), mean(CTfix(m)), max(CT(m)), max(CTfix(m)), res);
end
