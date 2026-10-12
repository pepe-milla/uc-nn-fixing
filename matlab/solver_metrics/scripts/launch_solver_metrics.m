
% Launches the solver-metrics campaigns unattended, as K independent MATLAB
% processes that share the days between them (see run_solver_metrics).
%
% On the laptop (4 cores) K = 2 is the right choice: two processes with the
% two CPLEX threads of the i9 campaigns each. Each process opens its own
% window and writes its own log in results\logs. Closing this MATLAB does
% not stop them. To stop them, close their windows; relaunching resumes
% from the cache.
%
% Order: the two thermal cases first (about 4 h each with K = 2), then the
% 61-day samples of the two hydrothermal cases (about 6 and 7 h).

cases = {'thermal_simpl', 'thermal_ext', 'hydro_simpl', 'hydro_ext'};
K = 2;

matlabExe = fullfile(matlabroot, 'bin', 'matlab.exe');
scriptDir = fileparts(mfilename('fullpath'));
caseList  = ['{' strjoin(strcat('''', cases, ''''), ',') '}'];
for k = 1:K
    cmd = sprintf('"%s" -batch "cd(''%s''); run_solver_metrics(%s, %d, %d)" &', ...
        matlabExe, scriptDir, caseList, k, K);
    system(cmd);
    fprintf('Lanzada la instancia %d de %d\n', k, K);
end
fprintf('Sigue el progreso en %s\n', fullfile(fileparts(scriptDir), 'results', 'logs'));
