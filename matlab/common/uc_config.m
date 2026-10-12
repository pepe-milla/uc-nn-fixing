function c = uc_config()
%UC_CONFIG  Local paths of this repository. Edit only this file.
%
%   c = UC_CONFIG() returns a struct with the paths used by the scripts:
%
%     gamsExe     the GAMS executable. The campaigns of the paper were run
%                 with GAMS 23.5 and CPLEX 12.2, both commercial and not
%                 distributed here
%     workerRoot  folder in which the parallel drivers create one copy of
%                 the model per worker (worker1, worker2, ...). It should be
%                 on a fast local disk and can be deleted afterwards
%     repoRoot    root of this repository, derived from the location of this
%                 file; the remaining paths hang from it
%
%   Only gamsExe and workerRoot depend on the machine.

c.gamsExe    = 'C:\Program Files\GAMS23.5\gams.exe';
c.workerRoot = fullfile(tempdir, 'uc_workers');

c.repoRoot     = fileparts(fileparts(fileparts(mfilename('fullpath'))));
c.gamsRoot     = fullfile(c.repoRoot, 'gams');
c.dataRoot     = fullfile(c.repoRoot, 'data');
c.matlabRoot   = fullfile(c.repoRoot, 'matlab');
end
