function runRefitJob(cond, nWorkers)
% runRefitJob  Runs RefitRefolding.m for one condition in its own pool, for a
% standalone MATLAB process that does not depend on an interactive session:
%   matlab -batch "runRefitJob('Relax')" -logfile refit_Relax.log
% Utility (not an entry point); the checkpoint lets a rerun resume.
if nargin < 2, nWorkers = 4; end
cd(fileparts(mfilename('fullpath')));
pool = parpool('Processes', nWorkers);
try
    RefitRefolding;
catch err
    fprintf(2, '%s\n', getReport(err));
end
delete(pool);
end
