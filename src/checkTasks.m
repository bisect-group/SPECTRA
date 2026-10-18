function results = checkTasks(model,tasks,printLevel)
% Check a model against a list of metabolic tasks.
%
% USAGE:
%
%    results = checkTasks(model,tasks,printLevel)
%
% INPUTS:
%    model:              COBRA model structure.
%    tasks:              Task structures from metabolicTask, as a struct
%                        array.
%
% OPTIONAL INPUTS:
%    printLevel:         0 prints nothing, 1 prints a line per task and a
%                        tally at the end. (Default: 1)
%
% OUTPUTS:
%    results:            Struct array, one entry per task, with the fields
%                        checkTask returns: taskID, feasible, ok and
%                        missing.
%
% NOTE:
%    Run this on the universal model before gap-filling towards tasks. A task
%    the universal model cannot perform can never be gap-filled towards,
%    because there is nothing to add that would help, and it will sit in the
%    results as a permanent failure.
%
% EXAMPLE:
%
%    results = checkTasks(model, tasks);
%    failed = {results(~[results.ok]).taskID};
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('printLevel', 'var') || isempty(printLevel)
    printLevel = 1;
end

n = numel(tasks);
results = repmat(struct('taskID','','feasible',false,'ok',false,'missing',{{}}), n, 1);

for k = 1:n
    results(k) = checkTask(model, tasks(k));
    if printLevel > 0
        if isempty(results(k).missing)
            note = '';
        else
            note = sprintf('  [%d metabolite(s) not in the model]', ...
                numel(results(k).missing));
        end
        fprintf('%-24s feasible=%d  as expected=%d%s\n', ...
            results(k).taskID, results(k).feasible, results(k).ok, note);
    end
end

if printLevel > 0 && n > 0
    fprintf('%d of %d tasks behaved as expected.\n', sum([results.ok]), n);
end
end
