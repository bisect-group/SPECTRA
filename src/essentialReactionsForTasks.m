function [core,perTask] = essentialReactionsForTasks(model,tasks,candidates,printLevel)
% Find the reactions essential to any of a list of metabolic tasks.
%
% USAGE:
%
%    [core,perTask] = essentialReactionsForTasks(model,tasks,candidates,printLevel)
%
% INPUTS:
%    model:              COBRA model structure to search in, normally the
%                        universal model being extracted from.
%    tasks:              Task structures from metabolicTask, as a struct
%                        array. Those marked shouldFail are skipped, since
%                        nothing should be forced in to support them.
%
% OPTIONAL INPUTS:
%    candidates:         Reaction IDs to restrict the search to. (Default:
%                        all)
%    printLevel:         0 prints nothing, 1 prints a line per task and the
%                        size of the union. (Default: 1)
%
% OUTPUTS:
%    core:               Union of the essential reactions, as a cell array
%                        of char, sorted. Hand this to spectraME as its
%                        core set, via findRxnIDs(model, core).
%    perTask:            Struct array with fields taskID and essential,
%                        giving each task's own set.
%
% NOTE:
%    This is the step that turns "the model must perform these tasks" into a
%    core reaction set that spectraME can enforce. Requiring a model to
%    perform a task is difficult to express as a constraint, but requiring it
%    to contain a reaction is exactly what spectraME already does.
%
%    It is a necessary condition rather than a sufficient one. A task with
%    two alternative routes has no essential reactions at all, since either
%    route can be removed on its own, so forcing the essential set does not
%    guarantee the task survives. It carries most of the weight cheaply;
%    checkTasks on the extracted model is what tells you whether anything is
%    left to repair.
%
% EXAMPLE:
%
%    [core,perTask] = essentialReactionsForTasks(universal, tasks);
%    extracted = spectraME(universal, findRxnIDs(universal,core), 1e-4);
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('candidates', 'var')
    candidates = {};
end
if ~exist('printLevel', 'var') || isempty(printLevel)
    printLevel = 1;
end

keep = ~[tasks.shouldFail];
tasks = tasks(keep);
n = numel(tasks);

perTask = repmat(struct('taskID','','essential',{{}}), n, 1);
core = {};

for k = 1:n
    found = essentialReactionsForTask(model, tasks(k), candidates);
    perTask(k).taskID = tasks(k).id;
    perTask(k).essential = found;
    core = [core; found]; %#ok<AGROW>
    if printLevel > 0
        fprintf('%-24s %3d essential\n', tasks(k).id, numel(found));
    end
end

core = unique(core);
core = core(:);

if printLevel > 0
    fprintf('union (the core set): %d reactions essential to at least one of %d tasks\n', ...
        numel(core), n);
end
end
