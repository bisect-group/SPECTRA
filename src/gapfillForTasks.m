function result = gapfillForTasks(draft,universal,tasks,tol,probType,freeRxns,solveTime,seed,essentialRxns,coreSolve,repair)
% Adds reactions until a draft can perform a list of metabolic tasks
%
% USAGE:
%
%    result = gapfillForTasks(draft,universal,tasks,tol,probType,freeRxns,solveTime,seed,essentialRxns,coreSolve,repair)
%
% INPUTS:
%    draft:             COBRA model structure that cannot yet perform every
%                       task.
%    universal:         COBRA model structure to borrow reactions from.
%    tasks:             Task structures from metabolicTask, as a struct
%                       array.
%
% OPTIONAL INPUTS:
%    tol:               Flux threshold handed to spectraME (Default: 1e-4).
%    probType:          Formulation for the core solve, any spectraME
%                       accepts (Default: 'minNetMILP').
%    freeRxns:          Further reaction IDs that cost nothing to include,
%                       beyond the draft's own. (Default: none)
%    solveTime:         Maximum runtime per MILP, in seconds (Default: 300).
%    seed:              Seed for the random coefficients. (Default: none)
%    essentialRxns:     Reaction IDs to use as the core set instead of
%                       computing them. Useful when
%                       essentialReactionsForTasks has already been run,
%                       since it is the expensive step. (Default: [])
%    coreSolve:         Whether to do the single extraction around the core
%                       set (Default: true).
%    repair:            Whether to follow up task by task on whatever still
%                       fails (Default: true).
%
% OUTPUTS:
%    result:            Structure with fields:
%
%                         * model - the gap-filled model
%                         * added - reaction IDs taken from the universal
%                           model, sorted
%                         * satisfied - IDs of the tasks the result
%                           performs
%                         * unsatisfied - IDs of the tasks it still does not
%                         * failed - struct mapping a task ID to why its
%                           solve failed
%                         * skipped - IDs of tasks the universal model
%                           cannot perform, which were left out entirely
%                         * summary - a readable multi-line summary, as char
%
% EXAMPLE:
%
%    result = gapfillForTasks(draft, universal, tasks);
%    result.unsatisfied   % the tasks still not performed
%
% NOTE:
%    Two steps, and the second is why this is not just an extraction. The
%    core solve forces in the reactions no task can do without, which is a
%    necessary condition but not a sufficient one: a task with two
%    alternative routes has no essential reactions at all, so forcing the
%    essential set does not guarantee it survives. The repair pass then
%    goes task by task over whatever still fails, inside each task's own
%    constraints, minimising the count of reactions not already present.
%
%    Repair always uses minNetMILP whatever probType the core solve used.
%    Repairing one task is a small problem and counting reactions exactly
%    is what that formulation is for.
%
%    Tasks the universal model itself cannot perform are skipped rather
%    than attempted. There is nothing to add that would help, so they would
%    sit in the results as permanent failures.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end
if ~exist('probType', 'var') || isempty(probType)
    probType = 'minNetMILP';
end
if ~exist('freeRxns', 'var') || isempty(freeRxns)
    freeRxns = {};
end
if ~exist('solveTime', 'var') || isempty(solveTime)
    solveTime = 300;
end
if ~exist('seed', 'var')
    seed = [];
end
if ~exist('essentialRxns', 'var')
    essentialRxns = {};
end
if ~exist('coreSolve', 'var') || isempty(coreSolve)
    coreSolve = true;
end
if ~exist('repair', 'var') || isempty(repair)
    repair = true;
end

draftIDs = draft.rxns(:);
draftIDs = draftIDs(ismember(draftIDs, universal.rxns));
freeIDs = unique([draftIDs; freeRxns(:)]);

% Only work on tasks the universal model can actually perform.
outcomes = checkTasks(universal, tasks, 0);
isUsable = [outcomes.ok]' & ~[tasks.shouldFail]';
usable = tasks(isUsable);
skipped = {tasks(~isUsable & ~[tasks.shouldFail]').id}';
if ~isempty(skipped)
    warning('gapfillForTasks:unusableTasks', ...
        ['%d task(s) the universal model cannot perform are skipped; ' ...
         'there is nothing to add that would help. The first is ''%s''.'], ...
        numel(skipped), skipped{1});
end

failedIDs = {};
failedWhy = {};
added = {};

% 1. the core solve: force in what no task can do without
if ~isempty(essentialRxns)
    core = essentialRxns(:);
elseif coreSolve
    core = essentialReactionsForTasks(universal, usable, {}, 0);
else
    core = {};
end

if ~isempty(core) && coreSolve
    weights = ones(numel(universal.rxns),1);
    weights(ismember(universal.rxns, freeIDs)) = 0;
    try
        filled = spectraME(universal, findRxnIDs(universal,core), tol, ...
            'stoichiometry', weights, 1, 'coreDirection', probType, ...
            solveTime, 0, {}, seed);
        added = unique([added; setdiff(filled.rxns, draftIDs)]);
    catch err
        failedIDs{end+1,1} = 'essentialReactions';
        failedWhy{end+1,1} = sprintf('%s: %s', err.identifier, err.message);
        warning('gapfillForTasks:coreSolveFailed', ...
            'the core solve failed (%s)', err.message);
    end
end

model = rebuildModel(universal, [draftIDs; added]);

% 2. repair whatever the core set could not express
if repair
    outcomes = checkTasks(model, usable, 0);
    for k = 1:numel(usable)
        if outcomes(k).ok
            continue
        end
        present = unique([draftIDs; added]);
        try
            constrained = taskConstraints(universal, usable(k));
        catch err
            failedIDs{end+1,1} = usable(k).id; %#ok<AGROW>
            failedWhy{end+1,1} = sprintf('%s: %s', err.identifier, err.message); %#ok<AGROW>
            continue
        end

        n = numel(constrained.rxns);
        weights = ones(n,1);
        weights(ismember(constrained.rxns, present)) = 0;
        % Inside the task's own constraints the model is already forced to
        % perform it, so minimising the count of reactions that are not
        % already present gives the cheapest repair directly.
        [reacInd,~,stat] = minNet(constrained, zeros(n,1), weights, tol, 1, ...
            'minNetMILP', solveTime);
        if isempty(reacInd) || (~isempty(stat) && stat ~= 1 && stat ~= 3)
            failedIDs{end+1,1} = usable(k).id; %#ok<AGROW>
            failedWhy{end+1,1} = sprintf('repair solve returned status %d', stat); %#ok<AGROW>
            continue
        end

        kept = constrained.rxns(reacInd);
        % Drop the task scaffolding, which is not part of the universal model.
        kept = kept(ismember(kept, universal.rxns));
        added = unique([added; setdiff(kept, draftIDs)]);
    end
    model = rebuildModel(universal, [draftIDs; added]);
end

% 3. report
outcomes = checkTasks(model, usable, 0);
isOK = [outcomes.ok]';
usableIDs = {usable.id}';

result.model       = model;
result.added       = sort(added);
result.satisfied   = usableIDs(isOK);
result.unsatisfied = usableIDs(~isOK);
result.skipped     = skipped;
if isempty(failedIDs)
    result.failed = struct();
else
    result.failed = cell2struct(failedWhy, matlab.lang.makeValidName(failedIDs), 1);
end

lines = {sprintf('%d reactions added to a %d-reaction draft', ...
    numel(result.added), numel(draftIDs))};
lines{end+1} = sprintf('  performs %d of %d usable tasks', sum(isOK), numel(usable));
if any(~isOK)
    still = usableIDs(~isOK);
    lines{end+1} = sprintf('  still fails: %s', strjoin(still(1:min(5,end))', ', '));
end
if ~isempty(skipped)
    lines{end+1} = sprintf('  %d task(s) skipped: the universal model cannot do them either', ...
        numel(skipped));
end
result.summary = strjoin(lines, newline);
end


function reduced = rebuildModel(universal,keepIDs)
% Returns the universal model cut down to the given reactions
%
% USAGE:
%
%    reduced = rebuildModel(universal,keepIDs)
%
% INPUTS:
%    universal:         COBRA model structure to cut down.
%    keepIDs:           Reaction IDs to retain.
%
% OUTPUTS:
%    reduced:           The reduced model.
%
% NOTE:
%    Orphaned metabolites stay: the draft may still need them, and they
%    cost nothing. COBRA's removeRxns drops them by default, hence metFlag
%    false.

drop = setdiff(universal.rxns, keepIDs);
if isempty(drop)
    reduced = universal;
else
    reduced = removeRxns(universal, drop, 'metFlag', false);
end
end
