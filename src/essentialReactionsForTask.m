function essential = essentialReactionsForTask(model,task,candidates)
% Find the reactions without which a metabolic task becomes infeasible.
%
% USAGE:
%   essential = essentialReactionsForTask(model,task,candidates)
%
% INPUTS:
%   model: COBRA model structure to search in.
%   task:  Task structure from metabolicTask. A task the model cannot
%          perform at all has no essential reactions.
%
% OPTIONAL INPUT:
%   candidates: Reaction IDs to restrict the search to, as a cell array of
%               char. (Default: every reaction in the model)
%
% OUTPUT:
%   essential: Reaction IDs without which the task is infeasible, as a cell
%              array of char, sorted.
%
% NOTE:
%   Only a reaction carrying flux in a feasible solution can be essential:
%   if some solution leaves it at zero, fixing it to zero keeps that
%   solution available. The search therefore starts from the reactions one
%   solution uses and shrinks that set by repeatedly minimising the total
%   flux through whatever is left, which rules out in a single LP every
%   reaction some other route can avoid. Only the survivors are then tested
%   one at a time, which is the expensive part.
%
%   The reactions this function adds to apply the task are never reported;
%   they are scaffolding, not part of the model.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

FLUX_EPSILON = 1e-9;

essential = {};

if ~isempty(missingMetabolites(model,task))
    return
end

[constrained,addedRxns] = taskConstraints(model,task);
constrained.c(:) = 0;

if ~exist('candidates', 'var') || isempty(candidates)
    pool = constrained.rxns;
else
    if ischar(candidates)
        candidates = {candidates};
    end
    pool = candidates(:);
    pool = pool(ismember(pool, constrained.rxns));
end
% Never report the scaffolding as essential.
pool = pool(~ismember(pool, addedRxns));

% A task the model cannot perform has no essential reactions: removing a
% reaction cannot make an already infeasible problem worse.
solution = optimizeCbModel(constrained);
if solution.stat ~= 1
    return
end

% Shrink the candidates. Anything a minimal-flux solution can leave at zero
% is not essential, and one LP rules out many at once.
previous = -1;
while previous ~= numel(pool) && ~isempty(pool)
    previous = numel(pool);
    solution = optimizeCbModel(constrained, 'max', 'one');
    if solution.stat ~= 1 || ~isfield(solution,'x') || isempty(solution.x)
        break
    end
    carried = abs(solution.x(findRxnIDs(constrained, pool))) > FLUX_EPSILON;
    pool = pool(carried);
end

% Test what survived, one reaction at a time.
pool = sort(pool);
isEssential = false(numel(pool),1);
for k = 1:numel(pool)
    idx = findRxnIDs(constrained, pool{k});
    knockout = constrained;
    knockout.lb(idx) = 0;
    knockout.ub(idx) = 0;
    solution = optimizeCbModel(knockout);
    isEssential(k) = (solution.stat ~= 1);
end

essential = pool(isEssential);
end
