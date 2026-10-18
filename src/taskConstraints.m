function [constrained,addedRxns] = taskConstraints(model,task)
% Apply a metabolic task to a model, returning the constrained copy.
%
% USAGE:
%
%    [constrained,addedRxns] = taskConstraints(model,task)
%
% INPUTS:
%    model:              COBRA model structure. It is not modified; MATLAB
%                        passes structs by value, so the caller's model is
%                        untouched.
%    task:               Task structure from metabolicTask.
%
% OUTPUTS:
%    constrained:        The model with every existing boundary reaction
%                        closed, supply and removal reactions opened for
%                        the task's inputs and outputs, and the task's
%                        equations added.
%    addedRxns:          IDs of the reactions this function added, as a
%                        cell array of char. They are all prefixed
%                        'spectra_task_'.
%
% NOTE:
%    The task defines its own medium, so the model's own exchanges are closed
%    first; otherwise they could supply something the task never granted.
%
%    An input is added as '--> met' (stoichiometric coefficient +1) rather
%    than as a conventional exchange written 'met <=>' (coefficient -1), so
%    its flux is the uptake rate directly rather than the negated rate an
%    exchange would carry. Outputs are added as 'met -->' for the same
%    reason. Both sides therefore want *positive* flux, which is why the
%    default bounds are [0 1000] on an input as well as an output. There is
%    no negative lower bound anywhere here: the direction is carried by the
%    sign of the coefficient, not by the sign of the bound.
%
%    checkDuplicate is passed explicitly even though COBRA already defaults
%    it to false, so this changes nothing today. It is pinned because these
%    reactions are exactly what a duplicate check would look at - the model
%    usually has an exchange on the same metabolite already - and because
%    every existing boundary reaction has just been closed. Were the check
%    ever to fire, the task would be left with no way to take anything up
%    and would fail with nothing to explain why.
%
%    A metabolite named as both an input and an output may flow either way:
%    both reactions are built. That is what naming it twice says.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

absent = missingMetabolites(model,task);
if ~isempty(absent)
    error('taskConstraints:missingMetabolites', ...
        'Task ''%s'' names %d metabolite(s) the model does not have, for example %s.', ...
        task.id, numel(absent), strjoin(absent(1:min(5,end))', ', '));
end

constrained = model;

% Close every existing boundary reaction.
isExc = findExcRxns(constrained);
constrained.lb(isExc) = 0;
constrained.ub(isExc) = 0;

addedRxns = {};

for k = 1:numel(task.inputs)
    rxnID = sprintf('spectra_task_in_%d_%s', k, task.inputs{k});
    constrained = addReaction(constrained, rxnID, ...
        'metaboliteList', task.inputs(k), ...
        'stoichCoeffList', 1, ...
        'lowerBound', task.inputBounds(k,1), ...
        'upperBound', task.inputBounds(k,2), ...
        'checkDuplicate', false, ...
        'printLevel', 0);
    addedRxns{end+1,1} = rxnID; %#ok<AGROW>
end

for k = 1:numel(task.outputs)
    rxnID = sprintf('spectra_task_out_%d_%s', k, task.outputs{k});
    constrained = addReaction(constrained, rxnID, ...
        'metaboliteList', task.outputs(k), ...
        'stoichCoeffList', -1, ...
        'lowerBound', task.outputBounds(k,1), ...
        'upperBound', task.outputBounds(k,2), ...
        'checkDuplicate', false, ...
        'printLevel', 0);
    addedRxns{end+1,1} = rxnID; %#ok<AGROW>
end

for k = 1:numel(task.equations)
    equation = task.equations(k);
    rxnID = sprintf('spectra_task_equ_%d', k);
    constrained = addReaction(constrained, rxnID, ...
        'metaboliteList', equation.mets, ...
        'stoichCoeffList', equation.coeffs, ...
        'lowerBound', equation.bounds(1), ...
        'upperBound', equation.bounds(2), ...
        'checkDuplicate', false, ...
        'printLevel', 0);
    addedRxns{end+1,1} = rxnID; %#ok<AGROW>
end
end
