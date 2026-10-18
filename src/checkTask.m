function result = checkTask(model,task)
% Check whether a model can perform a metabolic task.
%
% USAGE:
%
%    result = checkTask(model,task)
%
% INPUTS:
%    model:              COBRA model structure.
%    task:               Task structure from metabolicTask.
%
% OUTPUTS:
%    result:             Structure with fields:
%
%                          * taskID - the task's identifier
%                          * feasible - whether the constrained model had a
%                            solution
%                          * ok - whether that is the expected outcome,
%                            i.e. feasible unless the task was marked
%                            shouldFail
%                          * missing - metabolite IDs the task named that
%                            the model does not have. A task naming a
%                            missing metabolite is reported rather than
%                            checked, so one unusable task does not stop a
%                            list from being worked through.
%
% NOTE:
%    The objective is zeroed before solving. With nothing to optimise, an
%    optimal status means only that the constraints can be satisfied at once,
%    which is exactly the question a task asks.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

absent = missingMetabolites(model,task);
if ~isempty(absent)
    result.taskID   = task.id;
    result.feasible = false;
    result.ok       = task.shouldFail;
    result.missing  = absent;
    return
end

constrained = taskConstraints(model,task);
constrained.c(:) = 0;

solution = optimizeCbModel(constrained);
feasible = (solution.stat == 1);

result.taskID   = task.id;
result.feasible = feasible;
result.ok       = (feasible ~= task.shouldFail);
result.missing  = {};
end
