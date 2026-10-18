function absent = missingMetabolites(model,task)
% Return the metabolites a task names that a model does not have.
%
% USAGE:
%   absent = missingMetabolites(model,task)
%
% INPUTS:
%   model: COBRA model structure to check against.
%   task:  Task structure from metabolicTask.
%
% OUTPUT:
%   absent: The absent metabolite IDs, in the order inputs, outputs,
%           equations, with duplicates removed and that order preserved.
%
% NOTE:
%   A task naming a metabolite the model lacks cannot be checked against it.
%   checkTask reports that rather than raising, so one unusable task in a
%   list does not stop the rest from being checked.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

named = [task.inputs(:); task.outputs(:)];
for k = 1:numel(task.equations)
    named = [named; task.equations(k).mets(:)]; %#ok<AGROW>
end

if isempty(named)
    absent = {};
    return
end

% unique(...,'stable') keeps the inputs-outputs-equations order.
named = unique(named, 'stable');
absent = named(~ismember(named, model.mets));
absent = absent(:);
end
