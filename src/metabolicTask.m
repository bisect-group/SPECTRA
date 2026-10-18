function task = metabolicTask(id,varargin)
% Build a statement of something a metabolic network must be able to do.
%
% USAGE:
%
%    task = metabolicTask(id,'name',value,...)
%
% INPUTS:
%    id:                 A short identifier, as char.
%
% OPTIONAL INPUTS:
%    Given as name/value pairs.
%
%    'description':      What the task is for, in words. (Default: '')
%    'inputs':           Metabolite IDs that may be taken up, as a cell
%                        array of char. (Default: {})
%    'inputBounds':      Bounds on each uptake flux, as an n-by-2 matrix of
%                        [lower upper], or a single [lower upper] applied
%                        to all. (Default: [0 1000] each, i.e. allowed but
%                        not required) These are *not* exchange reaction
%                        bounds, so the lower bound is 0 rather than -1000.
%                        taskConstraints adds an input as '--> met' rather
%                        than as an exchange 'met <=>', so positive flux is
%                        uptake and the bounds read as uptake rates. See
%                        taskConstraints for why.
%    'outputs':          Metabolite IDs that may be produced. (Default: {})
%    'outputBounds':     Bounds on each production flux, in the same
%                        shapes. A positive lower bound *requires* the
%                        production. (Default: [0 1000] each)
%    'equations':        Reactions added for the duration of the check, as
%                        a struct array from taskEquation. (Default: [])
%    'shouldFail':       Whether the task is expected to be infeasible.
%                        (Default: false). These are the negative controls
%                        of a task list: a model that can do them is
%                        leaking.
%
% OUTPUTS:
%    task:               Structure with fields id, description, inputs,
%                        inputBounds, outputs, outputBounds, equations and
%                        shouldFail.
%
% NOTE:
%    A task with no positive lower bound anywhere is satisfied by the
%    all-zero flux vector, so it passes trivially. Force the task through an
%    output's lower bound or a taskEquation's. A warning says so if nothing
%    in the task requires flux.
%
% EXAMPLE:
%
%    hydrolysis = taskEquation({'atp_c','h2o_c','adp_c','pi_c','h_c'}, ...
%                              [-1 -1 1 1 1], [1 1000]);
%    task = metabolicTask('atp_aerobic', ...
%        'inputs',  {'glc__D_e','o2_e','pi_e','h2o_e','h_e'}, ...
%        'outputs', {'co2_e','h2o_e','h_e'}, ...
%        'equations', hydrolysis);
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~ischar(id) || isempty(id)
    error('metabolicTask:badId', 'id must be a non-empty char.');
end

p = inputParser;
p.addParameter('description','',@ischar);
p.addParameter('inputs',{});
p.addParameter('inputBounds',[]);
p.addParameter('outputs',{});
p.addParameter('outputBounds',[]);
p.addParameter('equations',[]);
p.addParameter('shouldFail',false);
p.parse(varargin{:});
r = p.Results;

[inputs,inputBounds]   = normaliseSide(r.inputs, r.inputBounds, 'inputs');
[outputs,outputBounds] = normaliseSide(r.outputs, r.outputBounds, 'outputs');

equations = r.equations;
if ~isempty(equations)
    needed = {'mets','coeffs','bounds'};
    if ~isstruct(equations) || ~all(isfield(equations,needed))
        error('metabolicTask:badEquations', ...
            'equations must be a struct array built by taskEquation.');
    end
    equations = equations(:);
end

task.id           = id;
task.description  = r.description;
task.inputs       = inputs;
task.inputBounds  = inputBounds;
task.outputs      = outputs;
task.outputBounds = outputBounds;
task.equations    = equations;
task.shouldFail   = logical(r.shouldFail);

if ~taskIsForced(task) && ~task.shouldFail
    warning('metabolicTask:notForced', ...
        ['Task ''%s'' has no positive lower bound on any output or ' ...
         'equation, so the all-zero flux vector satisfies it and it will ' ...
         'pass trivially.'], id);
end
end


function [mets,bounds] = normaliseSide(mets,bounds,side)
% Put one side of a task into cellstr plus an n-by-2 bounds matrix.
if ischar(mets)
    mets = {mets};
end
if isempty(mets)
    mets = {};
    bounds = zeros(0,2);
    return
end
if ~iscellstr(mets)
    error('metabolicTask:badSide', ...
        '%s must be a cell array of metabolite IDs.', side);
end
mets = mets(:);

if isempty(bounds)
    bounds = repmat([0 1000], numel(mets), 1);
elseif isequal(size(bounds),[1 2])
    bounds = repmat(bounds, numel(mets), 1);
elseif size(bounds,1) ~= numel(mets) || size(bounds,2) ~= 2
    error('metabolicTask:badSideBounds', ...
        ['%sBounds must be [lower upper] or an n-by-2 matrix with one ' ...
         'row per metabolite; got %d-by-%d for %d metabolites.'], ...
        side, size(bounds,1), size(bounds,2), numel(mets));
end
if any(bounds(:,1) > bounds(:,2))
    error('metabolicTask:emptySideBounds', ...
        'a lower bound in %sBounds is above its upper bound.', side);
end
end
