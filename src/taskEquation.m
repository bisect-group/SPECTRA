function equation = taskEquation(mets,coeffs,bounds)
% Build a reaction that a task adds to the model while it is being checked.
%
% USAGE:
%
%    equation = taskEquation(mets,coeffs,bounds)
%
% INPUTS:
%    mets:               Metabolite IDs, as a cell array of char.
%    coeffs:             Stoichiometric coefficient for each, as a numeric
%                        vector of the same length. Negative coefficients
%                        consume, positive produce.
%
% OPTIONAL INPUTS:
%    bounds:             Flux bounds as [lower upper] (Default: [-1000
%                        1000]). A positive lower bound is what forces the
%                        task to actually happen.
%
% OUTPUTS:
%    equation:           Structure with fields mets, coeffs and bounds.
%
% EXAMPLE:
%
%    % ATP hydrolysis, forced to run at a rate of at least 1:
%    equation = taskEquation({'atp_c','h2o_c','adp_c','pi_c','h_c'}, ...
%                            [-1 -1 1 1 1], [1 1000]);
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('bounds', 'var') || isempty(bounds)
    bounds = [-1000 1000];
end

if ischar(mets)
    mets = {mets};
end
if ~iscellstr(mets)
    error('taskEquation:badMets', ...
        'mets must be a cell array of metabolite IDs.');
end
mets = mets(:);
coeffs = coeffs(:);
if numel(mets) ~= numel(coeffs)
    error('taskEquation:lengthMismatch', ...
        'mets has %d entries but coeffs has %d.', numel(mets), numel(coeffs));
end
if numel(bounds) ~= 2
    error('taskEquation:badBounds', ...
        'bounds must be [lower upper].');
end
bounds = bounds(:)';
if bounds(1) > bounds(2)
    error('taskEquation:emptyBounds', ...
        'the lower bound %g is above the upper bound %g.', bounds(1), bounds(2));
end

equation.mets = mets;
equation.coeffs = coeffs;
equation.bounds = bounds;
end
