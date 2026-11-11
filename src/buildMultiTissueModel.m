function body = buildMultiTissueModel(model, tissues, pools, environment, varargin)
% Replicates one model across tissues and joins them into one network
%
% USAGE:
%
%    body = buildMultiTissueModel(model, tissues, pools, environment, varargin)
%
% INPUTS:
%    model:              The reference network every tissue starts from,
%                        typically a human genome-scale reconstruction. Give
%                        one COBRA model structure and every tissue is a
%                        copy of it, which is the usual case: what then
%                        distinguishes the tissues is their bounds, their
%                        core reactions and which pool they reach. Give a
%                        cell array of one model per tissue, in the order of
%                        tissues, when the tissues are genuinely different
%                        networks.
%    tissues:            Cell array of the tissue names, which become the
%                        tag on every reaction and metabolite.
%    pools:              Struct array with fields name and units, saying
%                        which tissues meet in each pool. A pool with no
%                        tissue of its own - blood, say - is simply a name
%                        that appears here and not in tissues.
%    environment:        Cell array of the pools that reach the
%                        environment. Usually just the blood.
%
% OPTIONAL INPUTS (name/value pairs):
%    anchorReactions:    The reaction each tissue's coupling gates on, as
%                        named in model: a char to use the same one for
%                        every tissue, or a cell array of one per tissue in
%                        the order of tissues. A tissue is rarely growing,
%                        so this is usually a maintenance or demand reaction
%                        rather than biomass, and it has to be said rather
%                        than guessed. (Default: each model's objective,
%                        else a reaction named biomass*)
%    varargin:           Everything else is passed to buildCommunityModel,
%                        which documents it.
%
% OUTPUTS:
%    body:               Struct describing the joined model, exactly as
%                        buildCommunityModel returns it. The tissues are
%                        reported in its organisms field.
%
% EXAMPLE:
%
%    ::
%
%        recon = readCbModel('Recon3D.mat');
%        pools = struct('name', {'blood', 'gut'}, ...
%                       'units', {{'liver','muscle','intestine'}, ...
%                                 {'intestine'}});
%        body = buildMultiTissueModel(recon, ...
%            {'liver','muscle','intestine'}, pools, {'blood'}, ...
%            'anchorReactions', 'biomass_maintenance', ...
%            'poolMedium', medium);
%
% NOTE:
%    Nothing in the construction distinguishes a tissue from an organism, so
%    this is buildCommunityModel with the replication done for you. Every
%    tissue is a copy of the same network, so what distinguishes them is
%    their bounds, their core reactions and which pool they can reach. The
%    copies are made here rather than by the caller because the reference
%    model is usually large and copying it correctly matters.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: buildCommunityModel

if nargin < 4
    error('buildMultiTissueModel:tooFewArguments', ...
        ['model, tissues, pools and environment must all be given. A pool ' ...
         'that merely links two tissues has to stay closed, so which of ' ...
         'them reaches the environment cannot be guessed.']);
end
if isempty(tissues) || ~iscell(tissues)
    error('buildMultiTissueModel:noTissues', ...
        'tissues must be a non-empty cell array of names.');
end
tissues = tissues(:);

% One model per tissue, or the same model handed over once per tissue. The
% builder reads its inputs and never modifies them, so handing it the same
% structure several times costs nothing.
if iscell(model)
    copies = model(:);
    if numel(copies) ~= numel(tissues)
        error('buildMultiTissueModel:modelCountMismatch', ...
            ['Got %d models for %d tissues. Pass one model per tissue, in ' ...
             'the same order, or a single model for every tissue to be a ' ...
             'copy of.'], numel(copies), numel(tissues));
    end
else
    copies = repmat({model}, numel(tissues), 1);
end

[anchors, passThrough] = takeAnchors(varargin, tissues);

body = buildCommunityModel(copies, tissues, 'pools', pools, ...
    'environment', environment, 'biomassReactions', anchors, passThrough{:});
end


function [anchors, rest] = takeAnchors(args, tissues)
% Pulls anchorReactions out of the arguments, leaving the rest alone.
anchors = {};
rest = args;

at = find(strcmpi(args(1:2:end), 'anchorReactions'), 1);
if isempty(at)
    return
end
at = 2 * at - 1;
given = args{at + 1};
rest = args([1:at-1, at+2:end]);

if isempty(given)
    return
end
if ischar(given)
    % The usual case when every tissue is a copy of one reconstruction:
    % they all gate on the same reaction of it.
    anchors = repmat({given}, numel(tissues), 1);
    return
end
if ~iscell(given)
    error('buildMultiTissueModel:badAnchors', ...
        ['anchorReactions must be a reaction name, or a cell array of one ' ...
         'per tissue in the order of tissues.']);
end
if numel(given) ~= numel(tissues)
    error('buildMultiTissueModel:anchorCountMismatch', ...
        'Got %d anchorReactions for %d tissues.', numel(given), numel(tissues));
end
anchors = given(:);
end
