function units = communityDecompose(community)
% Splits a community back into one model per unit
%
% USAGE:
%
%    units = communityDecompose(community)
%
% INPUTS:
%    community: Structure returned by buildCommunityModel, or by
%               gapfillCommunity after it has been gap-filled.
%
% OUTPUTS:
%    units:     Struct array with one entry per unit, holding:
%
%                 * name - the unit's name
%                 * model - its own network, with the tags stripped and an
%                   exchange reaction wherever the unit met a pool
%                 * anchor - its anchor reaction, untagged
%
% NOTE:
%    What was a transport into the pool becomes an exchange. Dropping the
%    pool side of ``x[e] -> x[u]`` leaves ``x[e] ->``, which is an exchange
%    carrying the transport's own bounds, so the rates measured for the
%    unit in the community survive the split. Where the units shared the
%    external compartment there are no transports, and an exchange is added
%    for every pool metabolite the unit touches instead.
%
%    This is how a community gap-fill is cashed out: the units are filled
%    together, so that a gap in one can be closed by another's secretion
%    rather than by a new reaction, and then taken apart again into models
%    that stand on their own.
%
% EXAMPLE:
%
%    ::
%
%        filled = gapfillCommunity(community);
%        units  = communityDecompose(filled.community);
%        writeCbModel(units(1).model, 'format', 'mat', ...
%            'fileName', [units(1).name '.mat']);
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: buildCommunityModel, gapfillCommunity, communitySubset

if ~isstruct(community) || ~isfield(community, 'units')
    error('communityDecompose:badCommunity', ...
        'community must be the structure buildCommunityModel returns.');
end

model = community.model;
poolMets = vertcat(community.poolMetabolites.mets);
if isempty(poolMets)
    poolMets = cell(0, 1);
end
shared = strcmp(community.mode, 'shared');

units = struct('name', {}, 'model', {}, 'anchor', {});
for k = 1:numel(community.units)
    organism = community.units(k).name;
    own = community.units(k).rxns;
    own = own(ismember(own, model.rxns));
    if isempty(own)
        error('communityDecompose:emptyUnit', ...
            'Unit %s has no reactions left, so there is nothing to split out.', ...
            organism);
    end

    sub = removeRxns(model, setdiff(model.rxns, own));
    dangling = poolMets(ismember(poolMets, sub.mets));
    if shared
        % No transports to turn into exchanges, so make one per pool
        % metabolite the unit touches.
        if ~isempty(dangling)
            sub = addExchangeRxn(sub, dangling, ...
                -1000 * ones(numel(dangling), 1), 1000 * ones(numel(dangling), 1));
        end
    elseif ~isempty(dangling)
        sub = removeMetabolites(sub, dangling, false);
    end

    sub = stripTag(sub, organism);
    sub.modelID = organism;
    sub.description = organism;

    units(k).name = organism;
    units(k).model = sub;
    units(k).anchor = untag(community.units(k).anchor, organism);
end
end


function model = stripTag(model, organism)
% Takes the unit's tag off everything that wears it.
model.rxns = untagAll(model.rxns, organism);
model.mets = untagAll(model.mets, organism);
for field = {'comps', 'genes'}
    name = field{1};
    if isfield(model, name) && iscell(model.(name))
        model.(name) = untagAll(model.(name), organism);
    end
end
for field = {'grRules', 'rxnNames', 'metNames'}
    name = field{1};
    if isfield(model, name) && iscell(model.(name))
        model.(name) = regexprep(model.(name), ['__' regexptranslate('escape', organism) '\>'], '');
    end
end
end


function ids = untagAll(ids, organism)
% Takes the tag off a cell array of identifiers.
for i = 1:numel(ids)
    ids{i} = untag(ids{i}, organism);
end
end


function id = untag(id, organism)
% Takes the tag off one identifier, whether it sits inside [] or at the end.
suffix = ['__' organism];
token = regexp(id, '^(.*)\[(.*)\]$', 'tokens', 'once');
if ~isempty(token)
    if endsWith(token{2}, suffix)
        id = sprintf('%s[%s]', token{1}, token{2}(1:end-numel(suffix)));
    end
    return
end
if endsWith(id, suffix)
    id = id(1:end-numel(suffix));
end
end
