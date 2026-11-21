function narrowed = communityWithModel(community, model)
% Carries a community's bookkeeping over to a reduced model
%
% USAGE:
%
%    narrowed = communityWithModel(community, model)
%
% INPUTS:
%    community: Structure returned by buildCommunityModel.
%    model:     COBRA model holding some subset of the community's
%                reactions, as a consistency check or an extraction leaves
%                behind.
%
% OUTPUTS:
%    narrowed:  The same structure over the new model, with every list of
%               reactions and metabolites narrowed to what survived, so
%               that a unit's reactions, its database candidates and the
%               pool exchanges all still name things that exist.
%
% NOTE:
%    biomassReactions is left alone even where an anchor did not survive.
%    Which reaction anchors a unit is a fact about the community, not about
%    this particular reduction of it, and silently forgetting it would hide
%    the more interesting news that the reduction dropped it.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: buildCommunityModel, communitySubset, gapfillCommunity

if ~isstruct(community) || ~isfield(community, 'units')
    error('communityWithModel:badCommunity', ...
        'community must be the structure buildCommunityModel returns.');
end

presentRxns = model.rxns;
presentMets = model.mets;

narrowed = community;
narrowed.model = model;
narrowed.communityExchanges = ...
    community.communityExchanges(ismember(community.communityExchanges, presentRxns));

for p = 1:numel(narrowed.poolMetabolites)
    alive = ismember(narrowed.poolMetabolites(p).mets, presentMets);
    narrowed.poolMetabolites(p).mets = narrowed.poolMetabolites(p).mets(alive);
    narrowed.poolMetabolites(p).keys = narrowed.poolMetabolites(p).keys(alive);
end

for k = 1:numel(narrowed.units)
    organism = narrowed.units(k).name;
    narrowed.units(k).rxns = keepPresent(narrowed.units(k).rxns, presentRxns);
    narrowed.units(k).draftRxns = keepPresent(narrowed.units(k).draftRxns, presentRxns);
    narrowed.units(k).databaseRxns = ...
        keepPresent(narrowed.units(k).databaseRxns, presentRxns);
    % draftReactions and databaseReactions hold the unit's own untagged
    % names, so they are filtered on what the tag would have made of them.
    narrowed.draftReactions{k} = keepTagged(community.draftReactions{k}, ...
        organism, presentRxns);
    narrowed.databaseReactions{k} = keepTagged(community.databaseReactions{k}, ...
        organism, presentRxns);
end
end


function kept = keepPresent(ids, present)
% Returns the identifiers that the reduced model still has.
if isempty(ids)
    kept = cell(0, 1);
    return
end
kept = ids(ismember(ids, present));
kept = kept(:);
end


function kept = keepTagged(ids, organism, present)
% Returns the untagged names whose tagged form the reduced model still has.
if isempty(ids)
    kept = cell(0, 1);
    return
end
kept = ids(ismember(strcat(ids(:), '__', organism), present));
kept = kept(:);
end
