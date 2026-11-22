function result = gapfillCommunity(community, varargin)
% Gap-fills every unit of a community at once, from their own databases
%
% USAGE:
%
%    result = gapfillCommunity(community, varargin)
%
% INPUTS:
%    community:        Structure returned by buildCommunityModel, whose
%                      units each hold a database of reactions they may
%                      draw on. That is what passing 'databases' to
%                      buildCommunityModel arranges: the unit then holds
%                      the whole database, and the reactions only the
%                      database has are the candidates to be filled in.
%
% OPTIONAL INPUTS (name/value pairs):
%    coreRxns:         Cell array of reactions the result must keep, named
%                      as the community model names them, that is tagged by
%                      their unit. (Default: none beyond the anchors)
%    anchorsAreCore:   Whether each unit's anchor reaction is required.
%                      This is what makes every unit grow, rather than only
%                      the ones that happen to be cheapest to grow.
%                      (Default: true)
%    weights:          Struct array with fields rxn and weight, overriding
%                      the cost of keeping particular reactions. The usual
%                      reason is sequence evidence: a reaction the
%                      organism's genome supports should cost less than one
%                      it does not. (Default: none)
%    exchangeWeight:   What a community exchange costs. Left out it follows
%                      the formulation: 1 for 'minNetMILP', which counts
%                      reactions, and 0 for 'minNetLP', which sums flux and
%                      would otherwise charge a unit for the uptake of
%                      everything that shares a pool with it. (Default: by
%                      formulation)
%    tol:              Minimum flux a core reaction must carry.
%                      (Default: 1e-4)
%    probType:         'minNetMILP' (Default) or 'minNetLP'. The LP
%                      minimises total weighted flux rather than a count,
%                      so it is a relaxation: quick, and usually larger.
%                      'tradeOff' is refused for the same reason
%                      gapfillForGrowth refuses it - it requires every
%                      included reaction to carry at least tol, and growth
%                      needs trace fluxes far below that.
%    consistencyCheck: Whether to drop the blocked reactions of the joined
%                      model before gap-filling. Worth its cost: a database
%                      reaction that cannot carry flux in this community is
%                      not a candidate, and removing it shrinks the
%                      problem. (Default: true)
%    keepDraft:        Whether to keep every reaction a unit already had,
%                      whether or not the solution uses it. This is what
%                      separates gap-filling from extraction: the
%                      formulation returns the smallest network meeting the
%                      requirement, which would otherwise throw away parts
%                      of the draft that carry no flux under this medium.
%                      (Default: true)
%    solveTime:        Maximum runtime for the solve, in seconds.
%                      (Default: 300)
%    seed:             Seed for spectraME's random coefficients. Leave it
%                      out to draw them unseeded. (Default: none)
%
% OUTPUTS:
%    result:           Structure with fields:
%
%                        * community - the gap-filled community, with its
%                          bookkeeping narrowed to the reactions that
%                          survived
%                        * added - cell array in the order of
%                          community.organisms, each holding the reactions
%                          that unit took from its database, named as the
%                          unit's own model names them
%                        * core - the reactions the extraction was required
%                          to keep
%                        * blockedCore - core reactions the consistency
%                          check found could not carry flux, and which were
%                          therefore not required
%                        * summary - a readable multi-line summary, as char
%
% EXAMPLE:
%
%    ::
%
%        medium = struct('met', {'glc-D[e]','nh4[e]','pi[e]','h2o[e]'}, ...
%                        'lb',  {-10, -1000, -1000, -1000}, ...
%                        'ub',  {1000, 1000, 1000, 1000});
%        community = buildCommunityModel({draftA, draftB}, {'a','b'}, ...
%            'databases', universal, 'poolMedium', medium);
%        result = gapfillCommunity(community);
%        disp(result.summary)
%        units = communityDecompose(result.community);
%
% NOTE:
%    Gap-filling the units together rather than one at a time lets a gap in
%    one be closed by another's secretion instead of by a new reaction,
%    which is why the joint result is the smaller one. What a unit already
%    has costs nothing, so the solver is free to keep all of it; what only
%    its database has costs one, so it adds as little as it can. Transports
%    between a unit and a pool cost nothing either, because in the other
%    pooling mode the units simply share the compartment, and charging for
%    cross-feeding in one mode but not the other would make the two
%    disagree.
%
%    Follows the community-scale gap-filling of S, P. K. et al. (2026),
%    bioRxiv 2026.04.02.716249, whose community formulation in turn follows
%    Giannari, D., Ho, C. H. and Mahadevan, R. (2021), PLOS Computational
%    Biology 17(6), e1009060.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: buildCommunityModel, communityDecompose, gapfillForGrowth

if ~isstruct(community) || ~isfield(community, 'units')
    error('gapfillCommunity:badCommunity', ...
        'community must be the structure buildCommunityModel returns.');
end

parser = inputParser();
parser.addParameter('coreRxns', {}, @(x) iscell(x) || ischar(x));
parser.addParameter('anchorsAreCore', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('weights', [], @(x) isstruct(x) || isempty(x));
parser.addParameter('exchangeWeight', [], @isnumeric);
parser.addParameter('tol', 1e-4, @isnumeric);
parser.addParameter('probType', 'minNetMILP', @ischar);
parser.addParameter('consistencyCheck', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('keepDraft', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('solveTime', 300, @isnumeric);
parser.addParameter('seed', [], @isnumeric);
parser.parse(varargin{:});
opt = parser.Results;

if ~ismember(opt.probType, {'minNetMILP', 'minNetLP'})
    error('gapfillCommunity:badProbType', ...
        ['probType must be ''minNetMILP'' or ''minNetLP'', not ''%s''. ' ...
         '''tradeOff'' requires every included reaction to carry at least ' ...
         'tol, while growth needs trace fluxes far below that, so it is ' ...
         'infeasible at every tolerance.'], opt.probType);
end

tol = opt.tol;
model = community.model;
nUnits = numel(community.organisms);

if all(cellfun(@isempty, community.databaseReactions))
    warning('gapfillCommunity:noDatabase', ...
        ['No unit of this community holds a database, so there is nothing ' ...
         'to fill a gap with and every reaction costs one. The result will ' ...
         'be an extraction of the community rather than a gap-fill. Pass ' ...
         '''databases'' to buildCommunityModel.']);
end

%% what the result must keep
core = opt.coreRxns;
if ischar(core)
    core = {core};
end
core = core(:);
if opt.anchorsAreCore
    core = [core; community.biomassReactions(:)];
end
core = unique(core, 'stable');
unknown = core(findRxnIDs(model, core) == 0);
if ~isempty(unknown)
    error('gapfillCommunity:unknownCore', ...
        ['These core reactions are not in the community model: %s. They ' ...
         'are named as the community model names them, that is tagged by ' ...
         'their unit.'], strjoin(unknown(1:min(5, end))', ', '));
end

%% what each reaction costs
exchangeWeight = opt.exchangeWeight;
if isempty(exchangeWeight)
    exchangeWeight = double(strcmp(opt.probType, 'minNetMILP'));
end
costOf = communityWeights(community, exchangeWeight, opt.weights);

%% drop what cannot carry flux
working = community;
blockedCore = cell(0, 1);
if opt.consistencyCheck
    % With coupling, because a community's whole point is that a unit
    % which cannot grow cannot run its metabolism either. A database
    % reaction that can only carry flux by ignoring that is not a
    % candidate.
    consistent = spectraCC(model, tol, 'stoichiometry', opt.seed, ...
        getCobraSolverParams('LP', 'feasTol'), true);
    keepIDs = model.rxns(consistent);
    fprintf(['gapfillCommunity: %d of %d reactions can carry flux in this ' ...
             'community.\n'], numel(keepIDs), numel(model.rxns));
    blockedCore = setdiff(core, keepIDs, 'stable');
    if ~isempty(blockedCore)
        warning('gapfillCommunity:blockedCore', ...
            ['%d core reactions cannot carry flux in this community and ' ...
             'were not required, for example %s.'], numel(blockedCore), ...
            strjoin(blockedCore(1:min(3, end))', ', '));
    end
    core = setdiff(core, blockedCore, 'stable');
    working = communityWithModel(community, rebuildModel(model, keepIDs));
end

%% the fill
present = working.model.rxns;
weights = zeros(numel(present), 1);
for i = 1:numel(present)
    weights(i) = costOf(present{i});
end
candidates = present(weights ~= 0);
if strcmp(opt.probType, 'minNetMILP')
    groups = candidates;
else
    % minNetLP has no binaries to group; a zero weight is how it is told
    % that a draft reaction is free.
    groups = {};
end

fprintf(['gapfillCommunity: filling %d units with %s; %d candidate ' ...
         'reactions, %d core.\n'], nUnits, opt.probType, numel(candidates), ...
        numel(core));
filled = spectraME(working.model, findRxnIDs(working.model, core), tol, ...
    'stoichiometry', weights, 1, 'coreDirection', opt.probType, ...
    opt.solveTime, false, {}, opt.seed, groups, [], [], true);

kept = filled.rxns;
if opt.keepDraft
    drafted = vertcat(working.units.draftRxns);
    if isempty(drafted)
        drafted = cell(0, 1);
    end
    fprintf('gapfillCommunity: keeping %d draft reactions the solution does not use.\n', ...
        numel(setdiff(drafted, kept)));
    kept = unique([kept(:); drafted(:)], 'stable');
end
result.community = communityWithModel(working, rebuildModel(working.model, kept));

%% what each unit took
added = cell(nUnits, 1);
for k = 1:nUnits
    organism = community.organisms{k};
    candidatesOfUnit = community.databaseReactions{k};
    if isempty(candidatesOfUnit)
        added{k} = cell(0, 1);
        continue
    end
    tagged = strcat(candidatesOfUnit(:), '__', organism);
    added{k} = sort(candidatesOfUnit(ismember(tagged, kept)));
    added{k} = added{k}(:);
end
result.added = added;
result.core = core;
result.blockedCore = blockedCore;
result.summary = buildSummary(result, community);
end


function costOf = communityWeights(community, exchangeWeight, overrides)
% Returns a map from every reaction of the community to what keeping it costs.
%
% What a unit already has costs nothing, so the fill is free to keep all of
% it; what only its database has costs one, so the solver adds as little as
% it can. Transports cost nothing because in 'shared' mode they do not
% exist, and charging for cross-feeding in one mode but not the other would
% make the two disagree.
rxns = community.model.rxns;
costOf = containers.Map(rxns, ones(numel(rxns), 1));

for k = 1:numel(community.units)
    for field = {'draftRxns', 'rxns'}
        ids = community.units(k).(field{1});
        for i = 1:numel(ids)
            if strcmp(field{1}, 'rxns') && ~startsWith(ids{i}, 'TR_')
                continue
            end
            if isKey(costOf, ids{i})
                costOf(ids{i}) = 0;
            end
        end
    end
end

for i = 1:numel(community.communityExchanges)
    id = community.communityExchanges{i};
    if isKey(costOf, id)
        costOf(id) = exchangeWeight;
    end
end

for i = 1:numel(overrides)
    id = overrides(i).rxn;
    if ~isKey(costOf, id)
        error('gapfillCommunity:unknownWeight', ...
            'No reaction %s in the community model, so it cannot be weighted.', id);
    end
    costOf(id) = overrides(i).weight;
end
end


function reduced = rebuildModel(model, keepIDs)
% Returns the model holding only the named reactions.
drop = setdiff(model.rxns, keepIDs);
if isempty(drop)
    reduced = model;
else
    reduced = removeRxns(model, drop);
end
end


function text = buildSummary(result, community)
% Returns a readable account of what the fill did.
total = sum(cellfun(@numel, result.added));
lines = {sprintf(['%d reactions across %d units; %d added, %d core ' ...
                  'required.'], numel(result.community.model.rxns), ...
                 numel(community.organisms), total, numel(result.core))};
if ~isempty(result.blockedCore)
    lines{end+1} = sprintf('  %d core reactions were blocked and not required.', ...
        numel(result.blockedCore));
end
for k = 1:numel(community.organisms)
    lines{end+1} = sprintf('  %-16s %d added', community.organisms{k}, ...
        numel(result.added{k})); %#ok<AGROW>
end
text = strjoin(lines, newline);
end
