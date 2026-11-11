function community = buildCommunityModel(models, organisms, varargin)
% Joins models into one network whose units trade metabolites via a pool
%
% USAGE:
%
%    community = buildCommunityModel(models, organisms, varargin)
%
% INPUTS:
%    models:             Cell array of COBRA model structures, one per unit.
%                        They are not modified.
%    organisms:          Cell array of short names, one per model, in the
%                        same order. These become the tag on every reaction
%                        and metabolite, so they must be unique and must not
%                        contain '__'.
%
% OPTIONAL INPUTS (name/value pairs):
%    mode:               How to build the pool, 'pooled' or 'shared'.
%                        'pooled' gives each unit a private external
%                        compartment and a transport reaction per
%                        exchangeable metabolite, which can be bounded per
%                        unit, and allows several pools. 'shared' has the
%                        units share the external compartment outright: one
%                        row per exchangeable metabolite for the whole
%                        community, no transports and only one pool. For a
%                        hundred units across a few thousand metabolites
%                        that saves a hundred thousand reactions, which is
%                        the difference between a model that can be built
%                        and one that cannot. (Default: 'pooled')
%    pools:              Struct array with fields name and units, saying
%                        which units meet in each pool. A pool is named by a
%                        string that is not a unit, and a unit may belong to
%                        more than one. Every unit must be in at least one:
%                        a unit in no pool can exchange nothing and cannot
%                        grow. (Default: one pool, named by
%                        sharedCompartment, holding every unit)
%    environment:        Cell array of the pools carrying exchange reactions.
%                        Required when there is more than one pool: a pool
%                        that merely links two units must stay closed, or
%                        they could bypass the rest of the system.
%                        (Default: the single pool)
%    databases:          A model of everything a unit may draw on. Give one
%                        model struct to use it for every unit, a cell array
%                        of one per unit in the order of organisms (with []
%                        where a unit has none), or a struct array with
%                        fields organism and model. The unit then holds the
%                        whole database, bounded by its own model where the
%                        two agree, and the extra reactions are recorded in
%                        databaseReactions as candidates to gap-fill in.
%                        (Default: none)
%    biomassReactions:   Cell array of the anchor reaction of each unit, in
%                        the order of organisms. (Default: the unit's
%                        objective, else a reaction named biomass*)
%    sharedMetabolites:  Cell array of metabolites to leave untagged, so
%                        that every unit uses the same one. Rarely wanted,
%                        but a database with a single 'biomass' metabolite
%                        drained by one community-wide exchange needs it.
%                        (Default: none)
%    poolMedium:         Struct array with fields met, lb and ub, bounding
%                        the exchange between a pool and the environment.
%                        met is the metabolite as the unit models spell it.
%                        Anything unnamed keeps its default bounds, which
%                        are export-only. (Default: none, so the community
%                        takes up nothing)
%    linkBounds:         Struct array bounding what may cross between a unit
%                        and a pool. This is where measured uptake and
%                        secretion rates go. Two ways to name the link:
%
%                          * field rxn - one transport, named
%                            TR_<metabolite>_<pool>__<unit>, which is the
%                            per-unit statement a measured rate makes
%                          * field met - the metabolite as the unit models
%                            spell it, which bounds every unit's transport
%                            of it alike
%
%                        A struct array may mix the two; rxn wins where both
%                        match. Only meaningful when mode is 'pooled', since
%                        otherwise there is nothing to bound. (Default:
%                        [-1000 1000] for every link)
%    sharedCompartment:  Name of the default pool. (Default: 'u')
%    externalCompartment: The compartment whose metabolites are
%                        exchangeable: a char naming one for every unit, or
%                        a cell array of one per unit in the order of
%                        organisms. Left out, each unit's is worked out from
%                        the model, which matters because the units need not
%                        agree - one may call it 'e' and another
%                        'extracellular'. (Default: worked out per unit)
%    metaboliteKey:      Function handle. Given a metabolite identifier,
%                        returns the name of the thing it *is*, so that two
%                        units naming the same compound differently still
%                        meet in the pool. The default strips the
%                        compartment, which resolves 'glc-D[e]' and
%                        'glc-D_e' to the same 'glc-D' but cannot reconcile
%                        identifiers from different namespaces. Supply this
%                        when joining models from collections that do not
%                        share one. (Default: strip the compartment)
%    couple:             Whether to tie each unit's reactions to its own
%                        anchor. Without this a unit that is not growing can
%                        still run its metabolism and feed the others.
%                        (Default: true)
%    couplingFactor:     Flux allowed per unit of anchor flux. (Default:
%                        1000, the published value for a biomass anchor)
%    couplingThreshold:  Flux allowed at zero anchor flux. (Default: 0.01)
%    check:              Whether to measure each anchor and report a coupling
%                        that would throttle its unit. One LP for the
%                        community, and one more for each unit that reads
%                        below the point where coupling stops biting.
%                        (Default: true)
%    modelId:            Identifier of the joined model. (Default: 'community')
%
% OUTPUTS:
%    community:          Struct describing the joined model, with fields:
%
%                          * model - the joined COBRA model
%                          * organisms - the unit names, in order
%                          * units - struct array with one entry per unit,
%                            holding name, anchor, anchorSource,
%                            anchorCapacity, externalCompartment, rxns,
%                            draftRxns and databaseRxns
%                          * biomassReactions - the anchor of each unit, in
%                            the order of organisms
%                          * communityExchanges - the pool exchange reactions
%                          * pools - the struct array used
%                          * poolMetabolites - struct array of name, mets
%                            and keys
%                          * mode - how the pool was built, 'pooled' or
%                            'shared'
%                          * draftReactions, databaseReactions - where a
%                            database was supplied, the reactions that came
%                            from the unit's own model and the ones only the
%                            database has, as cell arrays in the order of
%                            organisms
%                          * couplingFactor, couplingThreshold - as used, or
%                            empty when couple was false
%
% EXAMPLE:
%
%    ::
%
%        model = readCbModel('ecoli_core_model.mat');
%        medium = struct('met', {'glc-D[e]','o2[e]','nh4[e]','pi[e]'}, ...
%                        'lb',  {-10, -1000, -1000, -1000}, ...
%                        'ub',  {1000, 1000, 1000, 1000});
%        community = buildCommunityModel({model, model}, {'org1','org2'}, ...
%            'poolMedium', medium);
%        sol = optimizeCbModel(community.model);
%
%        % the same community, sharing the external compartment outright
%        compact = buildCommunityModel({model, model}, {'org1','org2'}, ...
%            'mode', 'shared', 'poolMedium', medium);
%
% NOTE:
%    In 'pooled' mode each unit keeps a private external compartment and a
%    transport reaction moves each exchangeable metabolite between it and
%    the pool. Because the transport belongs to one unit it can be bounded
%    per unit, and several pools can coexist, so the units need not all meet
%    in the same place. In 'shared' mode the units share the external
%    compartment itself and there are no transports, so there is nothing to
%    bound per unit and only one pool is possible.
%
%    The unit's own exchange reactions do not survive: each becomes the link
%    to the pool it corresponds to. Only the exchanges facing the external
%    compartment are replaced - a demand or sink inside the cell is also a
%    reaction with a single metabolite, and an anchor that only consumes
%    looks exactly like one, so dropping every such reaction would delete
%    it. Anything entering or leaving the whole system does so through a
%    pool exchange, and those are export-only until poolMedium says
%    otherwise. A community built without a medium therefore cannot grow,
%    which is deliberate: a joined model that silently supplies everything
%    is the more dangerous default.
%
%    The compartmentalised layout, the biomass coupling and its parameters
%    follow the joint FBA formulation used for minimal microbiomes in Raghu
%    et al. (2024), npj Systems Biology and Applications 10, 46.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

parser = inputParser();
parser.addParameter('mode', 'pooled', @ischar);
parser.addParameter('pools', [], @(x) isstruct(x) || isempty(x));
parser.addParameter('environment', {}, @(x) iscell(x) || ischar(x));
parser.addParameter('databases', [], @(x) isstruct(x) || iscell(x) || isempty(x));
parser.addParameter('biomassReactions', {}, @iscell);
parser.addParameter('sharedMetabolites', {}, @(x) iscell(x) || ischar(x));
parser.addParameter('poolMedium', [], @(x) isstruct(x) || isempty(x));
parser.addParameter('linkBounds', [], @(x) isstruct(x) || isempty(x));
parser.addParameter('sharedCompartment', 'u', @ischar);
parser.addParameter('externalCompartment', {}, @(x) iscell(x) || ischar(x));
parser.addParameter('metaboliteKey', [], @(x) isempty(x) || isa(x, 'function_handle'));
parser.addParameter('couple', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('couplingFactor', 1000, @isnumeric);
parser.addParameter('couplingThreshold', 0.01, @isnumeric);
parser.addParameter('check', true, @(x) islogical(x) || isnumeric(x));
parser.addParameter('modelId', 'community', @ischar);
parser.parse(varargin{:});
opt = parser.Results;

POOLED = 'pooled';
SHARED = 'shared';
if ~ismember(opt.mode, {POOLED, SHARED})
    error('buildCommunityModel:badMode', ...
        'mode must be ''%s'' or ''%s'', not ''%s''.', POOLED, SHARED, opt.mode);
end

%% units and their names
if ~iscell(models) || isempty(models)
    error('buildCommunityModel:noModels', ...
        'models must be a non-empty cell array of COBRA model structures.');
end
models = models(:);
if nargin < 2 || isempty(organisms)
    error('buildCommunityModel:noOrganisms', ...
        ['organisms must be given: the names become the tag on every ' ...
         'reaction and metabolite, so they cannot be guessed.']);
end
organisms = organisms(:);
if numel(organisms) ~= numel(models)
    error('buildCommunityModel:countMismatch', ...
        'Got %d models but %d organism names.', numel(models), numel(organisms));
end
if numel(unique(organisms)) ~= numel(organisms)
    error('buildCommunityModel:duplicateOrganisms', ...
        'Organism names must be unique.');
end
if any(contains(organisms, '__')) || any(cellfun(@isempty, organisms))
    error('buildCommunityModel:badOrganismName', ...
        ['Organism names must be non-empty and must not contain ''__'', ' ...
         'which separates a name from the identifier it tags.']);
end

nUnits = numel(models);

%% pools, and which of them reach the environment
pools = resolvePools(opt.pools, organisms, opt.sharedCompartment, opt.mode, SHARED);
environment = resolveEnvironment(opt.environment, pools);
poolsOf = cell(nUnits, 1);
for k = 1:nUnits
    poolsOf{k} = {pools(cellfun(@(u) ismember(organisms{k}, u), {pools.units})).name};
end

if strcmp(opt.mode, SHARED) && ~isempty(opt.linkBounds)
    error('buildCommunityModel:linkBoundsInShared', ...
        ['linkBounds bounds what crosses between a unit and a pool, but ' ...
         'mode ''%s'' has no such reactions to bound. Use mode ''%s'', or ' ...
         'put the bounds in poolMedium.'], SHARED, POOLED);
end

%% fold each unit's database in, before anything is tagged
% So that the exchange detection and the tagging see one model per unit.
databases = resolveDatabases(opt.databases, organisms);
unitModels = cell(nUnits, 1);
draftRxns = repmat({cell(0, 1)}, nUnits, 1);
databaseRxns = repmat({cell(0, 1)}, nUnits, 1);
for k = 1:nUnits
    if isempty(databases{k})
        unitModels{k} = models{k};
        continue
    end
    [unitModels{k}, draftRxns{k}, databaseRxns{k}] = ...
        mergeWithDatabase(models{k}, databases{k}, organisms{k});
end

%% the anchor each unit's coupling gates on
[anchors, anchorSources] = resolveAnchors(unitModels, organisms, opt.biomassReactions);

%% which compartment each unit exchanges through
externalOf = resolveExternalCompartments(unitModels, organisms, ...
    opt.externalCompartment);

%% the pool identity of every exchangeable metabolite
% Worked out before anything is built, so that a medium can be given in
% whatever identifiers the caller has and still be matched against the
% metabolites that exist.
isExternal = cell(nUnits, 1);
metKeys = cell(nUnits, 1);
keyOfMet = containers.Map('KeyType', 'char', 'ValueType', 'char');
examples = cell(nUnits, 1);
for k = 1:nUnits
    mets = unitModels{k}.mets;
    outside = externalOf{k};
    isExternal{k} = strcmp(metaboliteCompartments(unitModels{k}), outside);
    metKeys{k} = repmat({''}, numel(mets), 1);
    for i = find(isExternal{k})'
        if isempty(opt.metaboliteKey)
            metKeys{k}{i} = stripCompartment(mets{i}, outside);
        else
            metKeys{k}{i} = opt.metaboliteKey(mets{i});
        end
        keyOfMet(mets{i}) = metKeys{k}{i};
    end
    examples{k} = mets(isExternal{k});
end

medium = resolveBounds(opt.poolMedium, keyOfMet, 'poolMedium');
[linkByRxn, linkByKey] = resolveLinkBounds(opt.linkBounds, keyOfMet);

%% the exchanges each unit reaches the outside world through
% Replaced by a transport into the pool, or by one shared exchange on the
% pool itself, so they are set aside here rather than carried over.
exchangeMask = cell(nUnits, 1);
for k = 1:nUnits
    exchangeMask{k} = externalExchanges(unitModels{k}, isExternal{k});
end

%% tag every unit and stack them into one model
untagged = opt.sharedMetabolites;
if ischar(untagged)
    untagged = {untagged};
end
[unitMetIds, poolMets] = communityMetaboliteIds(unitModels, organisms, poolsOf, ...
    isExternal, metKeys, untagged, pools, opt.mode, SHARED);
[community.model, unitRxns] = stackUnits(unitModels, organisms, unitMetIds, ...
    exchangeMask, opt.modelId);

%% replace each unit's exchanges with links into the pools it meets
links = repmat({cell(0, 1)}, nUnits, 1);
touched = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:nUnits
    for p = 1:numel(poolsOf{k})
        for i = find(isExternal{k})'
            mark = sprintf('%s|%s', poolsOf{k}{p}, metKeys{k}{i});
            if isKey(touched, mark)
                touched(mark) = union(touched(mark), organisms(k));
            else
                touched(mark) = organisms(k);
            end
        end
    end
end
if strcmp(opt.mode, POOLED)
    [community.model, links] = addPoolLinks(community.model, nUnits, organisms, ...
        poolsOf, isExternal, metKeys, unitMetIds, poolMets, linkByRxn, linkByKey);
    % Now that the transports exist, a bound that named one of them has been
    % applied; anything left over named nothing. Say so, because a misspelt
    % transport silently leaves its bounds wide open.
    if ~isempty(linkByRxn)
        missed = setdiff({linkByRxn.rxn}, vertcat(links{:}));
        if ~isempty(missed)
            fprintf(['buildCommunityModel: %d of %d transport(s) named in ' ...
                     'linkBounds do not exist and were ignored, for example ' ...
                     '%s. A transport is named TR_<metabolite>_<pool>__<unit>, ' ...
                     'with the metabolite as the pool spells it.\n'], ...
                numel(missed), numel(linkByRxn), ...
                strjoin(missed(1:min(5, numel(missed))), ', '));
        end
    end
end

%% exchange reactions on the pools that reach the environment
[community.model, communityExchanges] = addPoolExchanges(community.model, ...
    poolMets, environment, medium);

reportSharing(pools, touched, organisms, examples);

%% the community objective is the total anchor flux
taggedAnchors = cell(nUnits, 1);
for k = 1:nUnits
    taggedAnchors{k} = sprintf('%s__%s', anchors{k}, organisms{k});
end
missing = findRxnIDs(community.model, taggedAnchors) == 0;
if any(missing)
    error('buildCommunityModel:anchorNotKept', ...
        ['The anchor reaction of %s is not in the community; it may have ' ...
         'been an exchange reaction, which is not kept.'], ...
        strjoin(organisms(missing)', ', '));
end
if isempty(communityExchanges) && ~any(sum(community.model.S ~= 0, 1) == 1)
    error('buildCommunityModel:closed', ...
        ['The community has no boundary reaction at all, so nothing can ' ...
         'enter or leave it and no unit can grow. Name a pool in ' ...
         'environment, or keep a demand or sink reaction inside a unit.']);
end
community.model.c(:) = 0;
community.model.c(findRxnIDs(community.model, taggedAnchors)) = 1;

%% measure what each anchor can do before coupling throttles it
capacity = nan(nUnits, 1);
% Nothing to measure until something can be taken up. A community built
% before its medium is applied is starved by construction, and reporting
% every unit as inert would be noise rather than news.
fed = any(community.model.lb(sum(community.model.S ~= 0, 1) == 1) < 0);
if opt.check && ~fed
    fprintf(['buildCommunityModel: no boundary reaction permits uptake, so ' ...
             'the anchors were not measured. Apply a medium and rebuild to ' ...
             'see them.\n']);
elseif opt.check
    % The anchor flux at which coupling stops biting. Anything at or above
    % it needs no second look.
    clearAt = 0;
    if opt.couple && opt.couplingFactor > 0
        clearAt = max(0, (1 - opt.couplingThreshold) / opt.couplingFactor);
    end
    capacity = measureAnchors(community.model, taggedAnchors, clearAt);
    reportThrottling(organisms, taggedAnchors, capacity, opt.couple, ...
        opt.couplingFactor, opt.couplingThreshold);
end

%% tie each unit's reactions to its own anchor
if opt.couple
    for k = 1:nUnits
        own = setdiff(unitRxns{k}, taggedAnchors{k});
        own = [own; links{k}];                                    %#ok<AGROW>
        community.model = coupleRxnList2Rxn(community.model, own, ...
            taggedAnchors{k}, opt.couplingFactor, opt.couplingThreshold);
    end
    community.couplingFactor = opt.couplingFactor;
    community.couplingThreshold = opt.couplingThreshold;
else
    community.couplingFactor = [];
    community.couplingThreshold = [];
end

%% bookkeeping
community.organisms = organisms;
community.biomassReactions = taggedAnchors;
community.communityExchanges = communityExchanges;
community.pools = pools;
community.poolMetabolites = poolMets;
community.mode = opt.mode;
community.draftReactions = draftRxns;
community.databaseReactions = databaseRxns;

taggedDraft = cell(nUnits, 1);
taggedDatabase = cell(nUnits, 1);
for k = 1:nUnits
    taggedDraft{k} = tagAll(draftRxns{k}, organisms{k});
    taggedDatabase{k} = tagAll(databaseRxns{k}, organisms{k});
end
units = struct('name', organisms, 'anchor', taggedAnchors, ...
    'anchorSource', anchorSources, 'anchorCapacity', num2cell(capacity), ...
    'externalCompartment', externalOf, ...
    'rxns', cellfun(@(r, l) [r; l], unitRxns, links, 'UniformOutput', false), ...
    'draftRxns', taggedDraft, 'databaseRxns', taggedDatabase);
community.units = units;
end


function pools = resolvePools(given, organisms, sharedCompartment, mode, SHARED)
% Returns the pools as a validated struct array of name and units.
if isempty(given)
    pools = struct('name', sharedCompartment, 'units', {organisms});
    return
end
pools = given(:);
if ~all(isfield(pools, {'name', 'units'}))
    error('buildCommunityModel:badPools', ...
        'pools must be a struct array with fields name and units.');
end
names = {pools.name}';
if numel(unique(names)) ~= numel(names)
    error('buildCommunityModel:duplicatePool', 'Pool names must be unique.');
end
clash = intersect(names, organisms);
if ~isempty(clash)
    error('buildCommunityModel:poolIsUnit', ...
        ['A pool must be named by something that is not a unit, but %s is ' ...
         'both.'], strjoin(clash', ', '));
end
if strcmp(mode, SHARED) && numel(pools) > 1
    error('buildCommunityModel:tooManyPoolsForShared', ...
        ['mode ''%s'' shares the external compartment itself, so there can ' ...
         'only be one pool; got %d. Use mode ''pooled'' for units that must ' ...
         'meet in different places.'], SHARED, numel(pools));
end
for p = 1:numel(pools)
    members = pools(p).units(:);
    if isempty(members)
        error('buildCommunityModel:emptyPool', ...
            'Pool ''%s'' holds no units, so nothing could meet in it.', ...
            pools(p).name);
    end
    unknown = setdiff(members, organisms);
    if ~isempty(unknown)
        error('buildCommunityModel:unknownUnit', ...
            'Pool ''%s'' names %s, which is not a unit.', ...
            pools(p).name, strjoin(unknown', ', '));
    end
    pools(p).units = members;
end
homeless = setdiff(organisms, vertcat(pools.units));
if ~isempty(homeless)
    error('buildCommunityModel:unitInNoPool', ...
        ['Every unit must be in at least one pool, but %s is in none. A ' ...
         'unit in no pool cannot exchange anything and cannot grow.'], ...
        strjoin(homeless', ', '));
end
end


function environment = resolveEnvironment(given, pools)
% Returns the pools that carry exchange reactions.
%
% Defaulting to every pool would be wrong. A pool that exists only to link
% two units is an interface, and giving it exchanges would let those two
% dump metabolites straight into the environment, bypassing whichever pool
% was supposed to be the only way out.
names = {pools.name}';
if isempty(given)
    if numel(pools) > 1
        error('buildCommunityModel:environmentRequired', ...
            ['environment must be given when there is more than one pool: ' ...
             'a pool that merely links two units has to stay closed, or ' ...
             'they could bypass the rest of the system.']);
    end
    environment = names;
    return
end
if ischar(given)
    given = {given};
end
if ~iscellstr(given)
    error('buildCommunityModel:badEnvironment', ...
        ['environment must be a cell array of pool names, as in ' ...
         '''environment'', {''blood''}.']);
end
environment = given(:);
unknown = setdiff(environment, names);
if ~isempty(unknown)
    error('buildCommunityModel:unknownPool', ...
        'environment names %s, which is not a pool.', strjoin(unknown', ', '));
end
end


function databases = resolveDatabases(given, organisms)
% Returns the database of each unit, as a cell array with [] for none.
nUnits = numel(organisms);
databases = repmat({[]}, nUnits, 1);
if isempty(given)
    return
end

if isstruct(given) && all(isfield(given, {'organism', 'model'}))
    for d = 1:numel(given)
        at = find(strcmp(organisms, given(d).organism), 1);
        if isempty(at)
            error('buildCommunityModel:unknownDatabaseUnit', ...
                'databases names ''%s'', which is not a unit.', given(d).organism);
        end
        databases{at} = given(d).model;
    end
    return
end

if isstruct(given)
    % One model for every unit.
    if ~all(isfield(given, {'S', 'rxns', 'mets'}))
        error('buildCommunityModel:badDatabase', ...
            ['databases must be a COBRA model, a cell array of one per ' ...
             'unit, or a struct array with fields organism and model.']);
    end
    databases = repmat({given}, nUnits, 1);
    return
end

if numel(given) ~= nUnits
    error('buildCommunityModel:databaseCountMismatch', ...
        'Got %d databases for %d units. Pass one per unit, or one for all.', ...
        numel(given), nUnits);
end
databases = given(:);
end


function [merged, draftIds, databaseIds] = mergeWithDatabase(draft, database, organism)
% Returns a unit holding a whole database, bounded by its own model.
%
% The draft contributes its bounds rather than its reactions: where an
% identifier appears in both, the draft's bounds win, because they are what
% is known about this organism. Reactions the draft has and the database
% does not are kept, which a draft carved out of the database will rarely
% have beyond its biomass reaction.
extra = find(~ismember(draft.rxns, database.rxns));
merged = database;

if ~isempty(extra)
    merged = addReactionsToModel(merged, draft, extra);
end

% The draft's bounds win wherever the two agree.
[~, atDraft, atMerged] = intersect(draft.rxns, merged.rxns);
merged.lb(atMerged) = draft.lb(atDraft);
merged.ub(atMerged) = draft.ub(atDraft);

merged.c = zeros(numel(merged.rxns), 1);
objective = draft.rxns(draft.c ~= 0);
if isscalar(objective)
    at = findRxnIDs(merged, objective{1});
    if at > 0
        merged.c(at) = 1;
    end
end

draftIds = sort(draft.rxns(:));
databaseIds = sort(setdiff(database.rxns, draft.rxns));
fprintf(['buildCommunityModel: %s holds %d reactions from the database and ' ...
         '%d from the draft (%d of them the database does not have).\n'], ...
    organism, numel(database.rxns), numel(draft.rxns), numel(extra));
end


function model = addReactionsToModel(model, source, which)
% Appends the named reactions of another model, metabolites and all.
fresh = setdiff(source.mets(any(source.S(:, which) ~= 0, 2)), model.mets);
if ~isempty(fresh)
    model.mets = [model.mets; fresh(:)];
    model.S = [model.S; sparse(numel(fresh), numel(model.rxns))];
    if isfield(model, 'b')
        model.b = [model.b; zeros(numel(fresh), 1)];
    end
    if isfield(model, 'csense')
        model.csense = [model.csense(:); repmat('E', numel(fresh), 1)];
    end
    model = padMetField(model, source, fresh, 'metNames');
    model = padMetField(model, source, fresh, 'metFormulas');
end

columns = sparse(numel(model.mets), numel(which));
rows = findMetIDs(model, source.mets);
for j = 1:numel(which)
    nz = find(source.S(:, which(j)) ~= 0);
    columns(rows(nz), j) = source.S(nz, which(j));                %#ok<SPRIX>
end
model = appendReactions(model, source.rxns(which), ...
    fieldOrBlank(source, 'rxnNames', which), columns, ...
    source.lb(which), source.ub(which));
end


function model = padMetField(model, source, fresh, field)
% Extends an optional per-metabolite field to cover newly added rows.
if ~isfield(model, field)
    return
end
if isfield(source, field)
    at = findMetIDs(source, fresh);
    model.(field) = [model.(field)(:); source.(field)(at)];
else
    model.(field) = [model.(field)(:); repmat({''}, numel(fresh), 1)];
end
end


function values = fieldOrBlank(model, field, which)
% Returns a per-reaction field for the given reactions, or blanks.
if isfield(model, field)
    values = model.(field)(which);
else
    values = repmat({''}, numel(which), 1);
end
end


function [anchors, sources] = resolveAnchors(models, organisms, given)
% Returns each unit's anchor reaction and how it was decided.
%
% Guessing is deliberately avoided. Coupling every reaction of a unit to the
% wrong one would throttle or free the whole unit without any error, and a
% human model's anchor may be growth, maintenance, or a demand reaction with
% no telling name.
nUnits = numel(models);
anchors = cell(nUnits, 1);
sources = cell(nUnits, 1);

for k = 1:nUnits
    model = models{k};
    if ~isempty(given)
        if numel(given) ~= nUnits
            error('buildCommunityModel:anchorCountMismatch', ...
                'Got %d biomassReactions for %d models.', numel(given), nUnits);
        end
        if findRxnIDs(model, given{k}) == 0
            error('buildCommunityModel:unknownAnchor', ...
                'Unit ''%s'' has no reaction ''%s''.', organisms{k}, given{k});
        end
        anchors{k} = given{k};
        sources{k} = 'given';
        continue
    end

    objective = model.rxns(model.c ~= 0);
    if isscalar(objective)
        anchors{k} = objective{1};
        sources{k} = 'objective';
        continue
    end

    named = model.rxns(startsWith(lower(model.rxns), 'biomass'));
    if isscalar(named)
        anchors{k} = named{1};
        sources{k} = 'name';
        continue
    end

    candidates = model.rxns(startsWith(lower(model.rxns), 'biomass') | ...
        contains(lower(model.rxns), 'maintenance') | ...
        ismember(upper(model.rxns), {'ATPM', 'DM_ATP_C_', 'DM_ATP_C'}));
    hint = '';
    if ~isempty(candidates)
        hint = sprintf(' Candidates look like %s.', ...
            strjoin(candidates(1:min(5, numel(candidates)))', ', '));
    end
    error('buildCommunityModel:noAnchor', ...
        ['Cannot tell what unit ''%s'' is anchored on: it has %d objective ' ...
         'reaction(s) and %d named biomass*. Pass biomassReactions. The ' ...
         'anchor is what coupling gates on, and a unit that is not growing ' ...
         'is usually anchored on maintenance or demand rather than ' ...
         'biomass, so it has to be said rather than guessed.%s'], ...
        organisms{k}, numel(objective), numel(named), hint);
end
end


function externalOf = resolveExternalCompartments(models, organisms, given)
% Returns the compartment whose metabolites each unit may exchange.
%
% Worked out per unit rather than once for the community, because the units
% need not agree: one model may call it 'e' and another 'extracellular'.
nUnits = numel(models);
externalOf = cell(nUnits, 1);

if ischar(given) && ~isempty(given)
    externalOf(:) = {given};
    return
end
if iscell(given) && ~isempty(given)
    if numel(given) ~= nUnits
        error('buildCommunityModel:externalCountMismatch', ...
            'Got %d externalCompartment names for %d units.', numel(given), nUnits);
    end
    externalOf = given(:);
    return
end

for k = 1:nUnits
    externalOf{k} = findExternalCompartment(models{k}, organisms{k});
end
distinct = unique(externalOf);
if numel(distinct) > 1
    fprintf('buildCommunityModel: external compartments differ between units: %s.\n', ...
        strjoin(distinct', ', '));
end
end


function outside = findExternalCompartment(model, organism)
% Returns the compartment a model exchanges through.
%
% A compartment with a recognised name, else the one carrying the most
% single-metabolite reactions, else 'e', which is what nearly every
% reconstruction uses.
known = {'e', 'ex', 'extracellular', 'extraorganism', 'out', ...
         'extracellular space', 'extra organism', 'extra cellular', ...
         'extra-organism', 'external', 'external medium', 'env', 's', 'u'};

compartments = metaboliteCompartments(model);
boundary = find(sum(model.S ~= 0, 1) == 1);
counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
for r = boundary(:)'
    at = find(model.S(:, r) ~= 0, 1);
    name = compartments{at};
    if isKey(counts, name)
        counts(name) = counts(name) + 1;
    else
        counts(name) = 1;
    end
end

if counts.Count == 0
    outside = 'e';
    fprintf(['buildCommunityModel: %s has no exchange reaction, so which ' ...
             'compartment is external could not be told; assuming ''e''. ' ...
             'Pass externalCompartment to be sure.\n'], organism);
    return
end

names = keys(counts);
tally = cell2mat(values(counts));
recognised = ismember(lower(names), known);
if any(recognised)
    names = names(recognised);
    tally = tally(recognised);
end
[~, best] = max(tally);
outside = names{best};
end


function compartments = metaboliteCompartments(model)
% Returns the compartment of every metabolite of a model.
%
% The model's own annotation if it carries one, else the '[c]' suffix
% nearly every COBRA model uses, else the trailing '_c' of the other
% convention. The underscore is tried last and only when no metabolite
% carries a bracket, because a name like 'glc__D' would otherwise be read
% as a compartment 'D'.
if isfield(model, 'metComps') && isfield(model, 'comps') && ...
        numel(model.metComps) == numel(model.mets)
    compartments = reshape(model.comps(model.metComps), [], 1);
    return
end

bracketed = regexp(model.mets, '^.*\[(\w+)\]$', 'tokens', 'once');
if any(~cellfun(@isempty, bracketed))
    compartments = cellfun(@firstTokenOrEmpty, bracketed, 'UniformOutput', false);
    return
end

underscored = regexp(model.mets, '^.*_(\w+)$', 'tokens', 'once');
compartments = cellfun(@firstTokenOrEmpty, underscored, 'UniformOutput', false);
end


function value = firstTokenOrEmpty(token)
% Returns the single captured token, or '' when nothing matched.
if isempty(token)
    value = '';
else
    value = token{1};
end
end


function base = stripCompartment(met, compartment)
% Returns a metabolite identifier without its compartment suffix.
%
% Both conventions are in use - 'glc-D_e' and 'glc-D[e]' - and getting this
% wrong does not fail loudly. Each unit would simply trade through its own
% private pool, named after its own spelling, and the community would never
% cross-feed.
for suffix = {sprintf('[%s]', compartment), sprintf('_%s', compartment)}
    if endsWith(met, suffix{1})
        base = met(1:end - numel(suffix{1}));
        return
    end
end
base = met;
end


function resolved = resolveBounds(bounds, keyOfMet, what)
% Returns bounds keyed by pool identity rather than by metabolite.
%
% Resolved by looking each identifier up among the metabolites that actually
% exist rather than by stripping a suffix off it, so a metabolite named in a
% medium but present in no unit is reported instead of silently doing
% nothing.
resolved = struct('key', {}, 'lb', {}, 'ub', {});
if isempty(bounds)
    return
end
if ~all(isfield(bounds, {'met', 'lb', 'ub'}))
    error('buildCommunityModel:badBounds', ...
        '%s must be a struct array with fields met, lb and ub.', what);
end
unmatched = {};
for i = 1:numel(bounds)
    if ~isKey(keyOfMet, bounds(i).met)
        unmatched{end+1, 1} = bounds(i).met;                      %#ok<AGROW>
        continue
    end
    resolved(end+1) = struct('key', keyOfMet(bounds(i).met), ...
        'lb', bounds(i).lb, 'ub', bounds(i).ub);                  %#ok<AGROW>
end
if ~isempty(unmatched)
    fprintf(['buildCommunityModel: %d of %d metabolites named in %s are in ' ...
             'no unit and were ignored, for example %s.\n'], ...
        numel(unmatched), numel(bounds), what, ...
        strjoin(unmatched(1:min(5, numel(unmatched)))', ', '));
end
end


function [byRxn, byKey] = resolveLinkBounds(bounds, keyOfMet)
% Splits linkBounds into the transport-keyed and metabolite-keyed parts.
%
% A transport carries one metabolite between one unit and one pool, so
% bounding it is the per-unit statement a measured rate makes. Keying by
% metabolite instead sets the same bound on every unit's transport of it,
% which is a coarser thing to say but often the thing that is known.
byRxn = struct('rxn', {}, 'lb', {}, 'ub', {});
byKey = struct('key', {}, 'lb', {}, 'ub', {});
if isempty(bounds)
    return
end
if ~all(isfield(bounds, {'lb', 'ub'}))
    error('buildCommunityModel:badLinkBounds', ...
        'linkBounds must be a struct array with fields lb, ub and rxn or met.');
end
hasRxn = isfield(bounds, 'rxn');
hasMet = isfield(bounds, 'met');
if ~hasRxn && ~hasMet
    error('buildCommunityModel:badLinkBounds', ...
        ['linkBounds must name each link either by rxn, the transport ' ...
         'TR_<metabolite>_<pool>__<unit>, or by met, the metabolite as the ' ...
         'unit models spell it.']);
end

unmatched = {};
for i = 1:numel(bounds)
    if hasRxn && ~isempty(bounds(i).rxn)
        byRxn(end+1) = struct('rxn', bounds(i).rxn, ...
            'lb', bounds(i).lb, 'ub', bounds(i).ub);              %#ok<AGROW>
        continue
    end
    if hasMet && ~isempty(bounds(i).met)
        if ~isKey(keyOfMet, bounds(i).met)
            unmatched{end+1, 1} = bounds(i).met;                  %#ok<AGROW>
            continue
        end
        byKey(end+1) = struct('key', keyOfMet(bounds(i).met), ...
            'lb', bounds(i).lb, 'ub', bounds(i).ub);              %#ok<AGROW>
    end
end
if ~isempty(unmatched)
    fprintf(['buildCommunityModel: %d metabolite(s) named in linkBounds are ' ...
             'in no unit and were ignored, for example %s.\n'], ...
        numel(unmatched), strjoin(unmatched(1:min(5, numel(unmatched)))', ', '));
end
end


function mask = externalExchanges(model, isExternal)
% Returns the reactions that carry one metabolite out of the unit.
%
% Only the exchanges facing the external compartment. A demand or sink
% inside the cell is also a reaction with a single metabolite, and an anchor
% that only consumes looks exactly like one, so dropping every such reaction
% would delete it.
lone = sum(model.S ~= 0, 1) == 1;
mask = false(numel(model.rxns), 1);
for r = find(lone)
    at = find(model.S(:, r) ~= 0, 1);
    mask(r) = isExternal(at);
end
end


function [unitMetIds, poolMets] = communityMetaboliteIds(units, organisms, ...
    poolsOf, isExternal, metKeys, untagged, pools, mode, SHARED)
% Returns the community identifier of every metabolite of every unit.
%
% A metabolite is tagged by its unit unless it is named in
% sharedMetabolites, or unless the units share the external compartment, in
% which case the external metabolites are community-wide and found by what
% they are rather than by how one unit spells them.
nUnits = numel(units);
unitMetIds = cell(nUnits, 1);
canonical = containers.Map('KeyType', 'char', 'ValueType', 'char');

poolMets = struct('name', {}, 'mets', {}, 'keys', {});
for p = 1:numel(pools)
    poolMets(p).name = pools(p).name;
    poolMets(p).mets = cell(0, 1);
    poolMets(p).keys = cell(0, 1);
end
poolAt = containers.Map({pools.name}, num2cell(1:numel(pools)));

for k = 1:nUnits
    mets = units{k}.mets;
    ids = cell(numel(mets), 1);
    for i = 1:numel(mets)
        if strcmp(mode, SHARED) && isExternal{k}(i)
            pool = poolsOf{k}{1};
            mark = sprintf('%s|%s', pool, metKeys{k}{i});
            if isKey(canonical, mark)
                ids{i} = canonical(mark);
            else
                % One row for the whole community, under the first
                % spelling of it that turns up.
                ids{i} = mets{i};
                canonical(mark) = ids{i};
                p = poolAt(pool);
                poolMets(p).mets{end+1, 1} = ids{i};
                poolMets(p).keys{end+1, 1} = metKeys{k}{i};
            end
            continue
        end
        if ismember(mets{i}, untagged)
            ids{i} = mets{i};
            continue
        end
        ids{i} = tagMetabolite(mets{i}, organisms{k});
    end
    unitMetIds{k} = ids;
end

if strcmp(mode, SHARED)
    return
end

% In pooled mode the pool holds a metabolite of its own, one per thing any
% of its units can exchange.
for k = 1:nUnits
    for q = 1:numel(poolsOf{k})
        p = poolAt(poolsOf{k}{q});
        for i = find(isExternal{k})'
            key = metKeys{k}{i};
            if ismember(key, poolMets(p).keys)
                continue
            end
            poolMets(p).mets{end+1, 1} = sprintf('%s[%s]', key, poolMets(p).name);
            poolMets(p).keys{end+1, 1} = key;
        end
    end
end
end


function tagged = tagMetabolite(met, organism)
% Tags a metabolite, keeping any [compartment] suffix at the end.
%
% The tag goes on the compartment rather than after it, so that the unit
% keeps a private external compartment and COBRA's own parsing of a
% trailing [c] still works.
token = regexp(met, '^(.*)\[(\w+)\]$', 'tokens', 'once');
if isempty(token)
    tagged = sprintf('%s__%s', met, organism);
else
    tagged = sprintf('%s[%s__%s]', token{1}, token{2}, organism);
end
end


function tagged = tagAll(names, organism)
% Tags every identifier in a cell array by its unit.
if isempty(names)
    tagged = cell(0, 1);
    return
end
tagged = strcat(names(:), '__', organism);
end


function rules = tagRules(rules, genes, organism)
% Rewrites gene-protein-reaction rules to name the tagged genes.
%
% The rules are split into tokens once and every token that names a gene is
% rewritten. Looping over the genes instead, with one regexprep over the
% whole rule list per gene, is thousands of passes over thousands of rules:
% on Recon3D that is 95 s per unit against under a second here.
if isempty(rules) || isempty(genes)
    return
end
tagged = containers.Map(genes, strcat(genes(:), ['__' organism]));

% A token is matched whole, so a gene whose name is a prefix of another
% cannot capture it and nothing needs sorting by length. Operators like
% 'and' and 'or' are looked up, found to be no gene, and left alone.
[tokens, between] = regexp(rules, '[\w.\-]+', 'match', 'split');
for i = 1:numel(rules)
    if isempty(tokens{i})
        continue
    end
    found = isKey(tagged, tokens{i});
    if ~any(found)
        continue
    end
    tokens{i}(found) = values(tagged, tokens{i}(found));
    pieces = cell(1, 2 * numel(tokens{i}) + 1);
    pieces(1:2:end) = between{i};
    pieces(2:2:end) = tokens{i};
    rules{i} = [pieces{:}];
end
end


function [model, unitRxns] = stackUnits(units, organisms, unitMetIds, ...
    exchangeMask, modelId)
% Stacks the tagged units into one model over a merged metabolite universe.
nUnits = numel(units);
unitRxns = cell(nUnits, 1);
keep = cell(nUnits, 1);
for k = 1:nUnits
    keep{k} = ~exchangeMask{k};
    unitRxns{k} = strcat(units{k}.rxns(keep{k}), '__', organisms{k});
end

% The metabolites, in first-appearance order. Two units may contribute the
% same row - a shared metabolite, or a shared external compartment - so the
% universe is the union rather than the concatenation.
counts = cellfun(@numel, unitMetIds);
offsets = [0; cumsum(counts(:))];
allIds = vertcat(unitMetIds{:});
[model.mets, firstAt] = unique(allIds, 'stable');
model.mets = model.mets(:);
firstAt = firstAt(:);
nMets = numel(model.mets);
rowOf = containers.Map(model.mets, num2cell((1:nMets)'));

rowIdx = cell(nUnits, 1);
for k = 1:nUnits
    rowIdx{k} = reshape(cell2mat(values(rowOf, unitMetIds{k})), [], 1);
end

model.rxns = vertcat(unitRxns{:});
nRxns = numel(model.rxns);

rows = cell(nUnits, 1);
cols = cell(nUnits, 1);
vals = cell(nUnits, 1);
colAt = 0;
for k = 1:nUnits
    [ii, jj, vv] = find(units{k}.S(:, keep{k}));
    rows{k} = reshape(rowIdx{k}(ii), [], 1);
    cols{k} = reshape(jj, [], 1) + colAt;
    vals{k} = reshape(vv, [], 1);
    colAt = colAt + sum(keep{k});
end
model.S = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), nMets, nRxns);

lbParts = cell(nUnits, 1);
ubParts = cell(nUnits, 1);
for k = 1:nUnits
    lbParts{k} = units{k}.lb(keep{k});
    ubParts{k} = units{k}.ub(keep{k});
end
model.lb = vertcat(lbParts{:});
model.ub = vertcat(ubParts{:});

model.c = zeros(nRxns, 1);
model.b = zeros(nMets, 1);
model.csense = repmat('E', nMets, 1);
model.osenseStr = 'max';
model.modelID = modelId;
model.description = modelId;

model = carryOverRxns(model, units, keep, 'rxnNames');
model = carryOverRxns(model, units, keep, 'subSystems');
model = carryOverGrRules(model, units, organisms, keep);
% The first unit to contribute a row names it, which matters only where a
% row is shared and the units disagree.
owner = arrayfun(@(at) find(at > offsets, 1, 'last'), firstAt);
local = firstAt - offsets(owner);
model = carryOverMets(model, units, owner, local, 'metNames');
model = carryOverMets(model, units, owner, local, 'metFormulas');

genes = cell(nUnits, 1);
for k = 1:nUnits
    if isfield(units{k}, 'genes') && ~isempty(units{k}.genes)
        genes{k} = strcat(units{k}.genes(:), '__', organisms{k});
    else
        genes{k} = cell(0, 1);
    end
end
model.genes = vertcat(genes{:});
end


function model = carryOverRxns(model, units, keep, field)
% Copies an optional per-reaction field across the stack.
if ~all(cellfun(@(u) isfield(u, field), units))
    return
end
parts = cell(numel(units), 1);
for k = 1:numel(units)
    value = units{k}.(field);
    parts{k} = value(keep{k});
    parts{k} = parts{k}(:);
end
model.(field) = vertcat(parts{:});
end


function model = carryOverGrRules(model, units, organisms, keep)
% Copies the gene-protein-reaction rules across, renaming the tagged genes.
if ~all(cellfun(@(u) isfield(u, 'grRules'), units))
    return
end
parts = cell(numel(units), 1);
for k = 1:numel(units)
    rules = units{k}.grRules(keep{k});
    if isfield(units{k}, 'genes')
        rules = tagRules(rules, units{k}.genes, organisms{k});
    end
    parts{k} = rules(:);
end
model.grRules = vertcat(parts{:});
end


function model = carryOverMets(model, units, owner, local, field)
% Copies an optional per-metabolite field, the first contributor winning.
if ~all(cellfun(@(u) isfield(u, field), units))
    return
end
value = cell(numel(owner), 1);
for i = 1:numel(owner)
    source = units{owner(i)}.(field);
    value{i} = source{local(i)};
end
model.(field) = value;
end


function [model, links] = addPoolLinks(model, nUnits, organisms, poolsOf, ...
    isExternal, metKeys, unitMetIds, poolMets, linkByRxn, linkByKey)
% Replaces each unit's exchanges with transports into the pools it meets.
links = cell(nUnits, 1);

% The pool metabolites exist only as names so far, so they get their rows
% here, all at once.
fresh = setdiff(vertcat(poolMets.mets), model.mets);
model = addEmptyMetabolites(model, fresh);
rowOf = containers.Map(model.mets, num2cell((1:numel(model.mets))'));

% Every lookup below is by name and runs once per transport, of which there
% may be hundreds of thousands, so each one is a map rather than a scan.
metOfPoolKey = containers.Map('KeyType', 'char', 'ValueType', 'char');
for p = 1:numel(poolMets)
    for i = 1:numel(poolMets(p).mets)
        metOfPoolKey(sprintf('%s|%s', poolMets(p).name, poolMets(p).keys{i})) = ...
            poolMets(p).mets{i};
    end
end
boundOfRxn = boundsMap(linkByRxn, 'rxn');
boundOfKey = boundsMap(linkByKey, 'key');

nLinks = 0;
for k = 1:nUnits
    nLinks = nLinks + numel(poolsOf{k}) * sum(isExternal{k});
end
if nLinks == 0
    links = repmat({cell(0, 1)}, nUnits, 1);
    return
end

ids = cell(nLinks, 1);
names = cell(nLinks, 1);
lb = zeros(nLinks, 1);
ub = zeros(nLinks, 1);
rows = zeros(2 * nLinks, 1);
cols = zeros(2 * nLinks, 1);
vals = repmat([-1; 1], nLinks, 1);

at = 0;
for k = 1:nUnits
    mine = zeros(numel(poolsOf{k}) * sum(isExternal{k}), 1);
    taken = 0;
    for q = 1:numel(poolsOf{k})
        pool = poolsOf{k}{q};
        for i = find(isExternal{k})'
            key = metKeys{k}{i};
            at = at + 1;
            ids{at} = sprintf('TR_%s_%s__%s', key, pool, organisms{k});
            names{at} = sprintf('%s transport, %s to %s', key, organisms{k}, pool);
            [lb(at), ub(at)] = linkBoundsFor(boundOfRxn, boundOfKey, ids{at}, key);
            rows(2*at-1) = rowOf(unitMetIds{k}{i});
            rows(2*at) = rowOf(metOfPoolKey(sprintf('%s|%s', pool, key)));
            cols(2*at-1) = at;
            cols(2*at) = at;
            taken = taken + 1;
            mine(taken) = at;
        end
    end
    links{k} = ids(mine(1:taken));
end

columns = sparse(rows, cols, vals, numel(model.mets), nLinks);
model = appendReactions(model, ids, names, columns, lb, ub);
end


function lookup = boundsMap(bounds, field)
% Returns a struct array of bounds as a map from its key field to [lb ub].
lookup = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1:numel(bounds)
    lookup(bounds(i).(field)) = [bounds(i).lb, bounds(i).ub];
end
end


function [lb, ub] = linkBoundsFor(byRxn, byKey, name, key)
% Returns the bounds on what may cross between a unit and a pool.
%
% A transport named outright wins over a bound given for the metabolite,
% which is the broader statement.
lb = -1000;
ub = 1000;
if isKey(byKey, key)
    bounds = byKey(key);
    lb = bounds(1);
    ub = bounds(2);
end
if isKey(byRxn, name)
    bounds = byRxn(name);
    lb = bounds(1);
    ub = bounds(2);
end
end


function [model, communityExchanges] = addPoolExchanges(model, poolMets, ...
    environment, medium)
% Adds the exchange reactions through which the community meets the world.
open = find(ismember({poolMets.name}, environment));
nNew = sum(arrayfun(@(p) numel(poolMets(p).mets), open));
if nNew == 0
    communityExchanges = cell(0, 1);
    return
end

model = addEmptyMetabolites(model, ...
    setdiff(vertcat(poolMets(open).mets), model.mets));
rowOf = containers.Map(model.mets, num2cell((1:numel(model.mets))'));
boundOfKey = boundsMap(medium, 'key');

ids = cell(nNew, 1);
names = cell(nNew, 1);
lb = zeros(nNew, 1);
ub = zeros(nNew, 1);
rows = zeros(nNew, 1);

at = 0;
for p = open(:)'
    for i = 1:numel(poolMets(p).mets)
        poolMet = poolMets(p).mets{i};
        at = at + 1;
        ids{at} = sprintf('EX_%s', poolMet);
        names{at} = sprintf('%s community exchange', poolMet);
        % Export-only until the medium says otherwise: a community that
        % silently takes up everything is the more dangerous default.
        [lb(at), ub(at)] = mediumBoundsFor(boundOfKey, poolMets(p).keys{i});
        rows(at) = rowOf(poolMet);
    end
end

columns = sparse(rows, (1:nNew)', -ones(nNew, 1), numel(model.mets), nNew);
model = appendReactions(model, ids, names, columns, lb, ub);
communityExchanges = ids;
end


function [lb, ub] = mediumBoundsFor(boundOfKey, key)
% Returns the bounds on a pool exchange, which are export-only by default.
lb = 0;
ub = 1000;
if isKey(boundOfKey, key)
    bounds = boundOfKey(key);
    lb = bounds(1);
    ub = bounds(2);
end
end


function model = addEmptyMetabolites(model, mets)
% Adds metabolite rows that no reaction touches yet.
if isempty(mets)
    return
end
mets = mets(:);
model.mets = [model.mets; mets];
model.S = [model.S; sparse(numel(mets), numel(model.rxns))];
if isfield(model, 'b')
    model.b = [model.b; zeros(numel(mets), 1)];
end
if isfield(model, 'csense')
    model.csense = [model.csense(:); repmat('E', numel(mets), 1)];
end
for field = {'metNames', 'metFormulas'}
    if isfield(model, field{1})
        model.(field{1}) = [model.(field{1})(:); repmat({''}, numel(mets), 1)];
    end
end
end


function model = appendReactions(model, ids, names, columns, lb, ub)
% Appends a block of reactions, which is far cheaper than one at a time.
if isempty(ids)
    return
end
clash = intersect(ids, model.rxns);
if ~isempty(clash)
    error('buildCommunityModel:duplicateReaction', ...
        'Reaction %s already exists in the community.', clash{1});
end
model.rxns = [model.rxns; ids(:)];
model.S = [model.S, columns];
model.lb = [model.lb; lb(:)];
model.ub = [model.ub; ub(:)];
model.c = [model.c; zeros(numel(ids), 1)];
if isfield(model, 'rxnNames')
    model.rxnNames = [model.rxnNames(:); names(:)];
end
if isfield(model, 'subSystems')
    model.subSystems = [model.subSystems(:); blankLike(model.subSystems, numel(ids))];
end
if isfield(model, 'grRules')
    model.grRules = [model.grRules(:); repmat({''}, numel(ids), 1)];
end
end


function filler = blankLike(existing, n)
% Returns n blank entries shaped like the ones already in a field.
if ~isempty(existing) && iscell(existing) && iscell(existing{1})
    filler = repmat({{''}}, n, 1);
else
    filler = repmat({''}, n, 1);
end
end


function reportSharing(pools, touched, organisms, examples)
% Warns when the units of a pool have nothing in common.
%
% The failure this catches is silent and total. Two models can both report
% an 'e' compartment and still share nothing, because one writes 'glc__D[e]'
% and the other 'glc-D[e]': every unit then trades through a pool of its own
% and the community never cross-feeds, with no error anywhere. Compartment
% detection cannot fix that - the identifiers belong to different namespaces
% - so it is reported, and metaboliteKey is how to resolve it.
marks = keys(touched);
for p = 1:numel(pools)
    members = pools(p).units;
    if numel(members) < 2
        continue
    end
    prefix = sprintf('%s|', pools(p).name);
    mine = marks(startsWith(marks, prefix));
    shared = sum(cellfun(@(m) numel(touched(m)) > 1, mine));
    if shared > 0
        fprintf(['buildCommunityModel: pool ''%s'' has %d metabolites ' ...
                 'reached by more than one unit.\n'], pools(p).name, shared);
        continue
    end
    at = find(ismember(organisms, members), 1);
    sample = examples{at};
    fprintf(['buildCommunityModel: pool ''%s'' is shared by %d units but ' ...
             'not one metabolite in it is reached by more than one of them, ' ...
             'so nothing can cross-feed. The usual cause is identifiers ' ...
             'from different namespaces; pass metaboliteKey to say how they ' ...
             'correspond. In %s they look like %s.\n'], ...
        pools(p).name, numel(members), organisms{at}, ...
        strjoin(sample(1:min(3, numel(sample)))', ', '));
end
end


function capacity = measureAnchors(model, anchors, recheckBelow)
% Returns how much flux each unit's anchor can carry.
%
% Measured on the community rather than on each unit alone: a unit that
% cannot grow in isolation may grow perfectly well on what the others
% secrete, and reporting it as dead would be a false alarm. Measured before
% coupling, because coupling a unit to an anchor that cannot carry flux pins
% the unit at zero, which would make the reading agree with itself.
%
% One solve maximising every anchor at once lands on an arbitrary point of
% the optimal face, so a unit can read zero merely for having lost a tie to
% another that uses the same nutrient. That reading is a lower bound, which
% is enough for any unit already clear of recheckBelow; the rest are
% maximised individually, which costs one solve per unit but only for the
% units in question.
nUnits = numel(anchors);
capacity = nan(nUnits, 1);
at = findRxnIDs(model, anchors);

probe = model;
probe.c(:) = 0;
probe.c(at) = 1;
solution = optimizeCbModel(probe);
if solution.stat ~= 1
    return
end
flux = fluxOf(solution);
if ~isempty(flux)
    capacity = flux(at);
end

for k = 1:nUnits
    % Zero is always worth a second look, whatever the threshold: it is the
    % one reading that is reported, and it is also what a unit shows when it
    % merely lost a tie to another on the same nutrient.
    if ~isnan(capacity(k)) && capacity(k) > 0 && capacity(k) >= recheckBelow
        continue
    end
    probe.c(:) = 0;
    probe.c(at(k)) = 1;
    alone = optimizeCbModel(probe);
    if alone.stat == 1 && (isnan(capacity(k)) || alone.f > capacity(k))
        capacity(k) = alone.f;
    end
end
end


function flux = fluxOf(solution)
% Returns the flux vector of a COBRA solution, whichever field holds it.
if isfield(solution, 'v') && ~isempty(solution.v)
    flux = solution.v;
elseif isfield(solution, 'x') && ~isempty(solution.x)
    flux = solution.x;
else
    flux = [];
end
end


function reportThrottling(organisms, anchors, capacity, couple, factor, threshold)
% Warns where the chosen coupling would hold a unit below what it can do.
for k = 1:numel(anchors)
    if isnan(capacity(k))
        fprintf(['buildCommunityModel: %s cannot carry any flux through ' ...
                 '%s in this community.\n'], organisms{k}, anchors{k});
        continue
    end
    if abs(capacity(k)) <= 1e-9
        tail = '';
        if couple
            tail = ' and, coupled, is pinned at zero';
        end
        fprintf(['buildCommunityModel: %s has an anchor %s that cannot ' ...
                 'carry any flux even when maximised on its own, so the ' ...
                 'unit is inert%s.\n'], organisms{k}, anchors{k}, tail);
        continue
    end
    if ~couple
        continue
    end
    allowed = factor * capacity(k) + threshold;
    if allowed < 1
        fprintf(['buildCommunityModel: at its best %s runs %s at %.4g, so ' ...
                 'coupling caps its other reactions at %.4g. Raise ' ...
                 'couplingFactor or set couple to false.\n'], ...
            organisms{k}, anchors{k}, capacity(k), allowed);
    end
end
end
