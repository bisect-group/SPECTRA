% Tests for joining models into a community.
%
% Run through runSpectraTests, which sets the solver and reports.

%% two organisms, the second living off the first's waste
producer = toyOrganism('producer', 'glc', 'ac', '[%s]');
consumer = toyOrganism('consumer', 'ac', 'co2', '[%s]');
pair = {producer, consumer};
glucose = struct('met', {'glc[e]'}, 'lb', {-10}, 'ub', {1000});

%% joining keeps every organism separate
community = buildCommunityModel(pair, {'P', 'C'}, 'poolMedium', glucose);
assert(all(ismember({'UP__P', 'CONV__P', 'biomass__P'}, community.model.rxns)), ...
    'the producer''s reactions are not tagged by its name');
assert(all(ismember({'UP__C', 'CONV__C', 'biomass__C'}, community.model.rxns)), ...
    'the consumer''s reactions are not tagged by its name');
assert(isequal(community.biomassReactions, {'biomass__P'; 'biomass__C'}), ...
    'the anchors were not reported as the tagged reactions');
assert(strcmp(community.mode, 'pooled'), 'the default mode is not pooled');

%% a reaction that only consumes is not mistaken for an exchange
% biomass eats x[c] and makes nothing, so it has a single metabolite just
% as an exchange does. Dropping it would delete the anchor.
assert(ismember('biomass__P', community.model.rxns), ...
    'the anchor was dropped as though it were an exchange');

%% the organisms' own exchanges are replaced by links into the pool
assert(~any(startsWith(community.model.rxns, 'EX_glc[e]')), ...
    'a unit kept its own exchange reaction');
assert(ismember('TR_glc_u__P', community.model.rxns), ...
    'the producer has no transport for glucose');
assert(ismember('EX_glc[u]', community.communityExchanges), ...
    'the pool has no exchange for glucose');

%% cross-feeding through the shared pool
% Only glucose is on offer, and only the producer eats it. The consumer
% grows on the acetate the producer secretes into the pool.
grown = optimizeCbModel(changeObjective(community.model, 'biomass__C'));
assert(grown.stat == 1 && grown.f > 1e-6, ...
    'the consumer cannot grow on the producer''s waste (f = %g)', grown.f);

%% a closed community cannot grow
closed = buildCommunityModel(pair, {'P', 'C'}, 'check', false);
starved = optimizeCbModel(closed.model);
assert(abs(starved.f) < 1e-6, ...
    'a community given no medium grew at %g; the pool exchanges should be export-only', ...
    starved.f);

%% coupling shuts an organism down with its anchor
% With the consumer's anchor pinned at zero, coupling must stop the rest of
% its network too.
pinned = community.model;
pinned = changeRxnBounds(pinned, 'biomass__C', 0, 'b');
pinned = changeObjective(pinned, 'CONV__C');
held = optimizeCbModel(pinned);
assert(held.f <= 0.01 + 1e-6, ...
    'a unit with a dead anchor still ran its metabolism at %g', held.f);

%% without coupling a dead organism still runs
loose = buildCommunityModel(pair, {'P', 'C'}, 'poolMedium', glucose, ...
    'couple', false, 'check', false);
assert(isempty(loose.couplingFactor), ...
    'couplingFactor was reported even though couple was false');
free = loose.model;
free = changeRxnBounds(free, 'biomass__C', 0, 'b');
free = changeObjective(free, 'CONV__C');
ran = optimizeCbModel(free);
assert(ran.f > 0.01, ...
    'without coupling a dead unit should still run its metabolism, got %g', ran.f);

%% organism names are validated
assertFails(@() buildCommunityModel(pair, {'P', 'P'}), 'unique', ...
    'duplicate organism names were accepted');
assertFails(@() buildCommunityModel(pair, {'P__x', 'C'}), '__', ...
    'an organism name containing the separator was accepted');
assertFails(@() buildCommunityModel(pair, {'P'}), 'organism names', ...
    'a name for only one of two models was accepted');

%% both compartment spellings meet in the same pool
% Model collections differ on this, and getting it wrong does not fail
% loudly: each unit would trade through its own pool named after its own
% spelling, and the community would never cross-feed.
for mode = {'pooled', 'shared'}
    bracketed = toyOrganism('a', 'glc', 'ac', '[%s]');
    underscored = toyOrganism('b', 'glc', 'ac', '_%s');
    mixed = buildCommunityModel({bracketed, underscored}, {'A', 'B'}, ...
        'mode', mode{1}, 'poolMedium', glucose, 'check', false);
    keys = sort(mixed.poolMetabolites(1).keys);
    assert(isequal(keys, {'ac'; 'glc'}), ...
        '%s mode: the two spellings did not meet, pool holds %s', ...
        mode{1}, strjoin(keys', ', '));
end

%% shared mode needs no transports and agrees with pooled
pooled = buildCommunityModel(pair, {'P', 'C'}, 'poolMedium', glucose);
shared = buildCommunityModel(pair, {'P', 'C'}, 'mode', 'shared', ...
    'poolMedium', glucose);
assert(strcmp(shared.mode, 'shared'), 'the mode was not reported');
assert(~any(startsWith(shared.model.rxns, 'TR_')), ...
    'shared mode built transport reactions');
assert(any(startsWith(pooled.model.rxns, 'TR_')), ...
    'pooled mode built no transport reactions');
assert(numel(shared.model.rxns) < numel(pooled.model.rxns), ...
    'shared mode is not smaller: %d reactions against %d', ...
    numel(shared.model.rxns), numel(pooled.model.rxns));

onePooled = optimizeCbModel(changeObjective(pooled.model, 'biomass__C'));
oneShared = optimizeCbModel(changeObjective(shared.model, 'biomass__C'));
assert(abs(onePooled.f - oneShared.f) < 1e-6, ...
    'the two modes disagree on what the consumer can grow: %g against %g', ...
    onePooled.f, oneShared.f);

%% shared mode allows only one pool
% Its pool is the external compartment, and there is only one of those.
assertFails(@() buildCommunityModel(pair, {'P', 'C'}, 'mode', 'shared', ...
    'pools', struct('name', {'a', 'b'}, 'units', {{'P'}, {'C'}}), ...
    'environment', {'a'}), 'can only be one pool', ...
    'shared mode accepted two pools');

%% shared mode rejects linkBounds
% There is no per-unit reaction left to bound.
assertFails(@() buildCommunityModel(pair, {'P', 'C'}, 'mode', 'shared', ...
    'linkBounds', struct('met', {'glc[e]'}, 'lb', {0}, 'ub', {1})), ...
    'no such reactions to bound', 'shared mode accepted linkBounds');

%% three units, two on the blood and one reachable only through B
a = toyOrganism('a', 'glc', 'ac', '[%s]');
b = toyOrganism('b', 'glc', 'ac', '[%s]');
c = toyOrganism('c', 'ac', 'co2', '[%s]');
three = {a, b, c};
threePools = struct('name', {'blood', 'gut'}, ...
    'units', {{'A', 'B'}, {'B', 'C'}});

%% several pools must say which of them reach the environment
assertFails(@() buildCommunityModel(three, {'A', 'B', 'C'}, ...
    'pools', threePools), 'environment must be given', ...
    'two pools were accepted without saying which reaches the environment');

%% an interface pool gets no exchanges
% Giving the gut exchanges would let B and C bypass the blood entirely.
body = buildCommunityModel(three, {'A', 'B', 'C'}, 'pools', threePools, ...
    'environment', {'blood'}, 'poolMedium', glucose, 'check', false);
assert(~any(contains(body.communityExchanges, '[gut]')), ...
    'the interface pool was given exchange reactions');
assert(any(contains(body.communityExchanges, '[blood]')), ...
    'the open pool was given no exchange reactions');

%% a unit off the open pool lives through its neighbour
% C never touches the blood. It grows on the acetate B secretes into the gut.
fed = optimizeCbModel(changeObjective(body.model, 'biomass__C'));
assert(fed.stat == 1 && fed.f > 1e-6, ...
    'the unit off the open pool cannot grow (f = %g)', fed.f);

%% every unit must be in some pool
assertFails(@() buildCommunityModel(three, {'A', 'B', 'C'}, ...
    'pools', struct('name', {'blood'}, 'units', {{'A', 'B'}}), ...
    'environment', {'blood'}), 'is in none', ...
    'a unit belonging to no pool was accepted');

%% a pool may not be named after a unit
assertFails(@() buildCommunityModel(pair, {'P', 'C'}, ...
    'pools', struct('name', {'P'}, 'units', {{'P', 'C'}})), 'both', ...
    'a pool named after a unit was accepted');

%% linkBounds given by metabolite bound every unit's transport alike
capped = buildCommunityModel(pair, {'P', 'C'}, 'poolMedium', glucose, ...
    'linkBounds', struct('met', {'glc[e]'}, 'lb', {-2}, 'ub', {0}), ...
    'check', false);
at = findRxnIDs(capped.model, 'TR_glc_u__P');
assert(capped.model.lb(at) == -2 && capped.model.ub(at) == 0, ...
    'a metabolite-keyed linkBound did not reach the transport (got [%g %g])', ...
    capped.model.lb(at), capped.model.ub(at));

%% linkBounds can name one transport rather than a metabolite
% Which is the per-unit statement a measured rate makes.
oneUnit = buildCommunityModel(pair, {'P', 'C'}, 'poolMedium', glucose, ...
    'linkBounds', struct('rxn', {'TR_ac_u__C'}, 'lb', {-3}, 'ub', {4}), ...
    'check', false);
mine = findRxnIDs(oneUnit.model, 'TR_ac_u__C');
theirs = findRxnIDs(oneUnit.model, 'TR_ac_u__P');
assert(oneUnit.model.lb(mine) == -3 && oneUnit.model.ub(mine) == 4, ...
    'the named transport was not bounded');
assert(oneUnit.model.lb(theirs) == -1000, ...
    'bounding one unit''s transport changed another''s');

%% poolMedium constrains the community exchanges
tight = buildCommunityModel(pair, {'P', 'C'}, ...
    'poolMedium', struct('met', {'glc[e]'}, 'lb', {-4}, 'ub', {1000}), ...
    'check', false);
at = findRxnIDs(tight.model, 'EX_glc[u]');
assert(tight.model.lb(at) == -4, ...
    'the medium did not reach the pool exchange (lb = %g)', tight.model.lb(at));
total = optimizeCbModel(tight.model);
assert(abs(total.f - 4) < 1e-6, ...
    'the community grew at %g on an uptake capped at 4', total.f);

%% the anchor is reported with how it was decided
assert(strcmp(community.units(1).anchorSource, 'objective'), ...
    'the anchor taken from the objective was not reported as such');
named = buildCommunityModel(pair, {'P', 'C'}, 'poolMedium', glucose, ...
    'biomassReactions', {'biomass', 'biomass'}, 'check', false);
assert(strcmp(named.units(1).anchorSource, 'given'), ...
    'an anchor passed in was not reported as given');

%% an ambiguous anchor is refused rather than guessed
% Coupling a unit to the wrong reaction would throttle or free the whole
% unit with no error anywhere.
blank = toyOrganism('blank', 'glc', 'ac', '[%s]');
blank = changeObjective(blank, 'biomass', 0);
blank.rxns{strcmp(blank.rxns, 'biomass')} = 'grow';
assertFails(@() buildCommunityModel({blank, consumer}, {'X', 'C'}), ...
    'Cannot tell what unit', 'an ambiguous anchor was guessed at');

%% the external compartment is found per unit
% The units need not agree: one may call it 'e' and another something else.
odd = toyOrganism('odd', 'glc', 'ac', '[%s]', 'extracellular');
mixed = buildCommunityModel({producer, odd}, {'P', 'O'}, ...
    'poolMedium', glucose, 'check', false);
assert(strcmp(mixed.units(1).externalCompartment, 'e'), ...
    'the first unit''s external compartment was read as ''%s''', ...
    mixed.units(1).externalCompartment);
assert(strcmp(mixed.units(2).externalCompartment, 'extracellular'), ...
    'the second unit''s external compartment was read as ''%s''', ...
    mixed.units(2).externalCompartment);
assert(ismember('TR_glc_u__O', mixed.model.rxns), ...
    'the unit with its own compartment name got no transport');

%% externalCompartment can be said outright
told = buildCommunityModel({producer, odd}, {'P', 'O'}, ...
    'externalCompartment', {'e', 'extracellular'}, 'check', false);
assert(strcmp(told.units(2).externalCompartment, 'extracellular'), ...
    'a stated external compartment was not used');

%% metaboliteKey reconciles two namespaces
% Compartment detection cannot fix this; the identifiers differ outright.
one = toyOrganism('one', 'glc__D', 'ac', '[%s]');
two = toyOrganism('two', 'glc_D', 'ac', '[%s]');
apart = buildCommunityModel({one, two}, {'A', 'B'}, 'check', false);
assert(numel(apart.poolMetabolites(1).keys) == 3, ...
    'the two spellings of glucose were expected to stay apart, got %d keys', ...
    numel(apart.poolMetabolites(1).keys));

together = buildCommunityModel({one, two}, {'A', 'B'}, 'check', false, ...
    'metaboliteKey', @(met) strrep(regexprep(met, '\[\w+\]$', ''), '__', '_'));
assert(numel(together.poolMetabolites(1).keys) == 2, ...
    'metaboliteKey did not bring the two spellings together, got %d keys', ...
    numel(together.poolMetabolites(1).keys));

%% sharedMetabolites are left untagged
% A database with a single biomass metabolite drained by one community-wide
% exchange needs this.
communal = buildCommunityModel(pair, {'P', 'C'}, 'poolMedium', glucose, ...
    'sharedMetabolites', {'x[c]'}, 'check', false);
assert(ismember('x[c]', communal.model.mets), ...
    'the shared metabolite was tagged after all');
assert(~ismember('x[c__P]', communal.model.mets), ...
    'the shared metabolite was also kept per unit');

%% a database is folded in and its extra reactions recorded
draft = toyOrganism('draft', 'glc', 'ac', '[%s]');
draft = removeRxns(draft, {'CONV'});
database = toyOrganism('database', 'glc', 'ac', '[%s]');
filled = buildCommunityModel({draft, consumer}, {'D', 'C'}, ...
    'databases', {database, []}, 'poolMedium', glucose, 'check', false);
assert(ismember('CONV', filled.databaseReactions{1}), ...
    'the reaction only the database has was not recorded as a candidate');
assert(ismember('CONV__D', filled.model.rxns), ...
    'the database reaction was not folded into the unit');
assert(ismember('UP', filled.draftReactions{1}), ...
    'a reaction the draft already had was not recorded as its own');
assert(~ismember('UP', filled.databaseReactions{1}), ...
    'a reaction the draft already had was recorded as a candidate');
assert(isempty(filled.databaseReactions{2}), ...
    'a unit given no database was recorded as having candidates');

%% the draft's bounds win where the two agree
bounded = toyOrganism('bounded', 'glc', 'ac', '[%s]');
bounded = changeRxnBounds(bounded, 'UP', 7, 'u');
held = buildCommunityModel({bounded, consumer}, {'B', 'C'}, ...
    'databases', {database, []}, 'check', false);
at = findRxnIDs(held.model, 'UP__B');
assert(held.model.ub(at) == 7, ...
    'the draft''s bound was overwritten by the database''s (ub = %g)', ...
    held.model.ub(at));

%% a medium naming a metabolite no unit has says so
% Silently doing nothing would leave the community starved with no clue why.
absent = evalc(['buildCommunityModel(pair, {''P'', ''C''}, ''check'', false, ' ...
    '''poolMedium'', struct(''met'', {''nope[e]''}, ''lb'', {-1}, ''ub'', {1}));']);
assert(contains(absent, 'in no unit'), ...
    'a medium naming an absent metabolite was accepted in silence');

%% units sharing nothing at all are reported
% Two models can both report an 'e' compartment and still share nothing,
% because one writes 'glc__D' and the other 'glc_D'. Every unit then trades
% through a pool of its own and the community never cross-feeds, with no
% error anywhere, so it is reported.
strangerA = toyOrganism('strangerA', 'glc__D', 'ac__L', '[%s]');
strangerB = toyOrganism('strangerB', 'glc_D', 'ac_L', '[%s]');
quiet = evalc(['buildCommunityModel({strangerA, strangerB}, ' ...
    '{''A'', ''B''}, ''check'', false);']);
assert(contains(quiet, 'cross-feed'), ...
    'a pool whose units share no metabolite was not reported');

%% the anchor capacity is measured and reported
measured = buildCommunityModel(pair, {'P', 'C'}, 'poolMedium', glucose);
assert(~isnan(measured.units(1).anchorCapacity), ...
    'the anchor capacity was not measured');
assert(measured.units(1).anchorCapacity > 1e-6, ...
    'the producer was measured as unable to grow (%g)', ...
    measured.units(1).anchorCapacity);

%% a unit that merely lost a tie to another is not called inert
% One solve maximising every anchor at once lands on an arbitrary point of
% the optimal face, so two units on the same nutrient can leave one reading
% zero. Each such reading is taken again on its own.
rival = toyOrganism('rival', 'glc', 'ac', '[%s]');
tied = buildCommunityModel({producer, rival}, {'P', 'R'}, ...
    'poolMedium', glucose);
assert(all([tied.units.anchorCapacity] > 1e-6), ...
    'a unit competing for the only nutrient was measured as inert (%g, %g)', ...
    tied.units(1).anchorCapacity, tied.units(2).anchorCapacity);


function model = toyOrganism(name, eats, makes, pattern, outsideName)
% Returns a one-step organism: eats one metabolite, secretes another.
%
% pattern says how to spell a compartment, '[%s]' or '_%s', so that the
% tests can mix the two conventions in use. outsideName is what this unit
% calls its external compartment, which the units need not agree on.
if nargin < 5 || isempty(outsideName)
    outsideName = 'e';
end
outside = @(met) [met sprintf(pattern, outsideName)];
inside = ['x' sprintf(pattern, 'c')];

model = createModel();
model = addReaction(model, 'UP', 'reactionFormula', ...
    sprintf('%s -> %s', outside(eats), inside), ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, 'CONV', 'reactionFormula', ...
    sprintf('%s -> %s', inside, outside(makes)), ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
% Growth consumes without producing, so it has a single metabolite just as
% an exchange does. The builder must not mistake it for one.
model = addReaction(model, 'biomass', 'reactionFormula', ...
    sprintf('%s -> ', inside), ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, ['EX_' eats], 'reactionFormula', ...
    sprintf('%s <=> ', outside(eats)), ...
    'lowerBound', -1000, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, ['EX_' makes], 'reactionFormula', ...
    sprintf('%s <=> ', outside(makes)), ...
    'lowerBound', -1000, 'upperBound', 1000, 'printLevel', 0);
model = changeObjective(model, 'biomass');
model.modelID = name;
model.description = name;
end


function assertFails(call, fragment, message)
% Asserts that a call fails with a message holding the given fragment.
try
    call();
catch err
    assert(contains(err.message, fragment), ...
        '%s: failed with ''%s'', which does not mention ''%s''', ...
        message, err.message, fragment);
    return
end
error('%s: the call was expected to fail and did not', message);
end
