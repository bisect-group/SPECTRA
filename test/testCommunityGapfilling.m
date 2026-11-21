% Tests for gap-filling every unit of a community at once.
%
% Run through runSpectraTests, which sets the solver and reports.

fullA = gapfillOrganism('a', 'glc', 'ac');
fullB = gapfillOrganism('b', 'ac', 'co2');
draftA = removeRxns(fullA, {'UP'});     % cannot take up glucose
draftB = removeRxns(fullB, {'UP'});     % cannot take up acetate
glucose = struct('met', {'glc[e]'}, 'lb', {-10}, 'ub', {1000});
databases = struct('organism', {'a', 'b'}, 'model', {fullA, fullB});

community = buildCommunityModel({draftA, draftB}, {'a', 'b'}, ...
    'databases', databases, 'poolMedium', glucose, 'check', false);

%% neither draft can grow before the fill
% The community holds each unit's whole database, so the gap is in the
% drafts rather than in the joined model: it is the draft reactions the
% fill is free to keep, and everything else it must pay for.
drafts = buildCommunityModel({draftA, draftB}, {'a', 'b'}, ...
    'poolMedium', glucose, 'check', false);
before = optimizeCbModel(drafts.model);
assert(before.f < 1e-6, ...
    'the drafts already grow, so there is no gap to fill (f = %g)', before.f);
assert(isequal(community.databaseReactions{1}, {'UP'}), ...
    'the uptake reaction is not among the candidates to be filled in');

%% the fill takes what each unit needs from its own database
result = gapfillCommunity(community, 'seed', 1);
assert(isequal(result.added{1}, {'UP'}), ...
    'the first unit did not take the uptake reaction it was missing');
assert(isequal(result.added{2}, {'UP'}), ...
    'the second unit did not take the uptake reaction it was missing');

%% the filled community grows
after = optimizeCbModel(result.community.model);
assert(after.stat == 1 && after.f > 1e-6, ...
    'the filled community still cannot grow (f = %g)', after.f);

%% both anchors are required, so both units grow rather than only the cheaper one
anchors = findRxnIDs(result.community.model, result.community.biomassReactions);
assert(all(anchors > 0), 'the fill dropped an anchor it was told to keep');
each = optimizeCbModel(changeObjective(result.community.model, ...
    result.community.biomassReactions{2}));
assert(each.stat == 1 && each.f > 1e-6, ...
    ['the second unit cannot grow on its own, so the fill satisfied only ' ...
     'the cheaper anchor (f = %g)'], each.f);

%% the bookkeeping is narrowed to what survived
assert(all(ismember(result.community.units(1).rxns, result.community.model.rxns)), ...
    'a unit still lists reactions the filled model does not have');
assert(numel(result.core) == 2, ...
    'expected the two anchors as core, got %d reactions', numel(result.core));
assert(ischar(result.summary) && contains(result.summary, 'added'), ...
    'the summary does not say what was added');

%% the result splits back into standalone models
units = communityDecompose(result.community);
assert(ismember('UP', units(1).model.rxns), ...
    'the reaction the fill added did not survive the split');

%% keeping the draft is what separates a fill from an extraction
% A draft reaction that carries no flux under this medium is kept anyway.
% Shown with minNetLP, which decides by flux: under minNetMILP a free
% reaction gets no binary at all, so nothing is deciding whether to keep
% it and it survives either way.
idle = addReaction(draftA, 'IDLE', 'reactionFormula', 'y[c] <=> z[c]', ...
    'lowerBound', -1000, 'upperBound', 1000, 'printLevel', 0);
withIdle = buildCommunityModel({idle, draftB}, {'a', 'b'}, ...
    'databases', databases, 'poolMedium', glucose, 'check', false);

kept = gapfillCommunity(withIdle, 'seed', 1, 'probType', 'minNetLP', ...
    'consistencyCheck', false);
assert(ismember('IDLE__a', kept.community.model.rxns), ...
    'a draft reaction carrying no flux was thrown away by the fill');

dropped = gapfillCommunity(withIdle, 'seed', 1, 'probType', 'minNetLP', ...
    'keepDraft', false, 'consistencyCheck', false);
assert(~ismember('IDLE__a', dropped.community.model.rxns), ...
    'keepDraft false still kept a draft reaction the solution does not use');

%% under minNetMILP a reaction that costs nothing has no binary, so it stays
% Nothing is deciding whether to keep it, so nothing discards it. Worth
% stating, because it means keepDraft changes nothing under that
% formulation.
milp = gapfillCommunity(withIdle, 'seed', 1, 'keepDraft', false, ...
    'consistencyCheck', false);
assert(ismember('IDLE__a', milp.community.model.rxns), ...
    'a free reaction with no binary was dropped by the MILP');

%% a weight overrides what a reaction costs
% The database offers two ways to take up glucose, so the fill is free to
% pick either. Pricing one out is a statement it must heed, and it is how
% sequence evidence enters: a reaction the genome supports costs less.
twoWays = addReaction(fullA, 'UP2', 'reactionFormula', 'glc[e] -> x[c]', ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
choice = buildCommunityModel({draftA, draftB}, {'a', 'b'}, ...
    'databases', struct('organism', {'a', 'b'}, 'model', {twoWays, fullB}), ...
    'poolMedium', glucose, 'check', false);
priced = gapfillCommunity(choice, 'seed', 1, ...
    'weights', struct('rxn', {'UP__a'}, 'weight', {1000}));
assert(ismember('UP2', priced.added{1}) && ~ismember('UP', priced.added{1}), ...
    ['the fill took the reaction that was priced out of reach rather than ' ...
     'the free alternative; it added %s'], strjoin(priced.added{1}', ', '));

%% a weight naming nothing is refused rather than ignored
assertGapfillFails(@() gapfillCommunity(community, ...
    'weights', struct('rxn', {'nosuchreaction'}, 'weight', {1})), ...
    'cannot be weighted', 'a weight on a reaction that does not exist was accepted');

%% a core reaction the community does not have is refused
assertGapfillFails(@() gapfillCommunity(community, 'coreRxns', {'UP'}), ...
    'tagged by their unit', ...
    'an untagged core reaction was accepted without explanation');

%% tradeOff is refused, for the reason gapfillForGrowth refuses it
assertGapfillFails(@() gapfillCommunity(community, 'probType', 'tradeOff'), ...
    'minNetMILP', 'tradeOff was accepted for a growth gap-fill');

%% a community with no database says so rather than quietly extracting
bare = buildCommunityModel({fullA, fullB}, {'a', 'b'}, ...
    'poolMedium', glucose, 'check', false);
previous = warning('error', 'gapfillCommunity:noDatabase');
cleanup = onCleanup(@() warning(previous));
assertGapfillFails(@() gapfillCommunity(bare, 'seed', 1), 'database', ...
    'a community with nothing to draw on was filled without comment');
clear cleanup


function model = gapfillOrganism(name, eats, makes)
% Returns a one-step organism: eats one metabolite, secretes another.
model = createModel();
model = addReaction(model, 'UP', 'reactionFormula', ...
    sprintf('%s[e] -> x[c]', eats), ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, 'CONV', 'reactionFormula', ...
    sprintf('x[c] -> %s[e]', makes), ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, 'biomass', 'reactionFormula', 'x[c] -> ', ...
    'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, ['EX_' eats], 'reactionFormula', ...
    sprintf('%s[e] <=> ', eats), ...
    'lowerBound', -1000, 'upperBound', 1000, 'printLevel', 0);
model = addReaction(model, ['EX_' makes], 'reactionFormula', ...
    sprintf('%s[e] <=> ', makes), ...
    'lowerBound', -1000, 'upperBound', 1000, 'printLevel', 0);
model = changeObjective(model, 'biomass');
model.modelID = name;
model.description = name;
end


function assertGapfillFails(call, fragment, message)
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
