% Tests for cutting a model down to the smallest network that still works.
%
% Run through runSpectraTests, which sets the solver and reports.

model = readCbModel('ecoli_core_model.mat');
M9 = {'EX_pi(e)','EX_h2o(e)','EX_h(e)','EX_nh4(e)','EX_o2(e)','EX_co2(e)'};
onCarbon = @(label, carbon) struct('label', label, ...
    'rxns', {[M9, {carbon}]}, ...
    'lb', {[-1000*ones(1, numel(M9)), -10]});
glucose = onCarbon('glucose', 'EX_glc(e)');
acetate = onCarbon('acetate', 'EX_ac(e)');

%% the result is a smaller model that still grows
half = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, 'seed', 0);
assert(numel(half.kept) < numel(model.rxns), ...
    'nothing was removed: %d of %d reactions', numel(half.kept), numel(model.rxns));
assert(numel(half.kept) + numel(half.removed) == numel(model.rxns), ...
    'the kept and removed reactions do not account for the model');
assert(isempty(half.unsatisfied), ...
    'the result does not meet: %s', strjoin(half.unsatisfied', ', '));
assert(half.growth.reached >= half.growth.required * (1 - 1e-6), ...
    'grew at %g against a requirement of %g', ...
    half.growth.reached, half.growth.required);

%% the model it returns really does grow on its own
grown = minimalReactome_grow(half.model, glucose);
assert(grown >= half.growth.required * (1 - 1e-6), ...
    'the returned model grew at %g, below the %g required', ...
    grown, half.growth.required);

%% preprocessing does not change the optimum
% Both of its passes are exact, so they can only change how many binaries
% the solver branches on.
with = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, ...
    'seed', 0, 'preprocess', true);
without = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, ...
    'seed', 0, 'preprocess', false);
assert(numel(with.kept) == numel(without.kept), ...
    ['preprocessing changed the answer from %d reactions to %d; both its ' ...
     'passes are exact and must not'], numel(without.kept), numel(with.kept));
assert(~isempty(with.essential), 'preprocessing proved nothing essential');
assert(isempty(without.essential), ...
    'reactions were reported essential with preprocessing off');

%% an essential reaction really is in the answer
assert(all(ismember(with.essential, with.kept)), ...
    'a reaction proved essential was not kept');

%% requiring more growth cannot need fewer reactions
lean = minimalReactome(model, 'media', glucose, 'growthFraction', 0.1, 'seed', 0);
full = minimalReactome(model, 'media', glucose, 'growthFraction', 1.0, 'seed', 0);
assert(numel(full.kept) >= numel(lean.kept), ...
    ['growing at the full rate needed %d reactions against %d for a tenth ' ...
     'of it, which cannot be'], numel(full.kept), numel(lean.kept));

%% two media are solved jointly, not one at a time
% The joint answer is at least as big as either alone and no bigger than
% their union. Were the conditions filled separately and unioned, it would
% be exactly that union; solving them together is what makes it smaller.
onGlucose = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, ...
    'seed', 0, 'preprocess', false);
onAcetate = minimalReactome(model, 'media', acetate, 'growthFraction', 0.5, ...
    'seed', 0, 'preprocess', false);
onBoth = minimalReactome(model, 'media', [glucose, acetate], ...
    'growthFraction', 0.5, 'seed', 0, 'preprocess', false);

assert(numel(onBoth.conditions) == 2, 'two media did not become two conditions');
assert(numel(onBoth.kept) >= max(numel(onGlucose.kept), numel(onAcetate.kept)), ...
    ['the joint answer has %d reactions, fewer than the %d one medium alone ' ...
     'needs'], numel(onBoth.kept), max(numel(onGlucose.kept), numel(onAcetate.kept)));
naive = union(onGlucose.kept, onAcetate.kept);
assert(numel(onBoth.kept) <= numel(naive), ...
    ['the joint answer has %d reactions against %d for the union of the two ' ...
     'solved separately; solving them together cannot do worse'], ...
    numel(onBoth.kept), numel(naive));

%% and it grows in both of them
assert(numel(onBoth.growth) == 2, 'growth was not reported for both media');
for i = 1:numel(onBoth.growth)
    assert(onBoth.growth(i).reached >= onBoth.growth(i).required * (1 - 1e-6), ...
        'the joint result grew at %g in %s, below the %g required', ...
        onBoth.growth(i).reached, onBoth.growth(i).label, ...
        onBoth.growth(i).required);
end

%% a product is required alongside growth, in the same flux distribution
secreting = minimalReactome(model, 'media', glucose, 'growthFraction', 0.1, ...
    'products', struct('id', 'EX_ac(e)', 'rate', 1.0), 'seed', 0);
assert(~isempty(secreting.production), 'no production was reported');
assert(secreting.production(1).carried >= 1.0 * (1 - 1e-6), ...
    'acetate was required at 1.0 but carries %g', secreting.production(1).carried);
assert(isempty(secreting.unsatisfied), 'the product requirement was not met');

%% a metabolite is resolved to the exchange the model already has
viaMetabolite = minimalReactome(model, 'media', glucose, ...
    'growthFraction', 0.1, 'products', struct('id', 'ac[e]', 'rate', 1.0), ...
    'seed', 0);
assert(strcmp(viaMetabolite.production(1).rxn, 'EX_ac(e)'), ...
    'the metabolite was resolved to ''%s''', viaMetabolite.production(1).rxn);

%% nothing is added to the network to measure a product
% A drain added here would make the answer the smallest subnetwork of a
% model the caller does not have, feasible only because of the reaction
% added to measure it.
assertFails(@() minimalReactome(model, 'media', glucose, ...
    'growthFraction', 0.1, 'products', struct('id', '13dpg[c]', 'rate', 1.0)), ...
    'no exchange reaction', 'an unexchanged metabolite was given a drain');

%% keepRxns are kept and made to carry flux
held = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, ...
    'keepRxns', {'SUCDi','FRD7'}, 'seed', 0, 'preprocess', false);
assert(all(ismember({'SUCDi','FRD7'}, held.kept)), ...
    'a reaction named in keepRxns was removed');
assert(numel(held.kept) >= numel(onGlucose.kept), ...
    'forcing two reactions in made the answer smaller');

%% keepRxns must name reactions the model has
assertFails(@() minimalReactome(model, 'media', glucose, ...
    'growthFraction', 0.5, 'keepRxns', {'NOPE'}), ...
    'does not have', 'keepRxns named a reaction the model lacks');

%% a weight of zero makes a reaction free to keep, not forced
freeToKeep = zeros(numel(model.rxns), 1);
freeToKeep(:) = 1;
freeToKeep(findRxnIDs(model, 'FRD7')) = 0;
weighted = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, ...
    'weights', freeToKeep, 'seed', 0, 'preprocess', false);
assert(numel(weighted.kept) <= numel(onGlucose.kept) + 1, ...
    'making one reaction free changed the answer by more than itself');

%% negative weights are refused
% They would pay the objective to keep a reaction, which is not what
% minimising a count means.
bad = ones(numel(model.rxns), 1);
bad(1) = -1;
assertFails(@() minimalReactome(model, 'media', glucose, ...
    'growthFraction', 0.5, 'weights', bad), 'negative weight', ...
    'a negative weight was accepted');

%% alternative solutions are distinct
three = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, ...
    'nSol', 3, 'seed', 0, 'preprocess', false);
assert(numel(three) > 1, 'nSol 3 returned %d solution(s)', numel(three));
signatures = arrayfun(@(r) strjoin(sort(r.kept)', '|'), three, ...
    'UniformOutput', false);
assert(numel(unique(signatures)) == numel(signatures), ...
    'the alternative solutions are not all distinct');
assert(numel(three(1).kept) <= numel(three(end).kept), ...
    'the first solution is not the smallest, but it is the minimal one');

%% a task becomes its own condition, with its own medium
HYDROLYSIS = {'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'};
AEROBIC = {'glc-D[e]','o2[e]','pi[e]','h2o[e]','h[e]'};
RESPIRED = {'co2[e]','h2o[e]','h[e]'};
equation = taskEquation(HYDROLYSIS, [-1 -1 1 1 1], [1 1000]);
atpTask = metabolicTask('aerobic_atp', 'inputs', AEROBIC, ...
    'outputs', RESPIRED, 'equations', equation);

withTask = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, ...
    'tasks', atpTask, 'seed', 0, 'preprocess', false);
assert(numel(withTask.conditions) == 2, ...
    'the task did not become a second condition');
assert(any(strcmp(withTask.conditions, 'aerobic_atp')), ...
    'the task was not named among the conditions');
assert(isempty(withTask.unsatisfied), ...
    'the result does not meet: %s', strjoin(withTask.unsatisfied', ', '));

%% and the result really can still perform it
outcome = checkTask(withTask.model, atpTask);
assert(outcome.ok, 'the minimal reactome cannot perform the task it was given');

%% requiring a task cannot need fewer reactions than not requiring it
assert(numel(withTask.kept) >= numel(onGlucose.kept), ...
    ['requiring a task alongside growth needed %d reactions against %d for ' ...
     'growth alone'], numel(withTask.kept), numel(onGlucose.kept));

%% the task's own scaffolding is never reported as network
assert(all(ismember(withTask.kept, model.rxns)), ...
    'a reaction the task added was reported as part of the network');

%% a task the model cannot perform is refused
impossible = metabolicTask('impossible', 'inputs', {'unobtainium[e]'}, ...
    'outputs', RESPIRED);
assertFails(@() minimalReactome(model, 'media', glucose, ...
    'growthFraction', 0.5, 'tasks', impossible), 'cannot perform task', ...
    'a task the model itself fails was required of a subset');

%% a negative control is skipped rather than required
control = metabolicTask('no_free_atp', 'outputs', RESPIRED, ...
    'equations', equation, 'shouldFail', true);
skipped = minimalReactome(model, 'media', glucose, 'growthFraction', 0.5, ...
    'tasks', control, 'seed', 0, 'preprocess', false);
assert(~any(strcmp(skipped.conditions, 'no_free_atp')), ...
    'a shouldFail task was required of the minimal reactome');

%% asking for nothing is refused
% A medium says what may be taken up, not what has to happen.
assertFails(@() minimalReactome(model, 'media', glucose, ...
    'growthFraction', 0), 'Nothing was asked', ...
    'a minimal reactome was built with no requirement');

%% the two growth requirements are mutually exclusive
assertFails(@() minimalReactome(model, 'media', glucose, ...
    'minGrowth', 0.1, 'growthFraction', 0.5), 'not both', ...
    'minGrowth and growthFraction were accepted together');

%% an unreachable growth rate is refused rather than blamed on the solve
assertFails(@() minimalReactome(model, 'media', glucose, 'minGrowth', 1e6), ...
    'cannot meet the requirement', ...
    'a growth rate the model itself cannot reach was accepted');

%% a medium naming an unknown exchange is refused
assertFails(@() minimalReactome(model, ...
    'media', struct('label','bad','rxns',{{'EX_nope(e)'}},'lb',{-10}), ...
    'growthFraction', 0.5), 'does not have', ...
    'a medium named an exchange the model lacks');

%% two media sharing a label are refused
% They would become one copy of the network rather than two.
same = [onCarbon('x','EX_glc(e)'), onCarbon('x','EX_ac(e)')];
assertFails(@() minimalReactome(model, 'media', same, 'growthFraction', 0.5), ...
    'share a label', 'two media with the same label were accepted');


function grown = minimalReactome_grow(model, medium)
% Returns what a model grows at in a medium, from its own bounds out.
isExc = findExcRxns(model);
model.lb(isExc) = 0;
present = medium.rxns(ismember(medium.rxns, model.rxns));
bounds = medium.lb(ismember(medium.rxns, model.rxns));
model.lb(findRxnIDs(model, present)) = bounds(:);
solution = optimizeCbModel(model);
grown = 0;
if solution.stat == 1
    grown = solution.f;
end
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
