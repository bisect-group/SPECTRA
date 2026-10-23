% Tests for spectraME and spectraCCME: what comes back, and what happens
% when the core set cannot be satisfied.
%
% Run through runSpectraTests, which sets the solver and reports.

model = readCbModel('ecoli_core_model.mat');
consistentIDs = spectraCC(model, 1e-3, 'stoichiometry', 0);
blocked = setdiff(1:numel(model.rxns), consistentIDs);
consistent = removeRxns(model, model.rxns(blocked), 'metFlag', false);
core = findRxnIDs(consistent, {'PGI'});

%% an extraction keeps its core and shrinks the model
extracted = spectraME(consistent, core, 1e-3, 'stoichiometry', [], 1, ...
    'coreDirection', 'minNetLP', [], 0, {}, 0);
assert(ismember('PGI', extracted.rxns), 'the core reaction was not kept');
assert(numel(extracted.rxns) < numel(consistent.rxns), ...
    'the extraction returned no smaller than the model it was given');

%% an infeasible core raises instead of looping
% A blocked core reaction makes the direction LP infeasible. forward and
% reverse answer that with an all-NaN flux vector, and NaN fails the
% comparison that clears temp_core, so without a guard the while loop never
% terminates. This asserts the guard, and the timeout asserts that it is
% reached before any spinning.
caught = false;
started = tic;
try
    spectraME(model, findRxnIDs(model, {'PGI','FRUpts2'}), 1e-3);
catch err
    caught = strcmp(err.identifier, 'spectraME:infeasibleCore');
end
assert(caught, 'an infeasible core set should raise spectraME:infeasibleCore');
assert(toc(started) < 60, ...
    'the infeasible core took %.0f s to report, which suggests it is looping', ...
    toc(started));

%% spectraCCME drops a blocked core reaction rather than raising
[fromCCME, blockedCore] = spectraCCME(model, ...
    findRxnIDs(model, {'PGI','FRUpts2'}), 1e-3, 'stoichiometry', [], 1, ...
    'coreDirection', 'minNetLP', [], 0, 0);
assert(ismember('FRUpts2', blockedCore), ...
    'the blocked core reaction was not reported');
assert(~ismember('FRUpts2', fromCCME.rxns), ...
    'a reaction reported as blocked is in the returned model');
assert(ismember('PGI', fromCCME.rxns), ...
    'the core reaction that was not blocked should have been kept');

%% the formulations disagree, and minNetMILP is the smallest
sizes = struct();
for p = {'minNetLP','minNetMILP','growthOptim','tradeOff'}
    ex = spectraME(consistent, core, 1e-3, 'stoichiometry', [], 1, ...
        'coreDirection', p{1}, [], 0, {}, 0);
    sizes.(p{1}) = numel(ex.rxns);
end
assert(sizes.minNetMILP <= sizes.minNetLP, ...
    'minNetMILP minimises the reaction count itself, so it cannot exceed the L1 relaxation');
assert(sizes.tradeOff > sizes.minNetLP, ...
    'tradeOff with default weights rewards inclusion, so it should keep more');

%% a weight large enough to change the cheapest route swaps subnetworks
oxygen = {'EX_o2(e)','O2t','CYTBD'};
assert(all(ismember(oxygen, extracted.rxns)), ...
    'the unweighted extraction was expected to be aerobic');
weights = ones(numel(consistent.rxns), 1);
weights(findRxnIDs(consistent, oxygen)) = 100;
anaerobic = spectraME(consistent, core, 1e-3, 'stoichiometry', weights, 1, ...
    'coreDirection', 'minNetLP', [], 0, {}, 0);
assert(~any(ismember(oxygen, anaerobic.rxns)), ...
    'the penalised oxygen route was kept anyway');
assert(numel(anaerobic.rxns) < numel(extracted.rxns), ...
    'abandoning respiration should give a smaller model, not a larger one');

%% alternative solutions differ from each other
sols = spectraME(consistent, core, 1e-3, 'stoichiometry', [], 3, ...
    'coreDirection', 'minNetLP', [], 0, {}, 0);
assert(iscell(sols) && numel(sols) == 3, ...
    'nSol>1 should return a cell array of that many models');
assert(~isequal(sort(sols{1}.rxns), sort(sols{2}.rxns)), ...
    'two alternative solutions came back identical');
for i = 1:numel(sols)
    assert(ismember('PGI', sols{i}.rxns), ...
        'alternative solution %d dropped the core reaction', i);
end

%% remGene changes the genes and nothing else
kept   = spectraME(consistent, core, 1e-3, 'stoichiometry', [], 1, ...
    'coreDirection', 'minNetLP', [], 0, {}, 0);
pruned = spectraME(consistent, core, 1e-3, 'stoichiometry', [], 1, ...
    'coreDirection', 'minNetLP', [], 1, {}, 0);
assert(isequal(sort(kept.rxns), sort(pruned.rxns)), ...
    'removing unused genes changed the reactions too');
assert(numel(pruned.genes) < numel(kept.genes), ...
    'remGene=1 should leave fewer genes than remGene=0');

%% a seeded extraction repeats
a = spectraME(consistent, core, 1e-3, 'stoichiometry', [], 1, ...
    'coreDirection', 'minNetLP', [], 0, {}, 0);
b = spectraME(consistent, core, 1e-3, 'stoichiometry', [], 1, ...
    'coreDirection', 'minNetLP', [], 0, {}, 0);
assert(isequal(sort(a.rxns), sort(b.rxns)), ...
    'two runs with the same seed returned different models');

%% consistencyType carries through to the extraction
% On the n=1 Meneco toy nothing can carry flux at steady state, so the only
% core reaction is blocked and the extraction comes back empty. Under
% accumulation the same call keeps the whole network. Tutorial_topology
% works through this.
toyIDs = {'r1','r2','r3','r4','r5','r6'};
formulas = {' -> S[c]', 'S[c] -> b[c] + 1 a[c]', 'b[c] + a[c] -> d[c]', ...
            'a[c] -> c[c]', 'd[c] + c[c] -> T1[c]', 'T1[c] -> '};
toy = createModel();
for i = 1:numel(toyIDs)
    toy = addReaction(toy, toyIDs{i}, 'reactionFormula', formulas{i}, ...
        'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
end

[byTopology, blockedByTopology] = spectraCCME(toy, findRxnIDs(toy, {'r6'}), ...
    1e-4, 'topology', [], 1, 'coreDirection', 'minNetLP', [], 0, 0);
assert(numel(byTopology.rxns) == 6, ...
    'under accumulation the whole toy network should survive, not %d reactions', ...
    numel(byTopology.rxns));
assert(isempty(blockedByTopology), ...
    'under accumulation no core reaction should be blocked');

[bySteadyState, blockedBySteadyState] = spectraCCME(toy, ...
    findRxnIDs(toy, {'r6'}), 1e-4, 'stoichiometry', [], 1, ...
    'coreDirection', 'minNetLP', [], 0, 0);
assert(ismember('r6', blockedBySteadyState), ...
    'at steady state r6 should be reported blocked');
assert(isempty(bySteadyState.rxns), ...
    'at steady state the extraction should be empty, not %d reactions', ...
    numel(bySteadyState.rxns));

%% an unknown consistency type is refused
caught = false;
try
    spectraME(consistent, core, 1e-3, 'guesswork');
catch err
    caught = strcmp(err.identifier, 'spectraME:badConsType');
end
assert(caught, 'an unknown consistencyType should be refused');
