% Tests for the SPECTRA consistency check.
%
% Run through runSpectraTests, which sets the solver and reports.

model = readCbModel('ecoli_core_model.mat');

%% the textbook model has 87 consistent reactions of 95
ids = spectraCC(model, 1e-4);
assert(numel(ids) == 87, ...
    'expected 87 consistent reactions, got %d', numel(ids));

%% calls written against the earlier signatures still work
% detectionCutoff was added after seed, and seed after consistencyType,
% precisely so that these keep working.
ids2 = spectraCC(model, 1e-4);
ids3 = spectraCC(model, 1e-4, 'stoichiometry');
assert(isequal(ids2, ids3), ...
    'the two- and three-argument calls disagree');

%% a seed makes the result repeat
a = spectraCC(model, 1e-4, 'stoichiometry', 0);
b = spectraCC(model, 1e-4, 'stoichiometry', 0);
assert(isequal(a, b), 'the same seed gave two different answers');

%% the detection cutoff is honoured
% Passing the default explicitly must not change anything.
withDefault = spectraCC(model, 1e-4, 'stoichiometry', 0, 0.99*1e-4);
assert(isequal(withDefault, a), ...
    'passing the default detectionCutoff changed the answer');

%% topology never finds fewer reactions than stoichiometry
steady = spectraCC(model, 1e-4, 'stoichiometry', 0);
topo = spectraCC(model, 1e-4, 'topology', 0);
assert(numel(topo) >= numel(steady), ...
    'topology found %d consistent reactions against stoichiometry''s %d; it is the weaker requirement and cannot find fewer', ...
    numel(topo), numel(steady));

%% the steady state rejects a network that accumulation accepts
% Figure 2a of the Meneco article, with the coefficient n on the second
% reaction. The mass balances force v2*(n-2) == 0, so at steady state the
% whole network is dead unless n is exactly 2, while the accumulation
% condition lets the surplus pile up and accepts every n. Tutorial_topology
% works through this.
toyIDs = {'r1','r2','r3','r4','r5','r6'};
for n = 1:3
    formulas = {' -> S[c]', sprintf('S[c] -> b[c] + %d a[c]', n), ...
                'b[c] + a[c] -> d[c]', 'a[c] -> c[c]', ...
                'd[c] + c[c] -> T1[c]', 'T1[c] -> '};
    % The compartment suffixes matter: a metabolite named plainly b
    % collides with addReaction's own 'b' name-value argument.
    toy = createModel();
    for i = 1:numel(toyIDs)
        toy = addReaction(toy, toyIDs{i}, 'reactionFormula', formulas{i}, ...
            'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
    end

    nSteady = numel(spectraCC(toy, 1e-4, 'stoichiometry', 0));
    nAccum = numel(spectraCC(toy, 1e-4, 'topology', 0));

    expected = 0;
    if n == 2
        expected = 6;
    end
    assert(nSteady == expected, ...
        'at n=%d the steady state should accept %d of 6 reactions, not %d', ...
        n, expected, nSteady);
    assert(nAccum == 6, ...
        'at n=%d accumulation should accept all 6 reactions, not %d', n, nAccum);
end

%% accumulation rescues a missing consumer but not a missing producer
% The two networks differ only in the direction of r3. Where C is produced
% and never consumed, letting it accumulate unblocks r3. Where C is consumed
% and never produced, the balance reads -v3 >= 0, which with v3 >= 0 still
% forces v3 == 0, so nothing can rescue it.
ids8 = {'r1','r2','r3','r4','r5','r6','r7','r8'};
variants = { ...
    {' -> A[c]','A[c] -> B[c]','B[c] -> C[c]','B[c] -> D[c]', ...
     'D[c] -> ','A[c] -> E[c]','E[c] -> F[c]','F[c] -> D[c]'}, ...
    {' -> A[c]','A[c] -> B[c]','C[c] -> B[c]','B[c] -> D[c]', ...
     'D[c] -> ','A[c] -> E[c]','E[c] -> F[c]','F[c] -> D[c]'}};
expectedAccum = [8 7];

for k = 1:2
    net = createModel();
    for i = 1:numel(ids8)
        net = addReaction(net, ids8{i}, 'reactionFormula', variants{k}{i}, ...
            'lowerBound', 0, 'upperBound', 1000, 'printLevel', 0);
    end

    nSteady = numel(spectraCC(net, 1e-4, 'stoichiometry', 0));
    accum = net.rxns(spectraCC(net, 1e-4, 'topology', 0));

    assert(nSteady == 7, ...
        'variant %d: the steady state should block r3 alone, leaving 7, not %d', ...
        k, nSteady);
    assert(numel(accum) == expectedAccum(k), ...
        'variant %d: accumulation should accept %d of 8 reactions, not %d', ...
        k, expectedAccum(k), numel(accum));
end

%% an unrecognised consistencyType is refused
caught = false;
try
    spectraCC(model, 1e-4, 'stoichometry');   % deliberate typo
catch err
    caught = strcmp(err.identifier, 'spectraCC:badConsistencyType');
end
assert(caught, 'a misspelt consistencyType was not refused');

%% the returned model is consistent in its own right
consistent = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), ids)), ...
    'metFlag', false);
again = spectraCC(consistent, 1e-4, 'stoichiometry', 1);
assert(numel(again) == numel(consistent.rxns), ...
    'the consistent model is not a fixed point: %d of %d survive a second check', ...
    numel(again), numel(consistent.rxns));
