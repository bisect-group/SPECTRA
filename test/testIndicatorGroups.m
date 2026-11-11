% Tests for sharing one binary between several reactions in a MILP.
%
% Run through runSpectraTests, which sets the solver and reports.

model = readCbModel('ecoli_core_model.mat');
n = numel(model.rxns);

% minNet and tradeOff are given a model the way spectraME hands it to them:
% reoriented, with rev set.
prepared = model;
backwards = prepared.ub <= 0;
swap = prepared.ub(backwards);
prepared.S(:, backwards) = -prepared.S(:, backwards);
prepared.ub(backwards) = -1 * prepared.lb(backwards);
prepared.lb(backwards) = -1 * swap;
prepared.rev = ones(n, 1);
prepared.rev(prepared.lb >= 0) = 0;

growing = changeRxnBounds(prepared, prepared.rxns{prepared.c ~= 0}, 0.1, 'l');
noDirections = zeros(n, 1);
singletons = num2cell((1:n)');

%% resolveIndicatorGroups gives every free reaction its own binary by default
[groupOf, nGroups, ungoverned] = resolveIndicatorGroups({}, true(n, 1), n);
assert(nGroups == n, 'expected %d binaries by default, got %d', n, nGroups);
assert(isequal(groupOf, (1:n)'), 'the default grouping is not one per reaction');
assert(~any(ungoverned), 'a reaction was left ungoverned by the default grouping');

%% a directed reaction gets no binary
free = true(n, 1);
free(1) = false;
[groupOf, nGroups, ~] = resolveIndicatorGroups({}, free, n);
assert(nGroups == n - 1, 'a directed reaction was still given a binary');
assert(groupOf(1) == 0, 'a directed reaction was assigned to a group');

%% a free reaction in no group is reported as ungoverned
[~, nGroups, ungoverned] = resolveIndicatorGroups({[2 3]}, true(n, 1), n);
assert(nGroups == 1, 'expected one group, got %d', nGroups);
assert(sum(ungoverned) == n - 2, ...
    'expected %d ungoverned reactions, got %d', n - 2, sum(ungoverned));
assert(~any(ungoverned([2 3])), 'a grouped reaction was reported as ungoverned');

%% the grouping is validated
assertFails(@() resolveIndicatorGroups({[1 2], [2 3]}, true(n, 1), n), ...
    'both group', 'a reaction in two groups was accepted');
assertFails(@() resolveIndicatorGroups({[]}, true(n, 1), n), ...
    'is empty', 'an empty group was accepted');
assertFails(@() resolveIndicatorGroups({n + 1}, true(n, 1), n), ...
    'not a reaction index', 'a group named a reaction that does not exist');
directed = true(n, 1);
directed(5) = false;
assertFails(@() resolveIndicatorGroups({5}, directed, n), ...
    'given a direction', 'a group named a reaction already forced in');

%% minNet: singleton groups reproduce the default exactly
% The default formulation is the one-binary-per-reaction case of the
% general one, so saying it explicitly must change nothing.
[byDefault, ~, statDefault, selDefault] = minNet(growing, noDirections, ...
    ones(n, 1), 1e-4, 1, 'minNetMILP', 120);
[spelledOut, ~, statSpelled, selSpelled] = minNet(growing, noDirections, ...
    ones(n, 1), 1e-4, 1, 'minNetMILP', 120, [], {}, singletons);
assert(statDefault == 1 && statSpelled == 1, ...
    'the singleton comparison did not solve (stat %d and %d)', ...
    statDefault, statSpelled);
assert(isequal(byDefault, spelledOut), ...
    ['spelling out one group per reaction changed the answer from %d to %d ' ...
     'reactions'], sum(byDefault), sum(spelledOut));
assert(isequal(selDefault, selSpelled), ...
    'the two formulations switched on different binaries');

%% minNet: a shared binary keeps or drops its members together
active = find(abs(optimizeCbModel(growing).x) >= 1e-3)';
pairs = cell(0, 1);
for i = 1:2:numel(active) - 1
    pairs{end+1, 1} = [active(i), active(i+1)];                   %#ok<SAGROW>
end
[kept, flux, stat, selected] = minNet(growing, noDirections, ...
    ones(numel(pairs), 1), 1e-4, 1, 'minNetMILP', 120, [], {}, pairs);
assert(stat == 1, 'the grouped problem did not solve (stat %d)', stat);
assert(numel(selected) == numel(pairs), ...
    'expected one binary per group, got %d for %d groups', ...
    numel(selected), numel(pairs));
for g = 1:numel(pairs)
    if ~selected(g)
        assert(all(abs(flux(pairs{g})) <= 1e-9), ...
            'group %d is off but a reaction it governs carries flux', g);
    end
end

%% minNet: a reaction with no binary is kept whatever it does
% Nothing is deciding whether to keep it, so nothing should discard it.
first10 = num2cell((1:10)');
ungovernedKept = minNet(growing, noDirections, ones(10, 1), 1e-4, 1, ...
    'minNetMILP', 120, [], {}, first10);
assert(all(ungovernedKept(11:end)), ...
    'a free reaction outside every group was dropped');

%% minNet: weights are looked up by group, not by reaction
assertFails(@() minNet(growing, noDirections, ones(n, 1), 1e-4, 1, ...
    'minNetMILP', 60, [], {}, {[1 2]}), 'one entry per group', ...
    'per-reaction weights were accepted alongside groups');

%% tradeOff: singleton groups reproduce the default exactly
[tradeDefault, ~, tradeStat] = tradeOff(growing, noDirections, ones(n, 1), ...
    1e-4, 1, 120);
[tradeSpelled, ~, tradeStatSpelled] = tradeOff(growing, noDirections, ...
    ones(n, 1), 1e-4, 1, 120, [], {}, singletons);
assert(tradeStat == 1 && tradeStatSpelled == 1, ...
    'the tradeOff singleton comparison did not solve');
assert(isequal(tradeDefault, tradeSpelled), ...
    ['spelling out one group per reaction changed tradeOff''s answer from ' ...
     '%d to %d reactions'], sum(tradeDefault), sum(tradeSpelled));

%% tradeOff: a shared binary forces every member to carry flux
% This is the stronger statement minNet does not make: tradeOff requires
% every included reaction to carry at least tol, so switching a group on
% obliges all of its members to work.
[tKept, tFlux, tStat, tSelected] = tradeOff(growing, noDirections, ...
    ones(numel(pairs), 1), 1e-4, 1, 120, [], {}, pairs);
assert(tStat == 1, 'the grouped tradeOff did not solve (stat %d)', tStat);
for g = 1:numel(pairs)
    members = pairs{g};
    if tSelected(g)
        assert(all(abs(tFlux(members)) >= 1e-4 - 1e-9), ...
            ['group %d is on but a reaction it governs carries less than ' ...
             'tol; tradeOff must force every member to work'], g);
    else
        assert(all(abs(tFlux(members)) <= 1e-9), ...
            'group %d is off but a reaction it governs carries flux', g);
    end
end
assert(numel(tKept) == n, 'reacInd is not one entry per reaction');

%% tradeOff refuses a grouping it cannot satisfy rather than quietly relaxing
% Pairing unrelated reactions demands both carry flux whenever either does.
% With flux forced through biomass the cascade of partners cannot be met,
% and an infeasible status is the honest answer. minNet, which only permits
% its members to carry flux, accepts the same grouping.
arbitrary = cell(0, 1);
for i = 1:floor(n/2)
    arbitrary{end+1, 1} = [2*i-1, 2*i];                           %#ok<SAGROW>
end
[~, ~, impossible] = tradeOff(growing, noDirections, ...
    ones(numel(arbitrary), 1), 1e-4, 1, 60, [], {}, arbitrary);
assert(impossible ~= 1, ...
    'tradeOff accepted a grouping that cannot have every member carrying flux');

[~, ~, permitted] = minNet(growing, noDirections, ...
    ones(numel(arbitrary), 1), 1e-4, 1, 'minNetMILP', 60, [], {}, arbitrary);
assert(permitted == 1, ...
    ['minNet should accept the same grouping, since it only permits its ' ...
     'members to carry flux (stat %d)'], permitted);

%% tradeOff: weights are looked up by group, not by reaction
assertFails(@() tradeOff(growing, noDirections, ones(n, 1), 1e-4, 1, 60, ...
    [], {}, {[1 2]}), 'one entry per group', ...
    'per-reaction weights were accepted alongside groups in tradeOff');


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
