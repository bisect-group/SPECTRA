% Tests for metabolic tasks: building them, checking them, and finding the
% reactions they cannot do without.
%
% Run through runSpectraTests, which sets the solver and reports.

model = readCbModel('ecoli_core_model.mat');

HYDROLYSIS = {'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'};
COEFFS = [-1 -1 1 1 1];
AEROBIC = {'glc-D[e]','o2[e]','pi[e]','h2o[e]','h[e]'};
RESPIRED = {'co2[e]','h2o[e]','h[e]'};

%% taskEquation keeps what it is given
equation = taskEquation(HYDROLYSIS, COEFFS, [1 1000]);
assert(isequal(equation.mets, HYDROLYSIS(:)), 'metabolites were not preserved');
assert(isequal(equation.coeffs, COEFFS(:)), 'coefficients were not preserved');
assert(isequal(equation.bounds, [1 1000]), 'bounds were not preserved');

%% taskEquation refuses mismatched lengths
caught = false;
try
    taskEquation({'atp[c]','h2o[c]'}, [-1 -1 1]);
catch err
    caught = strcmp(err.identifier, 'taskEquation:lengthMismatch');
end
assert(caught, 'a coefficient vector of the wrong length was accepted');

%% bounds default to [0 1000] on both sides
% The direction is carried by the sign of the stoichiometric coefficient,
% not by the sign of the bound, so an input is not a negative-lower-bound
% exchange.
task = metabolicTask('t', 'inputs', AEROBIC, 'outputs', RESPIRED, ...
    'equations', equation);
assert(all(task.inputBounds(:,1) == 0) && all(task.inputBounds(:,2) == 1000), ...
    'input bounds did not default to [0 1000]');
assert(isequal(task.inputBounds, task.outputBounds(1,:) .* ones(numel(AEROBIC),1)), ...
    'inputs and outputs got different defaults');

%% a task with no positive lower bound anywhere is flagged
% It would be satisfied by the all-zero flux vector and pass trivially.
unforced = metabolicTask('loose', 'inputs', AEROBIC, 'outputs', RESPIRED);
assert(~taskIsForced(unforced), 'a task with no demand was reported as forced');
assert(taskIsForced(task), 'a task with a forced equation was reported as unforced');

%% taskConstraints closes the model's own exchanges and adds its own
[constrained, added] = taskConstraints(model, task);
assert(numel(added) == numel(AEROBIC) + numel(RESPIRED) + 1, ...
    'expected one reaction per input, output and equation; got %d', numel(added));
assert(all(ismember(added, constrained.rxns)), ...
    'taskConstraints reported reactions it did not add');
assert(numel(constrained.rxns) == numel(model.rxns) + numel(added), ...
    'the constrained model has the wrong number of reactions');

%% the model it was given is untouched
assert(numel(model.rxns) == 95, ...
    'taskConstraints modified the caller''s model');

%% a feasible task and a negative control behave
control = metabolicTask('no_free_atp', 'outputs', RESPIRED, ...
    'equations', equation, 'shouldFail', true);
results = checkTasks(model, [task; control], 0);
assert(results(1).feasible, 'the aerobic ATP task should be feasible');
assert(~results(2).feasible, 'ATP from nothing should be infeasible');
assert(all([results.ok]), 'both tasks should have behaved as expected');

%% a task naming an absent metabolite is reported, not raised
absent = metabolicTask('absent', 'inputs', {'unobtainium[e]'}, ...
    'outputs', RESPIRED, 'equations', equation);
result = checkTask(model, absent);
assert(~isempty(result.missing), 'the absent metabolite was not reported');
assert(~result.feasible, 'a task naming an absent metabolite cannot be feasible');

%% essential reactions are found, and are really essential
essential = essentialReactionsForTask(model, task);
assert(~isempty(essential), 'the task should have essential reactions');
assert(all(ismember(essential, model.rxns)), ...
    'an essential reaction is not in the model; scaffolding leaked into the result');

% Knocking one out must break the task.
knocked = model;
idx = findRxnIDs(knocked, essential{1});
knocked.lb(idx) = 0;
knocked.ub(idx) = 0;
assert(~checkTask(knocked, task).feasible, ...
    'removing the "essential" reaction %s left the task feasible', essential{1});

%% the union over a list skips the shouldFail controls
[core, perTask] = essentialReactionsForTasks(model, [task; control], {}, 0);
assert(numel(perTask) == 1, ...
    'the shouldFail control should have been skipped, leaving one task');
assert(isequal(sort(core), sort(essential)), ...
    'the union over one task differs from that task''s own set');
