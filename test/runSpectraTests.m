function results = runSpectraTests(solver)
% Runs every SPECTRA test and reports which passed
%
% USAGE:
%
%    results = runSpectraTests(solver)
%
% OPTIONAL INPUTS:
%    solver:             Solver to run the tests against (Default: 'gurobi',
%                        falling back to 'glpk' if it is not installed).
%
% OUTPUTS:
%    results:            Struct array with fields name, passed, message and
%                        seconds, one entry per test file.
%
% EXAMPLE:
%
%    ::
%
%        results = runSpectraTests;
%        results = runSpectraTests('glpk');
%
% NOTE:
%    Each test is a script of assertions. The first failure in a script
%    stops that script but not the run, so one broken thing does not hide
%    the state of everything else.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(fileparts(here),'src')));
addpath(here);

if ~exist('solver', 'var') || isempty(solver)
    solver = 'gurobi';
    if ~changeCobraSolver(solver, 'all', 0)
        fprintf('gurobi is not available; falling back to glpk\n');
        solver = 'glpk';
    end
end
changeCobraSolver(solver, 'all', 0);
fprintf('running the SPECTRA tests against %s\n\n', solver);

tests = dir(fullfile(here,'test*.m'));
results = repmat(struct('name','','passed',false,'message','','seconds',0), ...
    numel(tests), 1);

for k = 1:numel(tests)
    name = erase(tests(k).name, '.m');
    results(k).name = name;
    started = tic;
    try
        evalin('base', sprintf('clearvars; run(''%s'');', ...
            fullfile(here, tests(k).name)));
        results(k).passed = true;
    catch err
        results(k).passed = false;
        results(k).message = err.message;
    end
    results(k).seconds = toc(started);

    if results(k).passed
        fprintf('  PASS  %-22s %5.1f s\n', name, results(k).seconds);
    else
        fprintf('  FAIL  %-22s %5.1f s\n          %s\n', ...
            name, results(k).seconds, results(k).message);
    end
end

nPassed = sum([results.passed]);
fprintf('\n%d of %d test file(s) passed in %.0f s\n', ...
    nPassed, numel(results), sum([results.seconds]));

if nPassed < numel(results)
    error('runSpectraTests:failures', '%d test file(s) failed.', ...
        numel(results) - nPassed);
end
end
