function V = LP7( J, model, epsilon )
%
% V = LP7( J, model, epsilon )
%
% CPLEX implementation of LP-7 for input set J (see FASTCORE paper)
%
% BENCHMARK BUILD - SOLVER BACKEND ONLY.
% The single change from the shipped rFASTCORMICS file is the last stanza:
% the LP is handed to COBRA's solveCobraLP instead of cplexlp. The problem
% itself is untouched - same objective, same equalities, same inequalities,
% same bounds, same variable ordering x = [v; z].
%
% Why: fasterCC's LP layer hard-wires cplexlp, the IBM CPLEX MATLAB API,
% which exists only up to CPLEX 12.10. That 2019-era MEX segfaults
% intermittently inside cplexlink12100.mexa64 under MATLAB 2023b (the oldest
% MATLAB on fir), taking the whole process down mid-run. Routing through
% solveCobraLP lets both FASTERCC and SPECTRA_CC run on Gurobi.
%
% Consequence for reporting: timings from this build are "FASTERCC's
% algorithm on a common Gurobi backend", NOT "FASTERCC as published". It also
% removes the per-LP API asymmetry that otherwise favours FASTERCC, since
% SPECTRA_CC already goes through solveCobraLP.
%
% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%     LCSB / LSRU, University of Luxembourg

nj = numel(J);
[m,n] = size(model.S);

% x = [v;z]

% objective
f = -[zeros(1,n), ones(1,nj)];

% equalities
Aeq = [model.S, sparse(m,nj)];
beq = zeros(m,1);

% inequalities
Ij = sparse(nj,n);
Ij(sub2ind(size(Ij),(1:nj)',J(:))) = -1;
Aineq = sparse([Ij, speye(nj)]);
bineq = zeros(nj,1);

% bounds
lb = [model.lb; zeros(nj,1)];
ub = [model.ub; ones(nj,1)*epsilon];

% ---- backend: solveCobraLP in place of cplexlp --------------------------
% cplexlp(f, Aineq, bineq, Aeq, beq, lb, ub) minimises f'x subject to
% Aineq*x <= bineq and Aeq*x == beq. The COBRA equivalent stacks the two
% constraint blocks and labels them with csense.
LPproblem.A      = [Aeq; Aineq];
LPproblem.b      = [beq; bineq];
LPproblem.c      = f(:);
LPproblem.lb     = lb;
LPproblem.ub     = ub;
LPproblem.osense = 1;                                    % minimise
LPproblem.csense = [repmat('E', m, 1); repmat('L', nj, 1)];

solution = solveCobraLP(LPproblem);
x = solution.full;

if isempty(x)
    % fasterCC relaxes bounds so that v = 0 is feasible, so LP7 should never
    % be infeasible. If it is, fail loudly rather than silently reporting
    % every reaction as blocked.
    error('LP7_rFASTCORMICS:noSolution', ...
          'LP7 returned no solution (stat=%d, origStat=%s).', ...
          solution.stat, num2str(solution.origStat));
end

V = x(1:n);
end
