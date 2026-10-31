function V = LP3( J, model )
%
% V = LP3( J, model )
%
% CPLEX implementation of LP-3 for input set J (see FASTCORE paper)
%
% BENCHMARK BUILD - SOLVER BACKEND ONLY.
% The single change from the shipped rFASTCORMICS file is the last stanza:
% the LP is handed to COBRA's solveCobraLP instead of cplexlp. The problem
% itself is untouched - same objective, same equalities, same bounds.
%
% See LP7_rFASTCORMICS.m in this folder for why, and for what it means when
% reporting the timings.
%
% (c) Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%     LCSB / LSRU, University of Luxembourg

[m,n] = size(model.S);

% objective
f = zeros(1,n);
f(J) = -1;

% equalities
Aeq = model.S;
beq = zeros(m,1);

% bounds
lb = model.lb;
ub = model.ub;

% ---- backend: solveCobraLP in place of cplexlp --------------------------
% cplexlp(f, [], [], Aeq, beq, lb, ub) minimises f'x subject to Aeq*x == beq.
% No inequality block here, so every row is an equality.
LPproblem.A      = Aeq;
LPproblem.b      = beq;
LPproblem.c      = f(:);
LPproblem.lb     = lb;
LPproblem.ub     = ub;
LPproblem.osense = 1;                     % minimise
LPproblem.csense = repmat('E', m, 1);

solution = solveCobraLP(LPproblem);
V = solution.full;

if isempty(V)
    error('LP3_rFASTCORMICS:noSolution', ...
          'LP3 returned no solution (stat=%d, origStat=%s).', ...
          solution.stat, num2str(solution.origStat));
end
end
