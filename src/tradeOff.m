function [reacInd,x,stat,selected] = tradeOff(model,direction,weights,tol,steadyState,solveTime,x0,prevSols,indicatorGroups,inclusionCutoff)
% USAGE:
%   [reacInd,x,stat,selected] = tradeOff(model,direction,weights,tol,steadystate,solveTime,x0,prevSols,indicatorGroups,inclusionCutoff)
%
% INPUTS:
%     model:       COBRA model structure.
%     direction:   A vector of size equal to number of reactions in the
%                  model. The unique elements in this vector has to be -1,0
%                  and 1 defining the flux directionality info.
%     weights:     weights for non-core reactions. Positive weighted reactions have more
%                  chance to get included and negative weighted reactions
%                  have less chance to get included
%     tol:         tolerance level (minimum absolute flux that has to be carried
%                  by a reaction for it to be defined as consistent)
%     steadyState: Boolean value indicating whether to assume steady state
%                  condition (S.v = 0) or accumulation condition (S.v >= 0)
%
% OPTIONAL INPUTS:
%     solveTime:   Upper limit of solving time (Default-7200)
%     x0:          solution vector to initialize with for MILP problem.
%                  This can be obtained from the LPforward and LPreverse
%     prevSols:    A cell of previously obtained solutions that needs to be
%                  excluded in the current solution. Each entry holds
%                  reaction indices, or group indices when indicatorGroups
%                  is given.
%     indicatorGroups: Cell array of index vectors, one per group, giving the
%                  reactions that share a single binary. They are then
%                  included or excluded together and score the group's
%                  weight once however many of them there are, and weights
%                  must hold one entry per group rather than one per
%                  reaction. A free reaction in no group gets no binary and
%                  is always kept. (Default: one binary per free reaction)
%
%                  Note that a shared binary is a stronger statement here
%                  than it is in minNet. tradeOff requires every included
%                  reaction to carry at least tol, so switching a group on
%                  forces *every* reaction in it to carry flux, whereas
%                  minNet only permits them to. Group reactions together
%                  here only when you mean that all of them must work.
%     inclusionCutoff: The absolute flux at which a reaction counts as used.
%                  (Default: tol*1e-7)
%
% OUTPUTS:
%     reacInd: Reaction IDs corresponding to reactions that has to be
%              present in the final extracted model
%     x:       Solution returned after optimization
%     stat:    Status of optimization
%     selected: Logical over the groups, saying which binaries the solver
%              switched on. With one binary per free reaction this is a
%              logical over those reactions, in model order.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras
%
% .. See also: minNet, resolveIndicatorGroups


if ~exist('solveTime', 'var') || isempty(solveTime)
   solveTime=7200;
end
if ~exist('indicatorGroups', 'var')
    indicatorGroups = {};
end
if ~exist('inclusionCutoff', 'var') || isempty(inclusionCutoff)
    inclusionCutoff = tol*1e-7;
end

[m,n] = size(model.S);
dir0 = direction==0; % boolean value indicating if the reaction is non-core (1) or not (0)

% Which binary governs which reaction. With no groups every free reaction
% gets its own, which is the usual formulation; with groups the members of
% one group share a binary and a free reaction in no group gets none and is
% always kept.
[groupOf,nGroups,ungoverned] = resolveIndicatorGroups(indicatorGroups,dir0,n);
governed = groupOf > 0;
govIrr = find(governed & ~model.rev); % governed irreversible reactions
govRev = find(governed & model.rev);  % governed reversible reactions
nIrr = numel(govIrr);
n_rev = numel(govRev);

if isempty(indicatorGroups)
    groupWeights = weights(dir0);
else
    if numel(weights) ~= nGroups
        error('tradeOff:groupWeightCount', ...
            ['With indicatorGroups, weights must hold one entry per group, ' ...
             'but got %d weights for %d groups. What a group scores is ' ...
             'looked up by the group, not by its members.'], ...
            numel(weights), nGroups);
    end
    groupWeights = weights(:);
end

% the decision variables are v,z,a,b. v defines the flux values. z define
% the binary integer for the noncore reactions. a and b are for
% reverse and forward direction of the noncore reversible reactions.

% objective
f = [zeros(n,1);groupWeights;zeros(n_rev*2,1)];

% stoichiometry or topology constraints
Aeq1 = [model.S, sparse(m,nGroups+(n_rev*2))];
beq1 = zeros(m,1);

if steadyState
    csenseeq1 = repmat('E',m,1); % equality (For consistency based model extraction)
else
    csenseeq1 = repmat('G',m,1); % greater than (For topology based model extraction)
end

% Picks out a governed reaction's own flux on the left, and its group's
% binary on the right. With one group per reaction these are the identity
% blocks the earlier formulation wrote out directly.
pickIrr = sparse(1:nIrr, govIrr, 1, nIrr, n);
shareIrr = sparse(1:nIrr, groupOf(govIrr), 1, nIrr, nGroups);
pickRev = sparse(1:n_rev, govRev, 1, n_rev, n);
shareRev = sparse(1:n_rev, groupOf(govRev), 1, n_rev, nGroups);

% relation between the flux of non-core irreversible reactions and the
% corresponding binary variables
Aineq1 = [-1*pickIrr, tol*shareIrr, sparse(nIrr,n_rev*2)];
bineq1 = zeros(nIrr,1);
csenseineq1 = repmat('L',nIrr,1); % lesser than

Aineq2 = [pickIrr, -1*spdiags(model.ub(govIrr),0,nIrr,nIrr)*shareIrr, ...
    sparse(nIrr,n_rev*2)];
bineq2 = zeros(nIrr,1);
csenseineq2 = repmat('L',nIrr,1); % lesser than

% constraints to ignore the previously obtained solutions
Aineq3=[];bineq3=[];csenseineq3=[];
if exist('prevSols', 'var') && ~isempty(prevSols)
    for id=1:numel(prevSols)
        if isempty(indicatorGroups)
            temp = find(dir0);
            temp = ismember(temp,prevSols{id});
        else
            temp = ismember(1:nGroups,prevSols{id});
        end
        temp =temp(:)';
        Aineq3 = [Aineq3;[zeros(1,n),temp,zeros(1,n_rev*2)]];     %#ok<AGROW>
        bineq3 = [bineq3;sum(temp)-1];                            %#ok<AGROW>
        csenseineq3 = [csenseineq3;'L'];                          %#ok<AGROW>
    end
end

% constraints for reversible non-core reactions to carry flux only in
% one direction. Either in forward direction of reverse direction.
% a + b = z
Aeq2 = [sparse(n_rev,n), -1*shareRev, speye(n_rev), speye(n_rev)];
beq2 = zeros(n_rev,1);
csenseeq2 = repmat('E',n_rev,1);

% constraint linking the binary variable, a and the flux through the
% reversible non-core reactions
Aineq4 = [pickRev, sparse(n_rev,nGroups), -1*spdiag(tol*ones(n_rev,1)), ...
    -1*spdiag(model.lb(govRev))];
bineq4 = zeros(n_rev,1);
csenseineq4 = repmat('G',n_rev,1);

% constraint linking the binary variable, b and the flux through the
% reversible non-core reactions
Aineq5 = [pickRev, sparse(n_rev,nGroups), -1*spdiag(model.ub(govRev)), ...
    spdiag(tol*ones(n_rev,1))];
bineq5 = zeros(n_rev,1);
csenseineq5 = repmat('L',n_rev,1);

% bounds
lb = model.lb;
lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
lb = [lb;zeros(nGroups+(2*n_rev),1)];
ub = model.ub;
ub(direction==-1)=-tol*ones(sum(direction==-1),1);
ub = [ub;ones(nGroups+(2*n_rev),1)];

% Set up MILP problem
MILPproblem.A=[Aeq1;Aeq2;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5];
MILPproblem.b=[beq1;beq2;bineq1;bineq2;bineq3;bineq4;bineq5];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=-1; % maximize
MILPproblem.vartype = [repmat('C',n,1);repmat('B',nGroups+(2*n_rev),1)];
MILPproblem.csense = [csenseeq1; csenseeq2; csenseineq1; csenseineq2; csenseineq3; csenseineq4; csenseineq5];
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
fprintf(['tradeOff: MILP feasTol set to 1e-9. This is a global COBRA ' ...
         'setting and stays in effect after this function returns.\n']);
if exist('x0', 'var')&&~isempty(x0)
    warm = false(nGroups,1);
    on = false(n,1);
    on(1:min(n,numel(x0))) = abs(x0(1:min(n,numel(x0))))>1e-12;
    warm(groupOf(governed & on)) = true;
    MILPproblem.x0  = [x0;double(warm);x0(govRev)>1e-12;x0(govRev)<-1e-12];
end
solution = solveCobraMILP(MILPproblem,'timeLimit', solveTime);
stat =solution.stat;
if stat==1||stat==3
    x = solution.cont;
    z = solution.int;
    selected = z(1:nGroups) >= 0.5;
    % A reaction with no binary is one nothing is deciding about, so
    % nothing should discard it.
    reacInd = abs(x(1:n))>=inclusionCutoff | ungoverned;
else
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('MILP solution may not be optimal')
    x=[];reacInd=[];selected=[];
end
end
