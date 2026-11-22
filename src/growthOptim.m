function [reacInd,x,stat] = growthOptim(model,direction,weights,tol,steadyState,useCoupling)
% USAGE:
%   [reacInd,x,stat] = growthOptim(model,direction,weights,tol,steadyState)
%
% INPUTS:
%     model:       COBRA model structure.
%     direction:   A vector of size equal to number of reactions in the
%                  model. The unique elements in this vector has to be -1,0
%                  and 1 defining the flux directionality info.
%     weights:     weights for all the reactions. 
%     tol:         tolerance level (minimum absolute flux that has to be carried
%                  by a reaction for it to be defined as consistent)
%     steadyState: Boolean value indicating whether to assume steady state
%                  condition (S.v = 0) or accumulation condition (S.v >= 0)
%
% OUTPUTS:
%     reacInd: Reaction IDs corresponding to reactions that has to be
%              present in the final extracted model
%     x:       Solution returned after optimization
%     stat:    Status of optimization
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras



[m,n] = size(model.S);

temp = weights;
weights = -1*weights;
weights(find(model.c)) = 0;
c = zeros(n,1);
c(find(model.c)) = temp(find(model.c));

% objective
f = [c;weights];

% equalities
Aeq = [model.S, sparse(m,n)];
beq = zeros(m,1);
if steadyState
    csenseeq = repmat('E',m,1); % equality (For consistency based gap filling)
else
    csenseeq = repmat('G',m,1); % greater than (For topology based gap filling)
end

% inequalities
temp1 = speye(n);
temp2 = speye(n);
Aineq1 = [temp1,temp2];
bineq1 = zeros(n,1);
csenseineq1 = repmat('G',n,1); % greater than

Aineq2 = [temp1,-1*temp2];
bineq2 = zeros(n,1);
csenseineq2 = repmat('L',n,1); % lesser than

% bounds
lb = model.lb;
lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
lb = [lb;zeros(n,1)];
ub = model.ub;
ub(direction==-1)=-tol*ones(sum(direction==-1),1);
ub = [ub;Inf(n,1)];

% coupling constraints, where the model has any
if ~exist('useCoupling','var'); useCoupling = false; end
[Acpl,bcpl,csensecpl] = couplingRows(model,n,useCoupling);

% Set up LP problem
LPproblem.A=[Aeq;Aineq1;Aineq2;Acpl];
LPproblem.b=[beq;bineq1;bineq2;bcpl];
LPproblem.lb=lb;
LPproblem.ub=ub;
LPproblem.c=f;
LPproblem.osense=-1;%maximize
LPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csensecpl];
solution = solveCobraLP(LPproblem);
stat =solution.stat;
if stat~=1
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('LP solution may not be optimal')
    x =[];reacInd=[];
else
    x=solution.full;
    reacInd = abs(x(1:n))>=tol*1e-7;
end
end