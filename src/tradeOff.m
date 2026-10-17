function [reacInd,x,stat] = tradeOff(model,direction,weights,tol,steadyState,solveTime,x0,prevSols)
% USAGE:
%   [reacInd,x,stat] = tradeOff(model,direction,weights,tol,steadystate,solveTime,x0,prevSols)
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
%                  excluded in the current solution.
%
% OUTPUTS:
%     reacInd: Reaction IDs corresponding to reactions that has to be
%              present in the final extracted model
%     x:       Solution returned after optimization
%     stat:    Status of optimization
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


if ~exist('solveTime', 'var') || isempty(solveTime)
   solveTime=7200;     
end

[m,n] = size(model.S);
dir0 = direction==0; % boolean value indicating if the reaction is non-core (1) or not (0)
n_ = sum(dir0); % number of non core reactions
n_rev = sum(dir0 & model.rev); % number of noncore reversible reactions

% the decision variables are v,z,a,b. v defines the flux values. z define
% the binary integer for the noncore reactions. a and b are for
% reverse and forward direction of the noncore reversible reactions. 

% objective
f = [zeros(n,1);weights(dir0);zeros(n_rev*2,1)];

% stoichiometry or topology constraints
Aeq1 = [model.S, sparse(m,n_+(n_rev*2))];
beq1 = zeros(m,1);

if steadyState
    csenseeq1 = repmat('E',m,1); % equality (For consistency based model extraction)
else
    csenseeq1 = repmat('G',m,1); % greater than (For topology based model extraction)
end


% relation between the flux of non-core irreversible reactions and the
% corresponding binary variables
temp1 = -1*speye(n);
temp2 = spdiag(tol*ones(n_,1));
% getting the ids of irrev reactions in the non-core reactions
temp = ~model.rev(dir0);
temp1 = temp1(~model.rev & dir0,:);
temp2 = temp2(temp,:);

Aineq1 = [temp1,temp2,sparse(sum(temp),n_rev*2)];
bineq1 = zeros(sum(temp),1);
csenseineq1 = repmat('L',sum(temp),1); % lesser than

temp2 = -1*spdiag(model.ub(dir0));
temp2 = temp2(temp,:);

Aineq2 = [-1*temp1,temp2,sparse(sum(temp),n_rev*2)];
bineq2 = zeros(sum(temp),1);
csenseineq2 = repmat('L',sum(temp),1); % lesser than

% constraints to ignore the previously obtained solutions
Aineq3=[];bineq3=[];csenseineq3=[];
if exist('prevSols', 'var') && ~isempty(prevSols)
    for id=1:numel(prevSols)
        temp = find(dir0);
        temp = ismember(temp,prevSols{id});
        temp =temp';
        Aineq3 = [Aineq3;[zeros(1,n),temp,zeros(1,n_rev*2)]];
        bineq3 = [bineq3;sum(temp)-1];
        csenseineq3 = [csenseineq3;'L'];
    end
end

% constraints for reversible non-core reactions to carry flux only in
% one direction. Either in forward direction of reverse direction. 
% a + b = z
temp1 = speye(n_);
temp = find(model.rev(dir0));
temp1 = temp1(temp,:);

Aeq2 = [sparse(numel(temp),n), -1*temp1, speye(n_rev), speye(n_rev)];
beq2 = zeros(numel(temp),1);
csenseeq2 = repmat('E',numel(temp),1); 

% constraint linking the binary variable, a and the flux through the
% reversible non-core reactions

temp1 = speye(n);
temp1 = temp1(dir0 & model.rev,:); % ids of non-core reversible reactions
temp = sum(dir0 & model.rev);
Aineq4 = [temp1, sparse(temp,n_), -1*spdiag(tol*ones(temp,1)), -1*spdiag(model.lb(dir0 & model.rev))];
bineq4 = zeros(temp,1);
csenseineq4 = repmat('G',temp,1);

% constraint linking the binary variable, b and the flux through the
% reversible non-core reactions

temp1 = speye(n);
temp1 = temp1(dir0 & model.rev,:); % ids of non-core reversible reactions
temp = sum(dir0 & model.rev);
Aineq5 = [temp1, sparse(temp,n_), -1*spdiag(model.ub(dir0 & model.rev)), spdiag(tol*ones(temp,1))];
bineq5 = zeros(temp,1);
csenseineq5 = repmat('L',temp,1);

% bounds
lb = model.lb;
lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
lb = [lb;zeros(n_+(2*n_rev),1)];
ub = model.ub;
ub(direction==-1)=-tol*ones(sum(direction==-1),1);
ub = [ub;ones(n_+(2*n_rev),1)];

% Set up MILP problem
MILPproblem.A=[Aeq1;Aeq2;Aineq1;Aineq2;Aineq3;Aineq4;Aineq5];
MILPproblem.b=[beq1;beq2;bineq1;bineq2;bineq3;bineq4;bineq5];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=-1; % maximize
MILPproblem.vartype = [repmat('C',n,1);repmat('B',n_+(2*n_rev),1)];
MILPproblem.csense = [csenseeq1; csenseeq2; csenseineq1; csenseineq2; csenseineq3; csenseineq4; csenseineq5];
changeCobraSolverParams('MILP', 'feasTol', 1e-9);
fprintf(['tradeOff: MILP feasTol set to 1e-9. This is a global COBRA ' ...
         'setting and stays in effect after this function returns.\n']);
if exist('x0', 'var')&&~isempty(x0)

    MILPproblem.x0  = [x0;abs(x0(dir0))>1e-12;x0(dir0&model.rev)>1e-12;x0(dir0&model.rev)<-1e-12]; 
end
solution = solveCobraMILP(MILPproblem,'timeLimit', solveTime);
stat =solution.stat;
if stat==1||stat==3
    x = solution.cont;
    z = solution.int;
    reacInd = abs(x(1:n))>=tol*1e-7;
else
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('MILP solution may not be optimal')
    x=[];z=[];reacInd=[];
end
end