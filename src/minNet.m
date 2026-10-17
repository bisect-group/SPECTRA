function [reacInd,x,stat] = minNet(model,direction,weights,tol,steadyState,probType,solveTime,x0,prevSols)
% USAGE:
%   [reacInd,x,stat] = minNet(model,direction,weights,tol,steadystate,probType,solveTime,x0,prevSols)
%
% INPUTS:
%     model:       COBRA model structure.
%     direction:   A vector of size equal to number of reactions in the
%                  model. The unique elements in this vector has to be -1,0
%                  and 1 defining the flux directionality info.
%     weights:     weights for non-core reactions. More the weights, lesser
%                  the chance to get included in the final model 
%     tol:         tolerance level (minimum absolute flux that has to be carried
%                  by a reaction for it to be defined as consistent)
%     steadyState: Boolean value indicating whether to assume steady state
%                  condition (S.v = 0) or accumulation condition (S.v >= 0)
%
% OPTIONAL INPUTS:
%     probType:    'minNetLP' (Default): Linear Programming, 'minNetMILP': Mixed Integer Linear
%                  Programming, 'minNetDC': Difference of convex functions
%     solveTime:   If probType is 'minNetMILP', solveTime refers to the upper limit
%                  of solving time
%     x0:          solution vector to initialize with for MILP problem. 
%                  This can be obtained from the LPforward and LPreverse
%     prevSols:    A cell of previously obtained solutions that needs to be
%                  exclude in the current solution. Note: This works only for
%                  the probType 'minNetMILP'
%
% OUTPUTS:
%     reacInd: Reaction IDs corresponding to reactions that has to be
%              present in the final extracted model
%     x:       Solution returned after optimization
%     stat:    Status of optimization
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if ~exist('probType', 'var') || isempty(probType)
    probType='minNetLP';     
end

if strcmp(probType,'minNetMILP')    
    if ~exist('solveTime', 'var') || isempty(solveTime)
        solveTime=7200;     
    end
end

[m,n] = size(model.S);
dir0 = direction==0;
n_ = sum(dir0);

% objective
f = [zeros(n,1);weights(dir0)];

% equalities
Aeq = [model.S, sparse(m,n_)];
beq = zeros(m,1);
if steadyState
    csenseeq = repmat('E',m,1); % equality (For consistency based gap filling)
else
    csenseeq = repmat('G',m,1); % greater than (For topology based gap filling)
end

if strcmp(probType,'minNetLP')
    % inequalities
    temp1 = speye(n);
    temp2 = speye(n_);
    Aineq1 = [temp1(dir0,:),temp2];
    bineq1 = zeros(n_,1);
    csenseineq1 = repmat('G',n_,1); % greater than

    Aineq2 = [temp1(dir0,:),-1*temp2];
    bineq2 = zeros(n_,1);
    csenseineq2 = repmat('L',n_,1); % lesser than

    % bounds
    lb = model.lb;
    lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
    lb = [lb;zeros(n_,1)];
    ub = model.ub;
    ub(direction==-1)=-tol*ones(sum(direction==-1),1);
    ub = [ub;Inf(n_,1)];

    % Set up LP problem
    LPproblem.A=[Aeq;Aineq1;Aineq2];
    LPproblem.b=[beq;bineq1;bineq2];
    LPproblem.lb=lb;
    LPproblem.ub=ub;
    LPproblem.c=f;
    LPproblem.osense=1;%minimise
    LPproblem.csense = [csenseeq; csenseineq1; csenseineq2];
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
elseif strcmp(probType,'minNetMILP')
   % inequalities
    temp1 = speye(n);
    temp2 = -1*spdiag(model.lb(dir0));
    Aineq1 = [temp1(dir0,:),temp2];
    bineq1 = zeros(n_,1);
    csenseineq1 = repmat('G',n_,1); % greater than

    temp2 = -1*spdiag(model.ub(dir0));
    Aineq2 = [temp1(dir0,:),temp2];
    bineq2 = zeros(n_,1);
    csenseineq2 = repmat('L',n_,1); % lesser than
    
    
    Aineq3=[];bineq3=[];csenseineq3=[];
    if exist('prevSols', 'var') && ~isempty(prevSols)
        for id=1:numel(prevSols)
            temp = find(dir0);
            temp = ismember(temp,prevSols{id});
            temp =temp';
            Aineq3 = [Aineq3;[zeros(1,n),temp]];
            bineq3 = [bineq3;sum(temp)-1];
            csenseineq3 = [csenseineq3;'L'];
        end
    end

    % bounds
    lb = model.lb;
    lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
    lb = [lb;zeros(n_,1)];
    ub = model.ub;
    ub(direction==-1)=-tol*ones(sum(direction==-1),1);
    ub = [ub;ones(n_,1)];

    % Set up MILP problem
    MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3];
    MILPproblem.b=[beq;bineq1;bineq2;bineq3];
    MILPproblem.lb=lb;
    MILPproblem.ub=ub;
    MILPproblem.c=f;
    MILPproblem.osense=1; % minimise
    MILPproblem.vartype = [repmat('C',n,1);repmat('B',n_,1)];
    MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2;csenseineq3];
    changeCobraSolverParams('MILP', 'feasTol', 1e-9);
    fprintf(['minNet: MILP feasTol set to 1e-9. This is a global COBRA ' ...
             'setting and stays in effect after this function returns.\n']);
    if exist('x0', 'var')&&~isempty(x0)
        MILPproblem.x0  = [x0;abs(x0(dir0))>1e-12]; 
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
    % reacInd = abs(x(1:n))~=0;    
elseif strcmp(probType,'minNetDC')
    % equalities
    Aeq = model.S;
    beq = zeros(m,1);
    if steadyState
        csenseeq = repmat('E',m,1); % equality (For consistency based gap filling)
    else
        csenseeq = repmat('G',m,1); % greater than (For topology based gap filling)
    end
    % bounds
    lb = model.lb;
    lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
    
    ub = model.ub;
    ub(direction==-1)=-tol*ones(sum(direction==-1),1);
   

    % Set up MILP problem
    DCproblem.A=Aeq;
    DCproblem.b=beq;
    DCproblem.csense = csenseeq;
    DCproblem.lb=lb;
    DCproblem.ub=ub;
    DCproblem.p = direction==0;
    DCproblem.q = false(n,1);
    DCproblem.r = ~DCproblem.p;
    DCproblem.c = zeros(n,1);
    DCproblem.osense=1;
    DCproblem.lambda0=1;
    DCproblem.k = weights;
    solution = optimizeCardinality(DCproblem);
    stat=solution.stat;
    if stat~=1
        fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
        fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
        warning('DC solution may not be optimal')
        x =[];reacInd=[];
    else
        x=zeros(n,1);
        x(direction==0) = solution.x;
        x(direction~=0) = solution.z;
        reacInd = abs(x(1:n))>=tol*1e-7;
    end
end
end