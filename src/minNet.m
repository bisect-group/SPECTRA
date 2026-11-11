function [reacInd,x,stat,selected] = minNet(model,direction,weights,tol,steadyState,probType,solveTime,x0,prevSols,indicatorGroups,inclusionCutoff)
% USAGE:
%   [reacInd,x,stat,selected] = minNet(model,direction,weights,tol,steadystate,probType,solveTime,x0,prevSols,indicatorGroups,inclusionCutoff)
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
%                  the probType 'minNetMILP'. Each entry holds reaction
%                  indices, or group indices when indicatorGroups is given.
%     indicatorGroups: Cell array of index vectors, one per group, giving the
%                  reactions that share a single binary. They are then kept
%                  or dropped together and cost the group's weight once
%                  however many of them there are, and weights must hold one
%                  entry per group rather than one per reaction. A free
%                  reaction in no group gets no binary and is always kept.
%                  Sharing a binary is what expresses a requirement holding
%                  in several conditions at once: replicate the network once
%                  per condition, group each reaction's copies together, and
%                  the objective counts reactions while the constraints are
%                  satisfied separately in every condition. That is what
%                  minimalReactome does. Works only for 'minNetMILP'.
%                  (Default: one binary per free reaction)
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

if ~exist('probType', 'var') || isempty(probType)
    probType='minNetLP';     
end

if strcmp(probType,'minNetMILP')    
    if ~exist('solveTime', 'var') || isempty(solveTime)
        solveTime=7200;     
    end
end
if ~exist('indicatorGroups', 'var')
    indicatorGroups = {};
end
if ~exist('inclusionCutoff', 'var') || isempty(inclusionCutoff)
    inclusionCutoff = tol*1e-7;
end
selected = [];

[m,n] = size(model.S);
dir0 = direction==0;
n_ = sum(dir0);

% equalities
% The auxiliary block is one column per binary, which is one per free
% reaction unless indicatorGroups says otherwise, so the MILP branch sizes
% Aeq and the objective for itself.
beq = zeros(m,1);
if steadyState
    csenseeq = repmat('E',m,1); % equality (For consistency based gap filling)
else
    csenseeq = repmat('G',m,1); % greater than (For topology based gap filling)
end

if strcmp(probType,'minNetLP')
    % objective
    f = [zeros(n,1);weights(dir0)];
    Aeq = [model.S, sparse(m,n_)];

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
        reacInd = abs(x(1:n))>=inclusionCutoff;
    end
elseif strcmp(probType,'minNetMILP')
    % Which binary governs which reaction. With no groups every free
    % reaction gets its own, which is the usual formulation; with groups the
    % members of one group share a binary and a free reaction in no group
    % gets none and is always kept.
    [groupOf, nGroups, ungoverned] = resolveIndicatorGroups(indicatorGroups, dir0, n);
    governed = find(groupOf > 0);
    nGov = numel(governed);

    if isempty(indicatorGroups)
        groupWeights = weights(dir0);
    else
        if numel(weights) ~= nGroups
            error('minNet:groupWeightCount', ...
                ['With indicatorGroups, weights must hold one entry per ' ...
                 'group, but got %d weights for %d groups. What a group ' ...
                 'costs to keep is looked up by the group, not by its ' ...
                 'members.'], numel(weights), nGroups);
        end
        groupWeights = weights(:);
    end

    % objective: minimise the weighted count of binaries that are on
    f = [zeros(n,1);groupWeights];

    % equalities
    Aeq = [model.S, sparse(m,nGroups)];

    % inequalities: lb_i z_g <= v_i <= ub_i z_g, so z_g = 0 pins every
    % reaction the group governs to zero at once. Picks out the governed
    % reaction on the left and its own group's binary on the right.
    pick = sparse(1:nGov, governed, 1, nGov, n);
    share = sparse(1:nGov, groupOf(governed), 1, nGov, nGroups);

    Aineq1 = [pick, -spdiags(model.lb(governed),0,nGov,nGov)*share];
    bineq1 = zeros(nGov,1);
    csenseineq1 = repmat('G',nGov,1); % greater than

    Aineq2 = [pick, -spdiags(model.ub(governed),0,nGov,nGov)*share];
    bineq2 = zeros(nGov,1);
    csenseineq2 = repmat('L',nGov,1); % lesser than

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
            Aineq3 = [Aineq3;[zeros(1,n),temp]];
            bineq3 = [bineq3;sum(temp)-1];
            csenseineq3 = [csenseineq3;'L'];
        end
    end

    % bounds
    lb = model.lb;
    lb(direction==1)=max([tol*ones(sum(direction==1),1),lb(direction==1)],[],2);
    lb = [lb;zeros(nGroups,1)];
    ub = model.ub;
    ub(direction==-1)=-tol*ones(sum(direction==-1),1);
    ub = [ub;ones(nGroups,1)];

    % Set up MILP problem
    MILPproblem.A=[Aeq;Aineq1;Aineq2;Aineq3];
    MILPproblem.b=[beq;bineq1;bineq2;bineq3];
    MILPproblem.lb=lb;
    MILPproblem.ub=ub;
    MILPproblem.c=f;
    MILPproblem.osense=1; % minimise
    MILPproblem.vartype = [repmat('C',n,1);repmat('B',nGroups,1)];
    MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2;csenseineq3];
    changeCobraSolverParams('MILP', 'feasTol', 1e-9);
    fprintf(['minNet: MILP feasTol set to 1e-9. This is a global COBRA ' ...
             'setting and stays in effect after this function returns.\n']);
    if exist('x0', 'var')&&~isempty(x0)
        warm = false(nGroups,1);
        on = false(n,1);
        on(1:min(n,numel(x0))) = abs(x0(1:min(n,numel(x0))))>1e-12;
        warm(groupOf(governed(on(governed)))) = true;
        MILPproblem.x0  = [x0;double(warm)];
    end
    solution = solveCobraMILP(MILPproblem,'timeLimit', solveTime);
    stat =solution.stat;
    if stat==1||stat==3
        x = solution.cont;
        z = solution.int;
        selected = z >= 0.5;
        % A reaction with no binary is one nothing is deciding about, so
        % nothing should discard it.
        reacInd = abs(x(1:n))>=inclusionCutoff | ungoverned;
    else
        fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
        fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
        warning('MILP solution may not be optimal')
        x=[];z=[];reacInd=[];selected=[];
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
        reacInd = abs(x(1:n))>=inclusionCutoff;
    end
end
end
