function [Model,LPS,selection] = spectraME(model,coreRxns,tol,consType,weights,nSol,altSolMethod,probType,solveTime,remGene,prevSols,seed,indicatorGroups,inclusionCutoff)
% USAGE:
%   [Model,LPS,selection] = spectraME(model,coreRxns,tol,consType,weights,nSol,altSolMethod,probType,solveTime,remGene,prevSols,seed,indicatorGroups,inclusionCutoff)
%
% INPUTS:
%   model:    COBRA model structure. The model has to be a consistent model
%             (either topologically or stoichiometrically). Run spectraCC
%             before running spectraME.
%   coreRxns: Indices of reactions that has to be present in the final
%             model
%
% OPTIONAL INPUTS:
%   tol:          Tolerance level (minimum absolute flux that has to be carried
%                 by all the reactions in the Model) (Default: 1e-4)
%   consType:     Type of constraints to apply. Either 'topology' or
%                 'stoichiometry' constraints. (Default:'stoichiometry')
%   weights:      Weights for the reactions in the model. (Default: ones)
%   nSol:         Number of alternative solutions required (Default: 1)
%   altSolMethod: Method to find the alternate solutions.
%                 accepted values: 'coreDirection', 'pathwayExclusion'.
%                 Note: 'pathwayExclusion' works only for MILP problems
%                 (Default: 'coreDirection')
%   probType:     Which method to use to extract the final model.
%                 Accepted values: 'minNetLP','minNetMILP','minNetDC','growthOptim','tradeOff'.
%                 (Default: 'minNetLP')
%   solveTime:    maximum runtime for solving MILP problem (Default: 7200s)
%   remGene:      Bool value indicating whether to remove the unused genes
%                 or not (Default: 0 (doesn't remove the unused genes))
%   prevSols:     A cell of previously obtained solutions that should not
%                 be a part of any new solutions (works only for MILP probType)
%   seed:         Seed for the random coefficients used by the direction
%                 phase, making the result reproducible across runs. Leave
%                 it out to draw them unseeded, as every earlier version of
%                 this function did. (Default: none)
%   indicatorGroups: Which reactions share a binary in the MILP, as either
%                 a cell array of reaction index vectors (one per group) or
%                 a cell array of reaction IDs (one binary each, every other
%                 free reaction kept outright). Naming a handful of
%                 reactions is what makes a selection problem tractable on a
%                 large model: one binary per organism of a community, say,
%                 rather than one per reaction. Grouping several together
%                 keeps or drops them as one, which is how a requirement
%                 holding in several conditions at once is expressed. Works
%                 only for 'minNetMILP' and 'tradeOff'. (Default: one binary
%                 per free reaction)
%   inclusionCutoff: Absolute flux at which a reaction counts as included.
%                 Handed to the formulation. (Default: its own, tol*1e-7)
%
% OUTPUTS:
%   Model:     The consistent model with no blocked reactions and has
%              all the core reactions in them (If nSol==1). A cell 
%              consisting of consistent models in them (If nSol > 1)
%   LPS:       Number of LPs used to get ConsModel, if probType is 'minNetLP'
%   selection: What the MILP's binaries decided, as a structure with fields
%              groups (the resolved groups, as index vectors), selected (a
%              logical, one entry per group), selectedRxns (the IDs of the
%              reactions in the groups that are on) and ungoverned (a
%              logical over the reactions, true where no binary governed the
%              reaction and it was therefore kept). Empty for a formulation
%              with no binaries. Reading the membership off the model is not
%              equivalent: a reaction whose binary is off can still carry a
%              trace of flux, and one that is on need not carry any.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


if ~exist('tol', 'var') || isempty(tol)
    tol=1e-4;     
end
if ~exist('weights', 'var') || isempty(weights)
    weights = ones(numel(model.rxns),1);  
end
if ~exist('nSol', 'var') || isempty(nSol)
    nSol=1;  
end
if ~exist('probType', 'var') || isempty(probType)
    probType='minNetLP';  
end
if ~exist('remGene', 'var') || isempty(remGene)
    remGene=0;  
end
if ~exist('prevSols', 'var') || isempty(prevSols)
    prevSols={};
end
if ~exist('consType', 'var') || isempty(consType)
    consType='stoichiometry';
end
if ~exist('altSolMethod', 'var') || isempty(altSolMethod)
    altSolMethod='coreDirection';
end
% Seeding is opt-in. Without a seed the coefficients are drawn exactly as
% before, so published results stay reproducible.
if exist('seed', 'var') && ~isempty(seed)
    rng(seed);
end
if ~exist('indicatorGroups', 'var')
    indicatorGroups = {};
end
if ~exist('inclusionCutoff', 'var')
    inclusionCutoff = [];
end

[~,n] = size(model.S);
coreRxns = ismember(1:n,coreRxns)';
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
% A coupling constraint is written over the same fluxes, so flipping a
% reaction's column in S without flipping it in C would silently change
% what the constraint says about that reaction.
if isfield(model,'C') && ~isempty(model.C)
    model.C(:,IrR) = -model.C(:,IrR);
end
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
rev = ones(n,1);
rev(model.lb>=0) = 0;
model.rev=rev;

temp_core = coreRxns;
flux = zeros(n,1);
LPS=0;

if strcmp(consType,'stoichiometry')
    steadystate = 1;
elseif strcmp(consType,'topology')
    steadystate = 0;
else
    error('spectraME:badConsType', ...
        ['consType must be ''stoichiometry'' or ''topology'', ' ...
         'not ''%s''.'], consType);
end

% Indicator groups belong to the MILP formulations: they are groupings of
% the binaries, and a formulation without binaries has nothing to group.
[indicatorGroups, milpWeights] = resolveGroupArgument(indicatorGroups, ...
    weights, model.rxns, probType);

stalls = 0;
while any(temp_core)
    before = sum(temp_core);
    LPS = LPS+1;
    [flux1,~] = forward(model,temp_core,tol,steadystate);
    assertCoreFeasible(flux1,sum(coreRxns),tol);
    if sum(abs(flux))==0
        flux = flux1;
    else
        c1=unifrnd(0.45,0.55,1);
        flux = (c1*flux)+((1-c1)*flux1);
    end
    temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
    if ~any(temp_core)
        break
    end
    LPS = LPS+1;
    [flux2,~] = reverse(model,temp_core,tol,steadystate);
    assertCoreFeasible(flux2,sum(coreRxns),tol);
    c1=unifrnd(0.45,0.55,1);
    flux = (c1*flux)+((1-c1)*flux2);
    temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
    stalls = assertCoreProgress(temp_core,before,stalls,tol);
end


direction = zeros(n,1);
direction(coreRxns==1&flux>0) = 1;
direction(coreRxns==1&flux<0) = -1;

selected = [];
if strcmp(probType,'minNetMILP')
    if ~exist('solveTime', 'var') || isempty(solveTime)
        solveTime=7200;     
    end
    % if the core reactions is empty. This will be further used as an initial
    % solution for the minNetMILP
    if sum(coreRxns)==0
        sol = optimizeCbModel(model);
        flux = sol.x;
    end
    [reacInd,x,~,selected] = minNet(model,direction,milpWeights,tol,steadystate,probType,solveTime,flux,prevSols,indicatorGroups,inclusionCutoff);
elseif strcmp(probType,'minNetLP')
    LPS = LPS+1;
    [reacInd,x] = minNet(model,direction,weights,tol,steadystate,probType);
elseif strcmp(probType,'minNetDC')
    [reacInd,x] = minNet(model,direction,weights,tol,steadystate,probType);
elseif strcmp(probType,'growthOptim')
    [reacInd,x] = growthOptim(model,direction,weights,tol,steadystate);
elseif strcmp(probType,'tradeOff')
    if ~exist('solveTime', 'var') || isempty(solveTime)
        solveTime=7200;     
    end
    % if the core reactions is empty. This will be further used as an initial
    % solution for the tradeoff
    if sum(coreRxns)==0
        sol = optimizeCbModel(model);
        flux = sol.x;
    end
    [reacInd,x,~,selected] = tradeOff(model,direction,milpWeights,tol,steadystate,solveTime,flux,prevSols,indicatorGroups,inclusionCutoff);
end

flux = x(1:numel(model.rxns));
Model = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
if remGene
    Model = removeUnusedGenes(Model);
end
selection = reportSelection(indicatorGroups,selected,model.rxns,direction);

if nSol>1
    Model = {Model};
    if strcmp(altSolMethod,'coreDirection')
        LPS = {LPS};    
        for j=2:nSol
            LPS_ = 0;
            temp_core = coreRxns;
            flux = zeros(n,1);
            stalls = 0;
            while any(temp_core)
                before = sum(temp_core);
                LPS_ = LPS_+1;
                r = randi(2);
                if r==1
                    [flux1,~] = forward(model,temp_core,tol,steadystate);
                else
                    [flux1,~] = reverse(model,temp_core,tol,steadystate);
                end
                assertCoreFeasible(flux1,sum(coreRxns),tol);
                if sum(abs(flux))==0
                    flux = flux1;
                else
                    c1=unifrnd(0.45,0.55,1);
                    flux = (c1*flux)+((1-c1)*flux1);
                end
                temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
                if ~any(temp_core)
                    break
                end
                LPS_ = LPS_+1;
                if r==1
                    [flux2,~] = reverse(model,temp_core,tol,steadystate);
                else
                    [flux2,~] = forward(model,temp_core,tol,steadystate);
                end
                assertCoreFeasible(flux2,sum(coreRxns),tol);
                c1=unifrnd(0.45,0.55,1);
                flux = (c1*flux)+((1-c1)*flux2);
                temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
                stalls = assertCoreProgress(temp_core,before,stalls,tol);
            end
        
            direction = zeros(n,1);
            direction(coreRxns==1&flux>0) = 1;
            direction(coreRxns==1&flux<0) = -1;
            selected = [];
            if strcmp(probType,'minNetMILP')
                if ~exist('solveTime', 'var') || isempty(solveTime)
                    solveTime=7200;
                end
                % if the core reactions is empty. This will be further used as an initial
                % solution for the minNetMILP
                if sum(coreRxns)==0
                    sol = optimizeCbModel(model);
                    flux = sol.x;
                end
                [reacInd,x,~,selected] = minNet(model,direction,milpWeights,tol,steadystate,probType,solveTime,flux,{},indicatorGroups,inclusionCutoff);
            elseif strcmp(probType,'minNetLP')
                LPS_ = LPS_+1;
                [reacInd,x] = minNet(model,direction,weights,tol,steadystate,probType);
            elseif strcmp(probType,'minNetDC')
                [reacInd,x] = minNet(model,direction,weights,tol,steadystate,probType);
            elseif strcmp(probType,'growthOptim')
                [reacInd,x] = growthOptim(model,direction,weights,tol,steadystate);
            elseif strcmp(probType,'tradeOff')
                if ~exist('solveTime', 'var') || isempty(solveTime)
                    solveTime=7200;
                end
                % if the core reactions is empty. This will be further used as an initial
                % solution for the tradeoff
                if sum(coreRxns)==0
                    sol = optimizeCbModel(model);
                    flux = sol.x;
                end
                [reacInd,x,~,selected] = tradeOff(model,direction,milpWeights,tol,steadystate,solveTime,flux,{},indicatorGroups,inclusionCutoff);
            end
    
            flux = x(1:numel(model.rxns));
            Mod = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
            if remGene
                Mod = removeUnusedGenes(Mod);
            end 
            Model = [Model;Mod];
            LPS = [LPS;LPS_];
        end
    elseif strcmp(altSolMethod,'pathwayExclusion') % the problem has to be an milp problem
        for j=2:nSol
            if isempty(indicatorGroups)
                prevSols = [prevSols;find(reacInd)];
            else
                prevSols = [prevSols;{find(selected)}];
            end
            if strcmp(probType,'minNetMILP')
                [reacInd,~,stat,selected] = minNet(model,direction,milpWeights,tol,steadystate,probType,solveTime,[],prevSols,indicatorGroups,inclusionCutoff);
            elseif strcmp(probType,'tradeOff')
                [reacInd,~,stat,selected] = tradeOff(model,direction,milpWeights,tol,steadystate,solveTime,[],prevSols,indicatorGroups,inclusionCutoff);
            end
            if stat~=1
                break
            end
            Mod = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
            if remGene
                Mod = removeUnusedGenes(Mod);
            end
            Model = [Model;Mod];
        end
    end
    if ~isempty(selected)
        selection = reportSelection(indicatorGroups,selected,model.rxns,direction);
    end
end


function [groups,milpWeights] = resolveGroupArgument(given,weights,rxns,probType)
% Turns the indicatorGroups argument into index vectors and group weights.
%
% Reaction IDs are accepted as well as indices, because naming the handful
% of reactions that get a binary is how the argument is usually used, and
% looking the indices up by hand is a step that only invites mistakes.
groups = {};
milpWeights = weights;
if isempty(given)
    return
end
if ~ismember(probType,{'minNetMILP','tradeOff'})
    error('spectraME:groupsWithoutBinaries', ...
        ['indicatorGroups groups the MILP''s binaries, and ''%s'' has ' ...
         'none. Use ''minNetMILP'' or ''tradeOff''.'], probType);
end
if ~iscell(given)
    error('spectraME:badIndicatorGroups', ...
        ['indicatorGroups must be a cell array, of reaction index ' ...
         'vectors or of reaction IDs.']);
end
named = all(cellfun(@(g) ischar(g) || (isstring(g) && isscalar(g)), given(:)'));
if named
    ids = cellstr(given(:));
    [found,at] = ismember(ids,rxns);
    if ~all(found)
        missing = ids(~found);
        error('spectraME:unknownIndicatorReaction', ...
            'No reaction %s in the model, so it cannot be given a binary.', ...
            missing{1});
    end
    groups = num2cell(at(:));
else
    groups = cellfun(@(g) g(:), given(:), 'UniformOutput', false);
end

% A weight per reaction is what every other caller passes, so take it and
% add up each group's. One binary per reaction leaves every weight as it
% was; the copies of one reaction across replicated conditions share a
% weight, so the sum rescales every group alike and leaves the argument of
% the minimum where it was.
if numel(weights) == numel(rxns) && numel(groups) ~= numel(weights)
    milpWeights = cellfun(@(g) sum(weights(g)), groups);
    milpWeights = milpWeights(:);
end


function selection = reportSelection(indicatorGroups,selected,rxns,direction)
% Packs what the binaries decided into a structure, or [] where there were none.
if isempty(selected)
    selection = [];
    return
end
[groupOf,nGroups,ungoverned] = resolveIndicatorGroups(indicatorGroups, ...
    direction == 0, numel(rxns));
groups = cell(nGroups,1);
for g = 1:nGroups
    groups{g} = find(groupOf == g);
end
on = false(numel(rxns),1);
for g = reshape(find(selected(:)),1,[])
    on(groups{g}) = true;
end
selection = struct('groups',{groups},'selected',{logical(selected(:))}, ...
    'selectedRxns',{rxns(on)},'ungoverned',{ungoverned});


function assertCoreFeasible(flux,nCore,tol)
% Stops a direction LP that came back infeasible from becoming a hang.
%
% forward and reverse return an all-NaN flux vector when their LP has no
% solution. NaN fails every comparison, so the caller's loop clears no core
% reaction, `while any(temp_core)` never exits, and the run spins emitting a
% warning per iteration instead of stopping. Checking here turns that into a
% message saying what to do about it.
if any(isnan(flux))
    error('spectraME:infeasibleCore', ...
        ['No flux distribution lets all %d core reactions carry at least ' ...
         'tol=%g at once. Either the core set is mutually inconsistent, or ' ...
         'some of its reactions are blocked in this model. spectraCCME ' ...
         'handles the latter by dropping them and returning them as ' ...
         'BlockedCoreRxns.'], nCore, tol);
end
