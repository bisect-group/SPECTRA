function [Model,LPS] = spectraME(model,coreRxns,tol,consType,weights,nSol,altSolMethod,probType,solveTime,remGene,prevSols,seed)
% USAGE:
%   [Model,LPS] = spectraME(model,coreRxns,tol,consType,weights,nSol,altSolMethod,probType,solveTime,remGene,prevSols,seed)
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
%
% OUTPUTS:
%   Model:     The consistent model with no blocked reactions and has
%              all the core reactions in them (If nSol==1). A cell 
%              consisting of consistent models in them (If nSol > 1)
%   LPS:       Number of LPs used to get ConsModel, if probType is 'minNetLP'
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

[~,n] = size(model.S);
coreRxns = ismember(1:n,coreRxns)';
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
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

while any(temp_core)
    LPS = LPS+1;
    [flux1,~] = forward(model,temp_core,tol,steadystate);
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
    c1=unifrnd(0.45,0.55,1);
    flux = (c1*flux)+((1-c1)*flux2);
    temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
end


direction = zeros(n,1);
direction(coreRxns==1&flux>0) = 1;
direction(coreRxns==1&flux<0) = -1;

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
    [reacInd,x] = minNet(model,direction,weights,tol,steadystate,probType,solveTime,flux,prevSols);
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
    [reacInd,x] = tradeOff(model,direction,weights,tol,steadystate,solveTime,flux,prevSols);
end

flux = x(1:numel(model.rxns));
Model = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
if remGene
    Model = removeUnusedGenes(Model);
end

if nSol>1
    Model = {Model};
    if strcmp(altSolMethod,'coreDirection')
        LPS = {LPS};    
        for j=2:nSol
            LPS_ = 0;
            temp_core = coreRxns;
            flux = zeros(n,1);
            while any(temp_core)
                LPS_ = LPS_+1;
                r = randi(2);
                if r==1
                    [flux1,~] = forward(model,temp_core,tol,steadystate);
                else
                    [flux1,~] = reverse(model,temp_core,tol,steadystate);
                end
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
                c1=unifrnd(0.45,0.55,1);
                flux = (c1*flux)+((1-c1)*flux2);
                temp_core(coreRxns==1 & abs(flux)>=tol*1e-1) = 0; 
            end
        
            direction = zeros(n,1);
            direction(coreRxns==1&flux>0) = 1;
            direction(coreRxns==1&flux<0) = -1;
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
                [reacInd,x] = minNet(model,direction,weights,tol,steadystate,probType,solveTime,flux);
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
                [reacInd,x] = tradeOff(model,direction,weights,tol,steadystate,solveTime,flux);
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
            prevSols = [prevSols;find(reacInd)];
            if strcmp(probType,'minNetMILP')
                [reacInd,~,stat] = minNet(model,direction,weights,tol,steadystate,probType,solveTime,[],prevSols);
            elseif strcmp(probType,'tradeOff')
                [reacInd,~,stat] = tradeOff(model,direction,weights,tol,steadystate,solveTime,[],prevSols);
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
end