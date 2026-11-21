function [Model,LPS] = spectraME2(model,core,tol,weights,nSol,remGene)
% USAGE:
%   [Model,LPS] = spectraME2(model,core,tol,weights,nSol,remGene)
%
% INPUTS:
%   model:   COBRA model structure. The model has to be a consistent model
%            (either topologically or stoichiometrically). Run spectraCC
%            before running spectraME2.
%   core:    Indices of reactions that has to be present in the final
%            model
%
% OPTIONAL INPUTS:
%   tol:     tolerance level (minimum absolute flux that has to be carried
%            by all the reactions in the model) (Default: 1e-4)
%   weights: weights for non-core reactions. (Default: ones)
%   nSol:    Number of alternative solutions required (Default: 1)
%   remGene: Bool value indicating whether to remove the unused genes
%            or not (Default: 0 (doesn't remove the unused genes))
%
% OUTPUTS:
%   Model:     The consistent model with no blocked reactions and has
%              all the core reactions in them (If nSol==1). A cell 
%              consisting of consistent models in them (If nSol >1)
%   LPS:       Number of LPs used to get Model 
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
if ~exist('remGene', 'var') || isempty(remGene)
    remGene=0;  
end

[~,n] = size(model.S);
core = ismember([1:n],core)';
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
rev = true(n,1);
rev(model.lb>=0) = false;
model.rev=rev;

% considering only thermodynamic constraints
model.lb(~rev)=0;
model.lb(rev)=-Inf;
model.ub(:)=Inf;

temp_core = core;
flux = zeros(n,1);
LPS=0;
stalls = 0;
while any(temp_core)
    before = sum(temp_core);
    LPS = LPS+1;
    [flux1,~] = forward(model,temp_core,tol,1);
    if sum(abs(flux))==0
        flux = flux1;
    else
        c1=unifrnd(0.45,0.55,1);
        flux = (c1*flux)+((1-c1)*flux1);
    end
    temp_core(core==1 & abs(flux)>=tol*1e-1) = 0; 
    if ~any(temp_core)
        break
    end
    LPS = LPS+1;
    [flux2,~] = reverse(model,temp_core,tol,1);
    c1=unifrnd(0.45,0.55,1);
    flux = (c1*flux)+((1-c1)*flux2);
    temp_core(core==1 & abs(flux)>=tol*1e-1) = 0; 
    stalls = assertCoreProgress(temp_core,before,stalls,tol);
end

direction = zeros(n,1);
direction(core==1&flux>0) = 1;
direction(core==1&flux<0) = -1;
LPS = LPS+1;
reacInd = minNet(model,direction,weights,tol,1); %LPminimal
Model = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
if remGene
    Model = removeUnusedGenes(Model);
end
if nSol>1
    Model = {Model};
    for j=2:nSol
        temp_core = core;
        flux = zeros(n,1);
        stalls = 0;
        while any(temp_core)
            before = sum(temp_core);
            LPS = LPS+1;
            r = randi(2);
            if r==1
                [flux1,~] = forward(model,temp_core,tol,1);
            else
                [flux1,~] = reverse(model,temp_core,tol,1);
            end
            if sum(abs(flux))==0
                flux = flux1;
            else
                c1=unifrnd(0.45,0.55,1);
                flux = (c1*flux)+((1-c1)*flux1);
            end
            temp_core(core==1 & abs(flux)>=tol*1e-1) = 0; 
            if ~any(temp_core)
                break
            end
            LPS = LPS+1;
            if r==1
                [flux2,~] = reverse(model,temp_core,tol,1);
            else
                [flux2,~] = forward(model,temp_core,tol,1);
            end
            c1=unifrnd(0.45,0.55,1);
            flux = (c1*flux)+((1-c1)*flux2);
            temp_core(core==1 & abs(flux)>=tol*1e-1) = 0; 
            stalls = assertCoreProgress(temp_core,before,stalls,tol);
        end

        direction = zeros(n,1);
        direction(core==1&flux>0) = 1;
        direction(core==1&flux<0) = -1;
        LPS = LPS+1;
        reacInd = minNet(model,direction,weights,tol,1); %LPminimal
        Mod = removeRxns(model, setdiff(model.rxns,model.rxns(reacInd)));
        if remGene
            Mod = removeUnusedGenes(Mod);
        end
        Model = [Model;Mod];
    end
end