clear
initCobraToolbox(0)
changeCobraSolver('gurobi','all')
changeCobraSolverParams('LP', 'feasTol', 1e-8);
load('./DepMapData')
context = geneExpression.context;
load('./consRecon3DGeneSymbol.mat')
bio_id = find(ismember(model.rxns,'biomass_reaction'));
ATP_id = find(ismember(model.rxns,'DM_atp_c_'));
coreRxn = [bio_id;ATP_id];
geneExpression.value=geneExpression.value';
[RxnImp,Contexts] = GiniReactionImportance(geneExpression,model,90,10,coreRxn);
RxnImp(isnan(RxnImp))=0;
model = rmfield(model,{'metCharges','metFormulas','metSmiles','metInChIString',...
    'metKEGGID','metPubChemID','rxnNotes','rxnECNumbers','rxnReferences',...
    'rxnKEGGID','metCHEBIID','metPdMap','rxnCOG','rxnKeggOrthology'});

% creating directories for saving models
if ~exist('./fastcore','dir')
    mkdir('./fastcore')
end
if ~exist('./swiftcore','dir')
    mkdir('./swiftcore')
end
if ~exist('./spectra','dir')
    mkdir('./spectra')
end

FC_times = []; SP_times=[]; SW_times=[];
for i=1:size(RxnImp,2)
    curr_core = RxnImp(:,i);
    % weights are for those non-core reactions more the weights, less the
    % probabilty of getting added
    weights = zeros(numel(curr_core),1);
    weights(curr_core<1)=1-curr_core(curr_core<1);
    % core reactions are those with RxnImp more than 1
    core = find(curr_core>=1);
    
    % FASTCORE (weighted) model using LocalGini non-core penalty weights
    tic
    A_fcw = fastCoreWeighted(core, model, weights, 1e-4);
    m = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), A_fcw)));
    t = toc;
    FC_times(i) = t;
    save(['./fastcore/',context{i}],"m")

    % SWIFTCORE model
    tic
    m = swiftcore(model,core,weights,1e-10,1);
    t=toc;
    SW_times(i)=t;
    save(['./swiftcore/',context{i}],"m")
    
    % SPECTRA model
    tic
    m =  spectraME(model,core,1e-4,[],weights);
    t=toc;
    SP_times(i)=t;
    save(['./spectra/',context{i}],"m")
end
clear m
save('Results_BuildDepMapModels')
