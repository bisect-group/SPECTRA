%% s02_build_models.m
% Reviewer-2 response, Phase 1b.
%
% Rebuilds the 20 SPECTRA (minNetLP) context-specific models for the 10
% cancer/non-cancer pairs. Follows BuildCancerModels/DepMapModels_Recon/
% BuildDepMapModels.m exactly, with one deliberate difference: only the 20
% cell lines of interest are extracted rather than all 1,479.
%
% NOTE ON THE GINI STEP: GiniReactionImportance computes, per gene, a
% percentile threshold taken ACROSS ALL CONTEXTS in the matrix it is handed.
% Subsetting to 20 cell lines before that call would therefore change every
% threshold and produce different models. The Gini step is run on the full
% 1,479-sample matrix and the 20 columns are taken afterwards.

clear
root = '/home/pavan/scratch/SPECTRA_review_2';
addpath('/home/pavan/scratch/cobratoolbox');
addpath(genpath('/home/pavan/scratch/spectra/src'));
addpath(fullfile(root,'scripts'));

initCobraToolbox(0);
changeCobraSolver('gurobi','all');
changeCobraSolverParams('LP','feasTol',1e-8);

run(fullfile(root,'scripts','config_pairs.m'));

modelDir = fullfile(root,'models','spectra');
if ~exist(modelDir,'dir'); mkdir(modelDir); end

impFile = fullfile(root,'data','RxnImp_20.mat');

if exist(impFile,'file')
    fprintf('Loading cached reaction importances from %s\n', impFile);
    load(impFile,'RxnImp20','sel','model');
else
    fprintf('Loading expression data ...\n');
    load(fullfile(root,'data','DepMapData.mat'),'geneExpression');
    load(fullfile(root,'data','consRecon3DGeneSymbol.mat'),'model');

    bio_id  = find(ismember(model.rxns,'biomass_reaction'));
    ATP_id  = find(ismember(model.rxns,'DM_atp_c_'));
    assert(~isempty(bio_id) && ~isempty(ATP_id), ...
        'biomass_reaction / DM_atp_c_ not found in the model.');
    coreRxn = [bio_id; ATP_id];

    geneExpression.value = geneExpression.value';   % -> genes x samples

    fprintf('Computing Gini reaction importances over all %d contexts ...\n', ...
        numel(geneExpression.context));
    t0 = tic;
    [RxnImp, Contexts] = GiniReactionImportance(geneExpression, model, 90, 10, coreRxn);
    RxnImp(isnan(RxnImp)) = 0;
    fprintf('  done in %.1f s\n', toc(t0));

    [tf, sel] = ismember(allModels, Contexts);
    assert(all(tf), 'Some of the 20 cell lines are missing from the expression matrix.');
    RxnImp20 = RxnImp(:, sel);

    model = rmfield(model,{'metCharges','metFormulas','metSmiles','metInChIString', ...
        'metKEGGID','metPubChemID','rxnNotes','rxnECNumbers','rxnReferences', ...
        'rxnKEGGID','metCHEBIID','metPdMap','rxnCOG','rxnKeggOrthology'});

    save(impFile,'RxnImp20','sel','model','allModels','-v7.3');
    clear RxnImp geneExpression
end

SP_times = nan(numel(allModels),1);
for i = 1:numel(allModels)
    outFile = fullfile(modelDir, [allModels{i} '.mat']);
    if exist(outFile,'file')
        fprintf('[%2d/%d] %s already built, skipping.\n', i, numel(allModels), allModels{i});
        continue
    end

    curr_core = RxnImp20(:,i);

    % weights apply to non-core reactions: higher weight, lower chance of
    % being pulled into the extracted model
    weights = zeros(numel(curr_core),1);
    weights(curr_core < 1) = 1 - curr_core(curr_core < 1);
    core = find(curr_core >= 1);

    fprintf('[%2d/%d] %s: %d core reactions -> spectraME ...\n', ...
        i, numel(allModels), allModels{i}, numel(core));
    t0 = tic;
    m = spectraME(model, core, 1e-4, [], weights);
    SP_times(i) = toc(t0);
    fprintf('        %d reactions, %d metabolites, %.1f s\n', ...
        numel(m.rxns), numel(m.mets), SP_times(i));

    save(outFile,'m');
end

save(fullfile(root,'models','build_times.mat'),'SP_times','allModels');
fprintf('All 20 models built.\n');
