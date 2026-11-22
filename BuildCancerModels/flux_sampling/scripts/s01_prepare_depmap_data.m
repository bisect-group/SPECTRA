%% s01_prepare_depmap_data.m
% Reviewer-2 response, Phase 1a.
%
% Rebuilds DepMapData.mat, the input to the context-specific model
% reconstruction. The original BuildCancerModels/DepMapModels_Recon/
% PrepareDepMapData.m read 'Expression_Public_23Q4.csv' (the DepMap portal
% custom-download layout, gene symbols as bare column headers). That file is
% not archived, so we start from the equivalent figshare release of the same
% DepMap 23Q4 data, OmicsExpressionProteinCodingGenesTPMLogp1.csv, whose
% headers are 'SYMBOL (ENTREZ)'. The entrez suffix is stripped to recover the
% symbols; everything else (2^x - 1 to undo the log2(TPM+1)) is unchanged.

clear
root = '/home/pavan/scratch/SPECTRA_review_2';
csvFile = fullfile(root,'data','OmicsExpressionProteinCodingGenesTPMLogp1.csv');
outFile = fullfile(root,'data','DepMapData.mat');

if exist(outFile,'file')
    fprintf('%s already exists, skipping.\n', outFile);
    return
end

fprintf('Reading %s ...\n', csvFile);
t0 = tic;
tbl = readtable(csvFile,'VariableNamingRule','preserve');
fprintf('  read in %.1f s: %d rows x %d cols\n', toc(t0), height(tbl), width(tbl));

% first column holds the DepMap ModelID (ACH-######) and has an empty header
tbl.Properties.VariableNames{1} = 'Samples';

geneExpression.value   = table2array(tbl(:,2:end));
geneExpression.value   = (2.^geneExpression.value) - 1;   % undo log2(TPM+1)
geneExpression.context = tbl.Samples;

% 'TSPAN6 (7105)' -> 'TSPAN6'
rawNames = tbl.Properties.VariableNames(2:end)';
geneExpression.genes = regexprep(rawNames,'\s*\(\d+\)\s*$','');
geneExpression.genes = strtrim(geneExpression.genes);

fprintf('genes: %d (%d unique)\n', numel(geneExpression.genes), ...
    numel(unique(geneExpression.genes)));
fprintf('contexts: %d\n', numel(geneExpression.context));

% sanity check against the model the paper used
S = load(fullfile(root,'data','consRecon3DGeneSymbol.mat'));
mdl = S.model;
nOverlap = sum(ismember(geneExpression.genes, mdl.genes));
fprintf('model genes: %d | present in expression matrix: %d (%.1f%%)\n', ...
    numel(mdl.genes), nOverlap, 100*nOverlap/numel(mdl.genes));
if nOverlap < 0.5*numel(mdl.genes)
    error(['Only %d of %d model genes matched the expression matrix. The gene ' ...
           'identifier convention does not line up; do not proceed.'], ...
           nOverlap, numel(mdl.genes));
end

save(outFile,'geneExpression','-v7.3');
fprintf('Saved %s\n', outFile);
