%% config_pairs.m
% The 10 cancer / non-cancer CCLE pairs analysed in the paper.
% Copied from BuildCancerModels/flux_results/models_name.m so that this
% review folder is self-contained.

contexts = {'Brain';'Kidney';'Liver';'Lung';'Eye';'Pancreas';'Prostate'; ...
            'Breast';'Ovary';'Biliary Tract'};

models_cancer = {'ACH-000376'   % Glioblastoma
                 'ACH-000649'   % Renal Clear Cell Carcinoma
                 'ACH-000739'   % Hepatoblastoma
                 'ACH-000681'   % Lung Adenocarcinoma
                 'ACH-001295'   % Retinoblastoma
                 'ACH-000108'   % Adenosquamous Carcinoma of the Pancreas
                 'ACH-000177'   % Prostate Small Cell Carcinoma
                 'ACH-000019'   % Invasive Breast Carcinoma
                 'ACH-000001'   % High-Grade Serous Ovarian Cancer
                 'ACH-001961'}; % Gallbladder Cancer

models_normal = {'ACH-003150'   % Immortalized Astrocytes
                 'ACH-003143'   % Immortalized Epithelial Cells, Renal
                 'ACH-003155'   % Immortalized Liver Cells
                 'ACH-000064'   % Immortalized Lung Cells
                 'ACH-002466'   % Immortalized Epithelial Cells, Retinal
                 'ACH-003161'   % Immortalized Pancreatic Stromal Cells
                 'ACH-001453'   % Hyperplasia, Prostate
                 'ACH-000642'   % Immortalized Breast Cells
                 'ACH-000494'   % Immortalized Ovarian Cells
                 'ACH-001854'}; % Immortalized Cholangiocyte, Highly Differentiated

allModels = [models_cancer; models_normal];
