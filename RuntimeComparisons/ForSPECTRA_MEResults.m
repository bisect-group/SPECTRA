clearvars -except rep seed
% rep and seed come from the caller. A plain `clear` wiped them, which is why
% this script could not previously be run as replicates.
if ~exist('rep','var'),  rep  = 1;    end
if ~exist('seed','var'), seed = 1000 + rep; end

% SEED EXPLICITLY. MATLAB starts every session from the same RNG state, so
% running this script in three separate processes would draw the IDENTICAL
% core sets and report zero spread - error bars that look real and are not.
rng(seed, 'twister');

S = load('./consRecon3D');
model = S.consRecon3D;   % the file stores consRecon3D, not model
% initCobraToolbox(false)
changeCobraSolver('gurobi','all')
% The ten core-set sizes the published figure plots, i.e. the points
% linspace(0,528,10) selects from the original 20:20:10600 sweep. Running
% only these keeps a replicate under two minutes.
n = [20 1180 2360 3540 4700 5880 7060 8220 9400 10580]; % number of core reactions
tol =1e-4; % tolerance level for fastcore and sprintcore

weights = ones(numel(model.rxns),1); %weights for swiftcore

% variables to store the reactions in context-specific model
sprxns={}; sp2rxns={}; sp3rxns={}; swRrxns={}; swWRrxns={}; fcrxns={};
% variabels to store the runtime
sptimes=[]; sp2times=[]; sp3times=[];swRtimes=[]; swWRtimes=[]; fctimes=[];
% variables to store the nLPs
spLPS=[]; sp2LPS=[];swRLPS=[]; swWRLPS=[]; fcLPS=[];
% to store the core reactions obtained by random-sampling at each iteration
cores={};

% This variable stores the information about the models that do not have
% core reactions (or) not consistent
ModelQuality = struct();
ModelQuality.name = {};
ModelQuality.fault = {};
for i=1:numel(n)
    core = sort(randsample(numel(model.rxns),n(i)));
    cores{i} = core;
    core_rxns = model.rxns(core);
    
    % fastcore
    tic
    % 5th output: COBRA's fastcore returns
    % [tissueModel, coreRxnBool, coreMetBool, coreCtrsBool], so the LP count
    % is appended as a fifth rather than displacing one of those - every
    % other caller keeps working. Added by the fastcore LP-counter patch.
    [m,~,~,~,LPs] = fastcore(model,core,tol);
    time = toc;
    fcrxns{i}=m.rxns;
    fctimes(i)=time;
    fcLPS(i)=LPs;
    model_name=['FC_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % spectraME
    tic
    [m,LPs] = spectraME(model,core,tol);
    time = toc;
    sprxns{i}=m.rxns;
    sptimes(i)=time;
    spLPS(i)=LPs;
    model_name=['SP_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % spectraME2
    tic
    [m,LPs] = spectraME2(model,core,tol);
    time = toc;
    sp2rxns{i}=m.rxns;
    sp2times(i)=time;
    sp2LPS(i)=LPs;
    model_name=['SP2_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % SP3 - the old sprintcore(...,'MILP') arm - is DISABLED.
    %
    % sprintcore is what spectraME was called before the rename, and the
    % argument lists are positionally identical, so the call would now read
    %     [m,~] = spectraME(model,core,tol,[],[],[],[],'MILP');
    % with 'MILP' landing on spectraME's 8th argument, probType.
    %
    % It stays off because this benchmark compares LP formulations: a MILP
    % runtime is not comparable with them, so the arm would be timed and
    % never plotted - which is exactly what happened in the original run,
    % where sp3times was computed but never written to
    % Results_SPECTRA_ME_Runtime.mat and never appears in
    % Runtime_SPECTRA_ME.ipynb. It is also unbounded in cost, being the only
    % MILP in the loop.
    
    % swiftcore with reduction
    tic
    [m,~,LPs] = swiftcore(model,core,weights,1e-10,1);
    time = toc;
    swRrxns{i}=m.rxns;
    swRtimes(i)=time;
    swRLPS(i)=LPs;
    model_name=['SWr_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
    
    % swiftcore without reduction
    tic
    [m,~,LPs] = swiftcore(model,core,weights,1e-10,0);
    time = toc;
    swWRrxns{i}=m.rxns;
    swWRtimes(i)=time;
    swWRLPS(i)=LPs;
    model_name=['SWwr_',num2str(i)];
    ModelQuality=checkModelQuality(ModelQuality,m,core_rxns,model_name,tol);
end
save(sprintf('Results_SPECTRA_ME_Runtime_rep%d', rep))