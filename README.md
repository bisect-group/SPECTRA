<h1 align="center">
<table border="0">
  <tr>
    <td>
      <img src="images/SPECTRA_logo.png" width="80" height="80">
    </td>
    <td>
      <div align="center">
      <b>SPECTRA: Scalable Platform for Extracting</b><br>
      <b>Constraint-based Top-down Reconstructions and Analysis</b>
      </div>
    </td>
  </tr>
</table>
</h1>

# Overview

SPECTRA is a unified, scalable framework for reconstructing metabolic models from multi-omics data. It integrates diverse optimization strategies (linear and mixed-integer linear programming) to extract metabolic networks at multiple biological scales-from minimal reactomes and context-specific models (CSMs) to microbial community models and multi-tissue models. 

# Installation & Setup

## Hardware requirements
`SPECTRA` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
The package has been tested on the following systems:
+ Windows: 10
+ Linux: Ubuntu 22.04

### Prerequisites
- MATLAB
- [COBRA toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html)
- CPLEX or Gurobi (for optimization)

All simulations are implemented through the COBRA Toolbox 3.0 interface and run against Gurobi or IBM CPLEX (with GLPK supported for small models and unit tests). By default, SPECTRA passes `1e-9` feasibility, optimality, and (for MILPs) integrality tolerances per solve without modifying global COBRA solver parameters; pass `solverTol = 'default'` to use the active COBRA session defaults instead.

### Installation

```bash
git clone https://github.com/bisect-group/SPECTRA.git
cd SPECTRA
```

# Reproducing the manuscript

This branch holds the SPECTRA source code. The scripts, models and data behind
the manuscript — the cancer models, the Recon3D+ work, the runtime comparisons,
the AGORA gap-filling, the community models and the figures — live on the
[`paper_codes`](https://github.com/bisect-group/SPECTRA/tree/paper_codes)
branch:

```bash
git checkout paper_codes
```

# Quick Start Guide

## Step 1: Prepare The Inputs

The omics data can be used to derive the following three inputs to SPECTRA.

| Parameter | Purpose |
|---|---|
| **Core reactions** | These reactions will be present in the final model |
| **Weights** | Numbers that determine the relative importance of the reactions to be included (or excluded) |
| **Box constraints** | Lower and upper bounds to the reactions in the universal model |

The following code demonstrates the same for a toy network model:

```MATLAB
  % loading the universal model
  model = three_pathway_toy_model(); % A toy model with eight metabolites and eleven reactions

  % defining the box constraints 
  model.ub(:) = 10; 
  model.ub([1,6,9]) = 1; % constraining the media reactions

  % defining the core reactions
  core = 5; % the reaction indices of the core reactions

  % defining the weights of the reactions (For this example, we are using tradeOff as the network inference method)
  weights = [1,1,-2,-1,1,1,2,-1,-2,0,1]'; % Same size as the number of reactions
```

## Step 2: Choose The Network Inference Formulation

| Formulation | Objective | Use Case | Input Weights |
|---|---|---|---|
| **minNetLP** | Minimize weighted flux (`"With-Chain"` + IRLS) | Fast LP extraction preserving core reactions | Non-negative weights |
| **minNetMILP** | Minimize reaction count | Exact minimal networks preserving core reactions | Non-negative weights |
| **tradeOff** | Maximize net reaction weights | Balancing positive rewards vs. negative penalties | Any signed real values |
| **optimBiomass** | Maximize biomass − weighted flux | Biological growth optimization | Non-negative weights |

## Step 3: Model reconstruction

```matlab
% For flux consistent universal models
consType = 'stoichiometry'; % 'stoichiometry' (S*v = 0) or 'topology' (S*v >= 0)
probType = 'tradeOff'; % 'minNetLP' (default), 'minNetMILP', 'tradeOff', or 'optimBiomass'
altSolMethod = []; % defaults to 'coreDirection' for LP and 'pathwayExclusion' for MILP
nSol = 1; % number of alternative solutions required
solveTime = 7200; % maximum time (s) for solving MILP problem
tol = 1e-5; % minimum absolute flux threshold for core reactions
[Model_tradeOff] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); 

% For flux inconsistent universal models (simultaneous consistency checking + extraction)
[Model_tradeOff, BlockedCoreRxns] = spectraCCME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); 
```

## Step 4: To Generate Alternative Solutions

```matlab
  % Pass [] to use the formulation's default ('pathwayExclusion' for MILP
  % formulations 'tradeOff'/'minNetMILP', 'coreDirection' for LP formulations
  % 'minNetLP'/'optimBiomass'), or specify altSolMethod explicitly:
  altSolMethod = 'pathwayExclusion';
  nSol = 5; % number of alternative solutions required
  [Model_tradeoff_altSols] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); 
```

# Beyond context-specific models

The same machinery answers several other questions. What changes between
them is the universal model, the evidence supplied and the objective
chosen, not the routine called.

| Question | Function |
|---|---|
| Which reactions of this model can carry flux at all? | `spectraCC` |
| What is the smallest network that keeps my core reactions? | `spectraME`, `spectraCCME` |
| Can this model do the things the organism can? | `checkTask`, `checkTasks`, `metabolicTask` |
| What must I add so the draft grows in each of these media? | `gapfillForGrowth` |
| What must I add so it can perform these tasks? | `gapfillForTasks` |
| What is the fewest reactions that still grow and perform the tasks? | `minimalReactome` |
| Join several organisms into one model that cross-feeds | `buildCommunityModel` |
| Join the tissues of a body, sharing blood and a lumen | `buildMultiTissueModel` |
| Gap-fill every member of a community at once | `gapfillCommunity` |
| Which members of this community are actually needed? | `minimalMicrobiome` |
| Did the extraction do what it was asked? | `checkExtraction` |

A model may carry constraints that are not in its stoichiometric matrix:
the `C`, `d` and `dsense` fields that `coupleRxnList2Rxn` writes, and that
`buildCommunityModel` uses to tie every reaction of a unit to that unit's
anchor. SPECTRA builds its problems from `S`, `lb` and `ub` alone, so those
play no part unless you pass `useCoupling` — flux consistency is
conventionally defined over the stoichiometry, which is what a comparison
against fastcc or FASTERCC has to answer too. Turn it on when the
parameterisation is the point, or when the result is going straight into
FBA, which applies them regardless. `gapfillCommunity` and
`minimalMicrobiome` turn it on themselves.

A community comes back as a structure describing the joined model, which
`communityOrganismOf`, `couplingSummary`, `communitySubset`,
`communityWithModel` and `communityDecompose` read and reshape;
`communityDecompose` is what turns a joint result back into models that
stand on their own.

Every one of these has a tutorial as a MATLAB live script under
[`docs/Tutorials`](docs/Tutorials): open the `.mlx` in MATLAB and step
through it. The reference documentation is under [`docs`](docs), and the
test suite is run with `runSpectraTests` from [`test`](test).

SPECTRA is also available as a Python package for cobrapy,
[spectra-cobra](https://github.com/bisect-group/spectra-cobra).
