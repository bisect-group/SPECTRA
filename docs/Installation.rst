Installation
============

Requirements
------------

- MATLAB
- `COBRA Toolbox <https://opencobra.github.io/cobratoolbox/stable/index.html>`_
- A solver. The LP formulations run on any solver the COBRA Toolbox
  supports, including the bundled glpk. The mixed-integer formulations
  (``minNetMILP``, ``tradeOff``) also work on glpk but are much faster with
  Gurobi or CPLEX.

Setting up
----------

.. code-block:: bash

   git clone https://github.com/bisect-group/SPECTRA.git

.. code-block:: matlab

   initCobraToolbox(false)
   addpath(genpath('path/to/SPECTRA/src'))

A note on solver tolerance
--------------------------

Before extracting at genome scale, set the solver's feasibility tolerance
to 1e-9, the lowest most solvers accept:

.. code-block:: matlab

   changeCobraSolverParams('LP', 'feasTol', 1e-9);

The extraction reads its answer off an LP whose mass balance holds only to
the solver's tolerance. At a looser setting that residual is enough to break
chains carrying flux of the order of ``tol``, which leaves core reactions
present in the result but unable to carry flux. ``checkExtraction`` will
tell you whether that has happened.
