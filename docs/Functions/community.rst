Communities and bodies
======================

Joins models into one network whose units trade metabolites through a pool.
The units may be organisms of a microbial community or tissues of a body;
nothing in the construction distinguishes them.

Each unit keeps its own copy of the network, tagged so that its reactions
and metabolites cannot be confused with another's, and the units meet in a
**pool**: a compartment several of them can reach. One unit's secretion
arrives in the pool and is available to the rest, which is what makes
cross-feeding possible. A pool named in ``environment`` also carries
exchange reactions, which is how anything enters or leaves the whole system.

There are two ways to build the pool, and both are in published use:

``'pooled'``
   Each unit keeps a private external compartment, and a transport reaction
   moves each exchangeable metabolite between it and the pool. Because the
   transport belongs to one unit, it can be bounded per unit -- which is
   where measured uptake and secretion rates go. Several pools can coexist,
   so the units need not all meet in the same place.

``'shared'``
   The units share the external compartment outright: one row of the
   stoichiometric matrix per exchangeable metabolite, for the whole
   community. There are no transport reactions and there is only ever one
   pool. For a hundred units across a few thousand metabolites this saves a
   hundred thousand reactions, which is the difference between a model that
   can be built and one that cannot.

The compartmentalised layout, the biomass coupling and its parameters follow
the joint FBA formulation used for minimal microbiomes in Raghu et al.
(2024), *npj Systems Biology and Applications* **10**, 46.

What comes back is a structure describing the joined model, not the model
alone: which reactions belong to which unit, which anchor each unit is
coupled to, what each pool holds, and -- where databases were supplied --
which reactions a unit already had and which are candidates to be filled
in. A handful of functions read and reshape it:

``communityOrganismOf``
   Which unit a reaction or metabolite belongs to.

``couplingSummary``
   A table of what each unit is coupled to and how hard, which is where to
   look when a unit will not grow.

``communitySubset``
   The community restricted to some of its units, keeping the pools so
   that the survivors still trade and still face the same environment.

``communityWithModel``
   The same community over a reduced model, with every list narrowed to
   what survived. This is what a consistency check or an extraction leaves
   behind.

``communityDecompose``
   The community split back into one model per unit, untagged, with an
   exchange wherever the unit met a pool. This is how a joint result is
   cashed out into models that stand on their own.

.. automodule:: src
   :members: buildCommunityModel, buildMultiTissueModel, communityOrganismOf, couplingSummary, communitySubset, communityWithModel, communityDecompose
