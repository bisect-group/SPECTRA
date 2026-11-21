Minimal microbiomes
===================

A large community contains many organisms that contribute nothing to the
function you care about. A *minimal microbiome* is the smallest subset that
still performs it -- still grows, and still makes the metabolite you are
interested in -- at a stated fraction of what the whole community managed.

This is the extraction machinery pointed at organisms rather than
reactions, and the mapping is exact. Give each organism's anchor reaction a
binary through ``indicatorGroups`` and weight it 1, weight everything else
0, and minimising the weighted count *is* minimising the number of
organisms present. It runs through :doc:`extraction` with
``probType = 'minNetMILP'``, the same routine as every other reconstruction
in the package, with a different objective pointed at a different thing.

The biomass coupling built by ``buildCommunityModel`` is what makes an
absent organism truly absent: an organism whose binary is off cannot grow,
and coupling then pins every one of its reactions to zero, so it cannot go
on feeding the rest as a ghost. A community built without coupling is
rejected in spirit if not in code, and says so.

The membership is read off the binaries, not off the reduced model. The
distinction matters: an organism's anchor can survive on a trace of flux
with its binary off, so reading presence from the fluxes would count an
absent organism as present.

The requirements are measured before they are imposed. The full community
is solved for its growth, then for how much it can produce while each
organism keeps ``growthOptimumFraction`` of its own growth, and the minimal
community is held to a fraction of each. Asking for a fraction of something
unmeasured would make the answer depend on the units of the biomass
reaction.

Follows Raghu, A. K., Palanikumar, I. and Raman, K. (2024), Designing
function-specific minimal microbiomes from large microbial communities,
*npj Systems Biology and Applications* **10**, 46.

.. automodule:: src
   :members: minimalMicrobiome
