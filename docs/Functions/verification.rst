Verification
============

An extraction promises that every core reaction is present and that the
result is flux consistent. ``checkExtraction`` independently verifies both
properties on any extracted model via FVA or ``spectraCC``.

Run it with the same arguments the extraction used, ``tol`` and
``useCoupling`` included: a check that answers a different question from
the extraction is not a check of it. Both of its methods are held to that
choice, which takes some doing on the ``fva`` side -- ``fluxVariability``
applies a model's coupling constraints whether or not you want it to, so
when ``useCoupling`` is false they are dropped from the copy it is given.

.. automodule:: src
   :members: checkExtraction
