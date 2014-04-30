.. Epitope_Mapping documentation master file, created by
   sphinx-quickstart on Sun Apr 27 11:09:01 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Epitope Mapping's documentation!
===========================================

Epitope Mapping is intended to predict several regions on an antigen for a known antibody-antigen pair
that are candidates for the true epitope of the antigen with respect to that antibody.

Epitope Mapping also contains methods to evaluate threshholds and values if the true epitope is known.
Thus, given a similar example, researchers could adjust default values and cutoffs to better tailor
the mapping software to their research needs. 

The classes involved in initial mapping include:

1.  PDBParse
2. SequentialPairs
3. ColorResidues
4. PDBUtil
5. MatrixReduce

The classes involved in assessing results at given threshholds include:

1. PrecisionRecall
2. Statistics
3. TestThreshold

=============================
Epitope Predition Contents:
=============================

.. toctree::
   :maxdepth: 2

.. automodule:: Driver

.. automodule:: PDBParse

.. autoclass:: PDBParse
    :members:

.. automodule:: SequentialPairs

.. autoclass:: SequentialPairs
    :members:

.. automodule:: PyMolColorSpheres

.. autoclass:: ColorResidues
    :members:

.. automodule:: PDBUtil

.. autoclass:: PDBUtil
    :members:

=============================
Method Assessment Contents:
=============================

.. automodule:: PrecisionRecall

.. autoclass:: PrecisionRecall
    :members:

.. automodule:: Statistics

.. autoclass:: Statistics
    :members:

.. automodule:: TestThreshhold

.. autoclass:: TestThreshhold
    :members:

.. automodule:: MatrixReduce

.. autoclass:: MatrixReduce
    :members:

======================
Indices and tables
======================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

