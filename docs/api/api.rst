API reference
=============
This page gives an overview of all public bioScience objects, functions and methods. All classes and functions exposed in ``bioscience.*`` namespace are public.

The following subpackages are public.

* ``bioscience.base``: It contains functions for the management of I/O operations, files and objects and methods that form the core of the library and are used by the rest of the library's subpackages.
* ``bioscience.preprocess``: This subpackage includes the functions associated with each of the preprocessing methods implemented in bioScience.
* ``bioscience.dataMining``: This subpackage contains all the source code for those data mining techniques that can be used in bioScience.

Within each subpackage there may be various functions and methods that may be private because they are often intermediate operations of the implemented methods.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   api
   

bioscience.base
---------------

.. automodule:: bioscience.base.files
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: bioscience.base.models
   :members:
   :undoc-members:
   :show-inheritance:

bioscience.preprocess
---------------------

.. automodule:: bioscience.preprocess.Standard
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: bioscience.preprocess.RnaSeq
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: bioscience.preprocess.Binarization
   :members:
   :undoc-members:
   :show-inheritance:

bioscience.dataMining
---------------------

Biclustering
^^^^^^^^^^^^

.. automodule:: bioscience.dataMining.biclustering.Biclustering
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: bioscience.dataMining.biclustering.BiBit
   :members:
   :undoc-members:
   :show-inheritance:

