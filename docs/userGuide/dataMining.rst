Data mining
===========
In this section you will find a number of examples related to data mining techniques that can be used in bioScience. 

In order to deliver the results of these techniques in the shortest possible time and in the most efficient way, each data mining technique implemented in this library is adapted to a High-Performance Computing (HPC) environment. Therefore, each of these techniques have multiple execution modes, from a simple sequential execution, to a parallel CPU-based architecture or even multi-GPU (multi-GPU) architectures. 

Biclustering methods
^^^^^^^^^^^^^^^^^^^^
The following is a list of the Biclustering methods that are implemented in the bioScience library.

BiBit algorithm
---------------

BiBit algorithm was specifically designed to work with binary matrices, where the matrix elements can only take two values, usually 0 and 1. This is common in various fields like computational biology, where matrices may represent the presence or absence of genes across different experimental conditions.

The basic idea behind BiBit is to find maximal submatrices in a binary matrix where all elements are 1. In terms of biclustering, this translates into finding subgroups of rows that show consistent activity patterns (all 1s) across a subset of columns.

The resulting biclusters from BiBit are useful for identifying local patterns that do not necessarily apply to the entire data matrix. In bioinformatics, for example, this could represent a group of genes that are only expressed together under certain experimental conditions.

BiBit is particularly efficient for binary data due to its straightforward design and the discrete nature of the data. However, it cannot be directly applied to non-binary data without a prior binarization stage. Also, it may struggle with noise and inconsistencies in the data, as it looks for exact patterns of 1s.

For more information, see the research paper related to this algorithm: `Domingo S. Rodríguez-Baena et al. <https://academic.oup.com/bioinformatics/article/27/19/2738/231788?login=false>`_

Note that the dataset must be binary for the execution of this algorithm. Therefore, the ``dataset`` object or the ``listDatasets`` list obtained must be previously binarised. In order to binarise a dataset, please refer to the :doc:`preprocessing section <preprocess>` of this user guide.

An example of how to **run this algorithm sequentially** is shown below:

    .. code-block:: python
      
        import bioscience as bs
        
        # Single dataset
        listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=1, debug = True)

        # List of datasets (if bs.binarizeLevels function is used)
        listModels = bs.bibit(listDatasets, cMnr=2, cMnc=2, mode=1, debug = True)

If you want to run this algorithm on a **parallel CPU architecture**, use the following source code:

    .. code-block:: python
      
        import bioscience as bs
        
        # Single dataset
        listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=2, debug = True)

        # List of datasets (if bs.binarizeLevels function is used)
        listModels = bs.bibit(listDatasets, cMnr=2, cMnc=2, mode=2, debug = True)

If the algorithm is executed in a **parallel environment with multiple GPU devices (multi-GPU)**, the following source code can be used:

    .. code-block:: python
      
        import bioscience as bs
        
        # Single dataset
        listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=3, deviceCount=1, debug = True)

        # List of datasets (if bs.binarizeLevels function is used)
        listModels = bs.bibit(listDatasets, cMnr=2, cMnc=2, mode=3, deviceCount=1, debug = True)

To understand the meaning of each attribute you can access the :doc:`API reference <../api/api>`.