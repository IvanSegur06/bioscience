Pre-processing
==============
This section shows several examples of the many methods that can be used in bioScience to perform data pre-processing. 

Suppose you have an object called ``dataset`` with your dataset loaded. To see how a data load is performed, please refer to the :doc:`Load data section <load>` of this user guide.

In these cases, the pre-processed dataset will be found modified in the ``dataset.data`` attribute, while the originally loaded dataset will be found stored in the ``dataset.original`` attribute. This is because the user may wish to display or use the original dataset data for further validation or visualisation processes.

Generic pre-processing methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following example shows different basic pre-processing options such as discretisations, standardisations, scaling and normal distributions based on quantiles. In addition, outlier treatment is possible.
    
    .. code-block:: python
      
        import bioscience as bs
        bs.discretize(dataset, n_bins= 2)
        bs.standardize(dataset)
        bs.scale(dataset)
        bs.normalDistributionQuantile(dataset)
        bs.outliers(dataset)
    
To understand the meaning of each attribute you can access the :doc:`API reference <../api/api>`.

RNA-Seq oriented pre-processing methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This subsection shows how to perform preprocessing with specific methods such as CPM, TPM, FPKM and DESEq2. 

    .. code-block:: python
    
        import bioscience as bs
        bs.tpm(dataset)
        bs.cpm(dataset)
        bs.fpkm(dataset)
        bs.deseq2Norm(dataset)

To understand the meaning of each attribute you can access the :doc:`API reference <../api/api>`.

Binarisation methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Currently, there are two ways to binarise a dataset. The ``binarize`` function performs a standard binarisation of the dataset, while the ``binarizeLevels`` function gets a list of binarised datasets. The latter function uses fuzzy logic to avoid noise that may be incorporated into the data in the binarisation process. 

Different examples of binarisation are shown below:

    .. code-block:: python

        import bioscience as bs
        bs.binarize(dataset)
        listDatasets = bs.binarizeLevels(dataset, inactiveLevel = 0.2, activeLevel=0.8, soc = 0)
        listDatasets = bs.binarizeLevels(dataset)

To understand the meaning of each attribute you can access the :doc:`API reference <../api/api>`.