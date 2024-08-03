Stats
=====
This section shows several examples of the many methods that can be used in bioscience to perform statistical applications.

Suppose you have an object called ``dataset`` with your dataset loaded. To see how a data load is performed, please refer to the :doc:`Load data section <load>` of this user guide.

In these cases, the pre-processed dataset will be found modified in the ``dataset.data`` attribute, while the originally loaded dataset will be found stored in the ``dataset.original`` attribute. This is because the user may wish to display or use the original dataset data for further validation or visualisation processes.

Correlation methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following examples show how to use correlation methods such as Kendall, Spearman, NMI, among others.
    
    .. code-block:: python
      
        import bioscience as bs
        resultsCorrelation = bs.kendall(dataset)
        resultsCorrelation = bs.spearman(dataset)
        resultsCorrelation = bs.nmi(dataset)
    
To understand the meaning of each attribute you can access the :doc:`API reference <../api/api>`.