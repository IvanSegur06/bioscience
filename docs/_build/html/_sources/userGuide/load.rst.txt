Load dataset
============
Depending on the type of dataset you have, bioScience offers the possibility to offer different types of dataset loading. 

Load gene co-expression dataset (microarray, synthetic and generic)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This option allows the loading of a microarray, synthetic or generic dataset.

    .. code-block:: python

        import bioscience as bs
        dataset = bs.load(path="datasets/synthetic.txt", index_gene=0, naFilter=True, head = 0)

To understand the meaning of each attribute in this load function you can access the :doc:`API reference <../api/api>`.


Load RNA-Seq dataset 
^^^^^^^^^^^^^^^^^^^^
The following source code allows the loading of a dataset of type RNA-Seq.

    .. code-block:: python

        import bioscience as bs
        dataset = load(path="datasets/rnaseq.txt", index_gene=0, index_lengths=1 ,naFilter=True, head = 0)

To understand the meaning of each attribute in this load function you can access the :doc:`API reference <../api/api>`.


Load binary dataset 
^^^^^^^^^^^^^^^^^^^^
bioScience also allows loading of binary datasets because certain data mining algorithms only support this type of data. To do so, the user can either perform a direct load of a previously binarised dataset or load his dataset and binarise it internally with the bioScience library.

    .. code-block:: python

        import bioscience as bs
        dataset = bs.load(path="datasets/binary.txt", index_gene=0, naFilter=False, head = 0)

To understand the meaning of each attribute in this load function you can access the :doc:`API reference <../api/api>`.