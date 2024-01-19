Single-Cell: Real use-case
==========================

The aim of this real use-case is to illustrate the power of bioScience against a real data set.

----

Dataset information
^^^^^^^^^^^^^^^^^^^
For this real case, the `"GSE246622" <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246622>` dataset from the `"NCBI GEO (Gene Expression Omnibus)" <https://www.ncbi.nlm.nih.gov/geo/>` database is used. This Single-cell dataset for the organism Homo sapiens has a total of 22593 rows and 533 samples.

To generate this dataset, healthy infants and infants with mild (RSV-infected, non-hospitalized), moderate (RSV-infected, hospitalized without mechanical ventilation), and severe (RSV-infected, hospitalized with mechanical ventilation) RSV disease were recruited during the acute RSV infection and at convalescence. Whole blood samples were analysed using Clariom GOScreen human assay in 384 well plate format. Infants with comorbidities excluded from the downstream comparative analyses.

While the majority of infants infected with respiratory syncytial virus (RSV) exhibit mild or no symptoms, approximately 3 million children under the age of five are hospitalized every year due to complications from RSV. This research sought to explore the biological processes and related biomarkers responsible for the varied manifestations of RSV disease in young infants. The goal is to pave the way for a more precise categorization of RSV-infected infants based on their medical requirements. Whole blood samples are collected from infant case-control cohort study, the RESCEU case-control cohort is a multinational, multicenter, observational study (clinical trial registration number: NCT03756766). Infants < 12 months old with RSV disease were recruited from the University Medical Center Utrecht (UMCU) in The Netherlands, Hospital Clínico Universitario de Santiago (SERGAS) in Spain, Imperial College (IMPERIAL) National Health Service Trust (NHS) and Oxford University Hospital NHS Trust (OXFORD) in the United Kingdom during the RSV seasons 2017-2018, 2018-2019, and 2019-2020. Healthy controls without underlying comorbidities were recruited outside of the RSV season. Eligibility criteria included hospitalization for less than 48 hours at enrolment or within 96 hours of disease onset, no previous receipt of medications to treat RSV infection, no prior exposure to an investigational RSV vaccine or medication, no previous receipt of immunoglobulins or monoclonal antibodies, and had not used montelukast or oral steroids within seven days before enrolment. Infants with co-morbidities were not evaluated in the manuscript. RSV was detected using RSV point-of-care test (POCT) by either a rapid antigen detection test (Alere I) (Alere Inc, Waltham, Massachusetts) or rapid RSV polymerase chain reaction (PCR) test at the hospital setting, or a RSV PCR test at the laboratory. Convalescence samples were collected 7 ±1 weeks after a positive RSV diagnostic test result. We used microarray to assist us to identify biomarkers for severe RSV disease.


Biclustering analysis
^^^^^^^^^^^^^^^^^^^^^
The aim of this analysis is to run a biclustering algorithm called BiBit to obtain those biclusters in which certain genes in the dataset show a similar pattern of behaviour in their gene expression values against a subset of samples.

Once this Biclustering algorithm has been executed, the execution times of each of the versions will be compared to demonstrate the power and usefulness of this library.

**1. Load single-cell dataset:** Load the single-cell dataset from the input file.

    .. code-block:: python

      import bioscience as bs

      dataset = bs.load(path="/home/user/Desktop/GSE246622.csv", index_lengths = 0, index_gene=0, naFilter=False, head = 0, separator=";")

**2. Preprocessing:** For this dataset it is not necessary to normalise the data as they are pre-processed before log2. However, this dataset must be binarised because the Biclustering algorithm used requires the dataset to contain binary values.

    .. code-block:: python

      import bioscience as bs

      bs.binarize(dataset)

**3. Biclustering:** In this phase, the BiBit Biclustering algorithm is run to generate useful knowledge from this dataset. This algorithm can be run in three modes: sequential, parallel CPU and multi-GPU.

For the sequential mode, the BiBit algorithm is processed **sequentially on the CPU processors** and for this mode to be executed a value of 1 must be specified for the ``mode parameter`` of the ``bibit function``. 

    .. code-block:: python

      import bioscience as bs

      listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=1)

If the user wishes to process the algorithm in **parallel on CPU processors**, the ``mode parameter`` must contain a value of 2.

    .. code-block:: python

      import bioscience as bs

      listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=2)

If the user wishes to execute it in a **multi-GPU environment**, the ``mode parameter`` must contain a value of 3. In addition, for environments with GPU devices, there is a ``deviceCount parameter`` in which the user indicates the number of GPU devices to be used for processing.

    .. code-block:: python

      import bioscience as bs

      listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=3, deviceCount=2)

Once the Biclustering algorithm is executed, it is detected that a total of 149 biclusters are generated.


**4. Results:** It is possible to save the name of the genes of each bicluster generated by BiBit:

   .. code-block:: python
      
      bs.saveGenes(path="/home/user/Desktop/", models=listModels, data=dataset) # Single dataset


Execution times
^^^^^^^^^^^^^^^
The execution times of each of the versions are compared to demonstrate the power and usefulness of this library.

This experiment was conducted on a system equipped with an Intel Xeon E5-2686 v4 processor featuring 18 cores operating at 2.30 GHz, 32 GB of RAM, and 8 NVIDIA K80 12 GB graphics cards, each offering a combined total of 2496 CUDA cores.

  ..  csv-table:: Numbers
    :header: "Sequential mode", "CPU Parallel", "GPU Parallel (1 GPU)", "GPU Parallel (2 GPU)"
    :widths: 25, 25, 25, 25

    "20431,81 s.","2667,31 s.","674,46 s.","334,12 s."

The run times shown in the table above are in seconds. As can be seen, the interest in the use of High Performance Computing (HPC) in the field of Bioinformatics is gaining more and more relevance due to the increasing volume of datasets and the complexity of data mining techniques to extract useful knowledge.
