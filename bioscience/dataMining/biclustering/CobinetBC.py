from bioscience.base import *


def processCobinetBC(dataset, deviceCount=1, mode=1, debug=False):
    """
    Main processing function for the CoBiNet biclustering algorithm.
    Follows the execution pattern established in BCCA:
      - mode=1: Sequential execution (implemented)
      - mode=2: CPU-parallel execution (to be developed)
      - mode=3: GPU-parallel execution (to be developed)

    :param dataset: The dataset object storing the data of the input file.
    :type dataset: :class:`bioscience.base.models.Dataset`
    :param deviceCount: Number of GPU devices to execute, defaults to 1.
    :type deviceCount: int, optional
    :param mode: Type of execution of the algorithm.
                 1 = sequential, 2 = CPU parallel, 3 = GPU parallel.
    :type mode: int, optional
    :param debug: Run the algorithm in debug mode, defaults to False.
    :type debug: bool, optional

    :return: A BiclusteringModel object containing the biclusters.
    :rtype: :class:`bioscience.base.models.BiclusteringModel`
    """
    oModel = None

    sMode = ""
    if mode == 2:  # NUMBA: CPU Parallel mode
        # To be developed
        sMode = "NUMBA - CPU Parallel mode (to be developed)"
    elif mode == 3:  # NUMBA: GPU Parallel mode
        # To be developed
        sMode = "NUMBA - GPU Parallel mode (to be developed)"
    else:  # Sequential mode
        oModel = __cobinetSequential(dataset, debug)
        deviceCount = 0
        sMode = "CPU Sequential"

    return oModel


def __cobinetSequential(dataset, debug=False):


    print("Hello world")

    oModel = BiclusteringModel()
    return oModel
