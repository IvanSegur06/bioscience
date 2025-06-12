from .Cobinet import *
from bioscience.base import *

def cobinet(dataset, deviceCount = 1, mode = 1, debug = False):
    listModels = set()

    if dataset is not None:
        
        if isinstance(dataset, Dataset):
            oModel = processCobinet(dataset, deviceCount, mode, debug)
            listModels.add(oModel)
        
        if isinstance(dataset, set):
            iLevel = 1
            for oDataset in dataset:
                print("\nLevel: ", str(iLevel))
                oModel = processCobinet(oDataset, deviceCount, mode, debug)
                listModels.add(oModel)
                iLevel += 1
            
    return listModels