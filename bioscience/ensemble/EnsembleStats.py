from bioscience.base import *
import bioscience as bs

import math
import sys
import os
import threading
import warnings
import numpy as np
import time

def ensembleStats(dataset, methods, thresholds, deviceCount = 0, mode = 1, debug = False):
    
    finalStats = []
    for index, methodValue in enumerate(methods):
        thresholdValue = thresholds[index]
        resultsStats = None
        
        if isinstance(methodValue, str):                     
            if methodValue.lower()=="kendall": # Kendall                
                resultsStats = bs.kendall(dataset)
                results = resultsStats.results
                
                # Transform correlation values to range from 0 to 1.
                for indexCorr, valueCorr in enumerate(results):
                    if valueCorr != None and not np.isnan(valueCorr):                        
                        if valueCorr < 0:
                            resultsStats.results[indexCorr] = resultsStats.results[indexCorr] * -1
                        
                        if resultsStats.results[indexCorr] < thresholdValue:
                            resultsStats.results[indexCorr] = None
                            resultsStats.geneInteractionsIndex[indexCorr] = None
                
            elif methodValue.lower()=="spearman": # Spearman              
                resultsStats = bs.spearman(dataset)
                results = resultsStats.results
                
                # Transform correlation values to range from 0 to 1.
                for indexCorr, valueCorr in enumerate(results):
                    if valueCorr != None and not np.isnan(valueCorr):                        
                        if valueCorr < 0:
                            resultsStats.results[indexCorr] = resultsStats.results[indexCorr] * -1
                        
                        if resultsStats.results[indexCorr] < thresholdValue:
                            resultsStats.results[indexCorr] = None
                            resultsStats.geneInteractionsIndex[indexCorr] = None                
                
            elif methodValue.lower()=="nmi":  # NMI              
                resultsStats = bs.nmi(dataset)                
                results = resultsStats.results
                
                # Transform correlation values to range from 0 to 1.
                for indexCorr, valueCorr in enumerate(results):
                    if valueCorr != None and not np.isnan(valueCorr):                        
                        if valueCorr < 0:
                            resultsStats.results[indexCorr] = resultsStats.results[indexCorr] * -1
                        
                        if resultsStats.results[indexCorr] < thresholdValue:
                            resultsStats.results[indexCorr] = None
                            resultsStats.geneInteractionsIndex[indexCorr] = None
            else:
                resultsStats = -1
        
        if resultsStats == None:
            print(f"Some problems have arisen in the execution of the {methodValue} method.")
        elif resultsStats == -1:
            print(f"Method {methodValue} not found.")
        else:
            finalStats.append(resultsStats)
    
    for index, resultsStats in enumerate(finalStats):
        print(f"Indice: {index} - Results: {resultsStats.results}")