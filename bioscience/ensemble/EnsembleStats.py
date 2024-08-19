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
    # 1) Select methods based on categories. The category is based on the type of data (ordinal, continuous, dichotomous, mixed, distance), the nature of the data (parametric, non-parametric) or the type of relationship (linear, non-linear, rank-based, mixed, similarity/distance).
    if isinstance(methods, str):        
        if methods == "ordinal":
            methods = np.array(["kendall","spearman","hoeffdingsD","wrc","tabwilr"])
        elif methods == "continuous":
            methods = np.array(["pearson","distcorr","mcd","q","median","pc","biweight","cca"])
        elif methods == "dicotomic":
            methods = np.array(["mcc","pbc","phi","log-odds"])
        elif methods == "mixed":
            methods = np.array(["mi","nmi","mifs","chi-squared","adjusted-rand","contigency"])
        elif methods == "distance":
            methods == np.array(["euclidean","manhattan","mahalanobis","jaccard","cos","taba","tabwil","tabwilr"])
        elif methods == "parametric":
            methods = np.array(["pearson","phi","cca","log-odds","pc","euclidean","mahalanobis","chi-squared"])
        elif methods == "non-parametric":
            methods = np.array(["kendall","spearman","tabwilr","hoeffdingsD","distcorr","mcd","q","median","mi","nmi","mifs","biweight","adjusted-rand","contingency","manhattan","cos","jaccard"])
        elif methods == "lineal":
            methods = np.array(["pearson","phi","pc","cca","log-odds"])
        elif methods == "non-lineal":
            methods = np.array(["distcorr","mi","nmi","mifs","hoeffdingsD"])
        elif methods == "ranks":
            methods = np.array(["kendall","spearman","tabwilr","wrc"])
        elif methods == "dataMixed":
            methods = np.array(["chi-squared","adjusted-rand","contingency"])
        elif methods == "similarity":
            methods = np.array(["euclidean","manhattan","mahalanobis","cos","jaccard"])
        else:
            methods = np.array(["kendall","spearman","nmi"])
            
    if isinstance(thresholds, float) or isinstance(thresholds, int):
        thresholds = np.full(methods.shape[0], thresholds)
    
    
    if isinstance(methods,(np.ndarray, np.generic)) and isinstance(thresholds,(np.ndarray, np.generic)):        
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
                
                elif methodValue.lower()=="hoeffdingsd": 
                    resultsStats = bs.hoeffdingsD(dataset)
                    results = resultsStats.results
                    
                    # Replace values between -0.5 and 1.
                    resultsStats.results = np.clip(results, -0.5, 1)
                    
                    # Transform correlation values to range from 0 to 1.
                    for indexCorr, valueCorr in enumerate(results):
                        if valueCorr != None and not np.isnan(valueCorr):   
                            resultsStats.results[indexCorr] = (resultsStats.results[indexCorr] + 0.5) / 1.5
                                                    
                            if resultsStats.results[indexCorr] < thresholdValue:
                                    resultsStats.results[indexCorr] = None
                                    resultsStats.geneInteractionsIndex[indexCorr] = None
                
                elif methodValue.lower()=="pearson":  # NMI              
                    resultsStats = bs.pearson(dataset)                
                    results = resultsStats.results
                    
                    # Transform correlation values to range from 0 to 1.
                    for indexCorr, valueCorr in enumerate(results):
                        if valueCorr != None and not np.isnan(valueCorr):                        
                            if valueCorr < 0:
                                resultsStats.results[indexCorr] = resultsStats.results[indexCorr] * -1
                            
                            if resultsStats.results[indexCorr] < thresholdValue:
                                resultsStats.results[indexCorr] = None
                                resultsStats.geneInteractionsIndex[indexCorr] = None
                
                elif methodValue.lower()=="mi":  # MI              
                    resultsStats = bs.mi(dataset)                
                    results = resultsStats.results
                    
                    # Transform correlation values to range from 0 to 1.
                    for indexCorr, valueCorr in enumerate(results):
                        if valueCorr != None and not np.isnan(valueCorr):                        
                            if valueCorr < 0:
                                resultsStats.results[indexCorr] = resultsStats.results[indexCorr] * -1
                            
                            if resultsStats.results[indexCorr] < thresholdValue:
                                resultsStats.results[indexCorr] = None
                                resultsStats.geneInteractionsIndex[indexCorr] = None
                                
                elif methodValue.lower()=="median":  # Median
                    resultsStats = bs.median(dataset)                
                    results = resultsStats.results
                    
                    minValue = np.min(results)
                    maxValue = np.max(results)
                    
                    # Transform correlation values to range from 0 to 1.
                    for indexCorr, valueCorr in enumerate(results):
                        if valueCorr != None and not np.isnan(valueCorr):
                                                        
                            normValue = (valueCorr - minValue) / (maxValue - minValue)
                            
                            if normValue < thresholdValue:
                                resultsStats.results[indexCorr] = None
                                resultsStats.geneInteractionsIndex[indexCorr] = None
                
                elif methodValue.lower()=="q":  # Quadrant
                    resultsStats = bs.q(dataset)                
                    results = resultsStats.results
                    
                    # Transform correlation values to range from 0 to 1.
                    for indexCorr, valueCorr in enumerate(results):
                        if valueCorr != None and not np.isnan(valueCorr):
                                                        
                            if valueCorr < 0:
                                resultsStats.results[indexCorr] = resultsStats.results[indexCorr] * -1
                            
                            if resultsStats.results[indexCorr] < thresholdValue:
                                resultsStats.results[indexCorr] = None
                                resultsStats.geneInteractionsIndex[indexCorr] = None
                
                elif methodValue.lower()=="distcorr":  # Distance Correlation (distcorr)
                    resultsStats = bs.distcorr(dataset)                
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
        
        return np.array(finalStats)