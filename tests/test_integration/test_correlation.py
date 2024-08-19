import bioscience as bs
import numpy as np

###################
# 1) Load dataset 
###################
# 1.1) Binary dataset load
dataset = bs.load(path="/home/principalpc/git-repositories/bioscience/datasets/binaryTest3.txt", index_gene=0, naFilter=False, head = 0)

# 1.2) Non-binary dataset load
#dataset = bs.load(path="/home/principalpc/git-repositories/bioscience/datasets/synthetic3.txt", index_gene=0, naFilter=True, head = 0)

# 1.3.) RNA-Seq dataset load
#dataset = bs.load(path="/home/principalpc/git-repositories/bioscience/datasets/rnaseq.txt", index_gene=0, index_lengths=1 ,naFilter=True, head = 0)

###################
# 2) Preprocessing
###################
# 2.1) Standard preprocessing
#bs.discretize(dataset, n_bins= 2)
#bs.standardize(dataset)
#bs.scale(dataset)
#bs.normalDistributionQuantile(dataset)
#bs.outliers(dataset)

# 2.2) RNA-Seq preprocessing
#bs.tpm(dataset)
#bs.cpm(dataset)

# 2.3) Binary preprocessing
#bs.binarize(dataset, threshold=0.6)
#listDatasets = bs.binarizeLevels(dataset, inactiveLevel = 0.2, activeLevel=0.8, soc = 0)
#listDatasets = bs.binarizeLevels(dataset)

##################
# 3) Correlation 
##################
"""resultsCorrelation = bs.kendall(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.spearman(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.nmi(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.hoeffdingsD(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.pearson(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.mi(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.median(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.q(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.distcorr(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)

resultsCorrelation = bs.mcc(dataset)
print(resultsCorrelation.results)
#print(resultsCorrelation.geneInteractionsIndex)"""

#stats = bs.ensembleStats(dataset, methods=np.array(["kendall","spearman","hoeffdingsD","nmi","pearson","mi","median","q","distcorr","mcc"]),thresholds=np.array([0.5,0,0,0.1,0,0,0.6,0,0,0]))
#stats = bs.ensembleStats(dataset, methods="ordinal",thresholds=np.array([0.5,0,0]))
#stats = bs.ensembleStats(dataset, methods="ordinal",thresholds=0)

"""for index, resultsStats in enumerate(stats):
    print(f"Indice: {index} - Name: {resultsStats.name} - Results: {resultsStats.results}")"""

###################
# 3) Data mining 
####################
# BiBit algorithm
#listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=1, debug = True)

