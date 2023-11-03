import bioscience as bs

###################
# 1) Load dataset 
###################
# 1.1) Binary dataset load
#dataset = bs.load(path="/home/principalpc/git-repositories/benchmarkingBiBit/src/benchmarkingBiBit/datasets/binaryTest3.txt", index_gene=0, naFilter=False, head = 0)

# 1.2) Non-binary dataset load
#dataset = bs.load(path="/home/principalpc/git-repositories/benchmarkingBiBit/src/benchmarkingBiBit/datasets/synthetic3.txt", index_gene=0, naFilter=True, head = 0)

# 1.3.) RNA-Seq dataset load
#dataset = load(path="/home/principalpc/git-repositories/benchmarkingBiBit/src/benchmarkingBiBit/datasets/rnaseq.txt", index_gene=0, index_lengths=1 ,naFilter=True, head = 0)

###################
# 2) Preprocessing
###################
# 2.1) Standard preprocessing
#bs.discretize(dataset, n_bins= 2)
#bs.standardize(dataset)
#bs.scale(dataset)
#bs.normalDistributionQuantile(dataset)

# 2.2) RNA-Seq preprocessing
#bs.tpm(dataset)
#bs.cpm(dataset)

# 2.3) Binary preprocessing
#bs.binarize(dataset)
#listDatasets = bs.binarizeLevels(dataset, inactiveLevel = 0.2, activeLevel=0.8, soc = 0)
#listDatasets = bs.binarizeLevels(dataset)

###################
# 3) Data mining 
####################
# BiBit algorithm
# Single dataset
#listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=2, debug = True)

# List of datasets (bs.binarizeLevels function)
#listModels = bs.bibit(listDatasets, cMnr=2, cMnc=2, mode=2, debug = True)

###################
# 4) Save results
###################
# 4.1) Save genes
#bs.saveGenes(path="/home/principalpc/git-repositories/bioScience/", models=listModels, data=dataset) # Single dataset