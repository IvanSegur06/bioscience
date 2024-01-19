import bioscience as bs

###################
# 1) Load dataset 
###################
# Microarray (GSE26910)
dataset = bs.load(path="/home/principalpc/git-repositories/bioscience/datasets/bibit/GSE26910.csv", index_gene=0, naFilter=False, head = 0)

###################
# 2) Preprocessing
###################

bs.scale(dataset)

# Binarize
bs.binarize(dataset, threshold = 0.2)

# Handling outliers
#bs.outliers(dataset)

###################
# 3) Data mining 
####################
# BiBit algorithm Sequential
#listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=1, debug = True)

# BiBit algorithm CPU Parallel
listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=2, debug = True)

# BiBit algorithm GPU Parallel
#listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=3, deviceCount=1, debug = True)