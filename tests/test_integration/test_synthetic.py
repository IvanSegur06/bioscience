import bioscience as bs

###################
# 1) Load dataset 
###################
# Binary synthetic dataset with 4000 rows and 4000 columns
dataset = bs.load(path="/home/principalpc/git-repositories/bioscience/datasets/bibit/synthetic.csv", separator=",")

####################
# 3) Data mining 
####################
# BiBit algorithm Sequential
#listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=1, debug = True)

# BiBit algorithm CPU Parallel
#listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=2, debug = True)

# BiBit algorithm GPU Parallel
listModels = bs.bibit(dataset, cMnr=2, cMnc=2, mode=3, deviceCount=1, debug = True)