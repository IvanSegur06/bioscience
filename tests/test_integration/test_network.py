import bioscience as bs

###################
# 1) Load dataset 
###################
dataset = bs.loadNetwork(path="C:/Users/pablo/OneDrive/Escritorio/Repositorio Investigación/bioscience/datasets/network1.csv", separator = ",", index_nodeA = 0, index_nodeB = 2, index_weight = 8, skipr = 0, head = 0)

print("Dataset data:")
print(dataset.data)