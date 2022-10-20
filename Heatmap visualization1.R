#Reading an input cancer data file
data_file=read.csv("C:\\Users\\Admino\\Documents\\cancer genomics\\GSE169750_rawCounts.csv",sep=,header=T,row.names = 1)
print(data_file)

#Create a heatmap
library(ComplexHeatmap)
Heatmap(pmat)