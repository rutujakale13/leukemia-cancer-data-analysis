#Reading an input cancer data file
data_file=read.csv("C:\\Users\\Admino\\Documents\\cancer genomics\\GSE169750_rawCounts.csv",sep=,header=T,row.names = 1)
print(data_file)

data=source("data normalization.R")


mat=matrix(NA,ncol=4,nrow = nrow(logcpm))
rownames(mat)=rownames(logcpm)
colnames(mat)=c('meanTumor','meanControl','pvalue','log2FC')

for(i in 1:nrow(logcpm)){
  vector1 = as.numeric(logcpm[i, 1:7])
  vector2 = as.numeric(logcpm[i, 8:11])
  
  res=t.test(vector1, vector2, paired = F, alternative = "two.sided")
  mat[i,1]=res$estimate[[1]]
  mat[i,2]=res$estimate[[2]]
  mat[i,3]=res$p.value
  mat[i,4]=mat[i,1]-mat[i,2]
  
}

mat=as.data.frame(mat)
num=which(is.nan(mat$pvalue))
mat[num,'pvalue']=1

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)
EnhancedVolcano(mat,lab = rownames(mat),x = 'log2FC' ,y ='pvalue')