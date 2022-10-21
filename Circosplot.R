#Reading an input cancer data file
data_file=read.csv("C:\\Users\\Admino\\Documents\\cancer genomics\\GSE169750_rawCounts.csv",sep=,header=T,row.names = 1)
print(data_file)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RCircos")

library("RCircos") 

data("UCSC.HG38.Human.CytoBandIdeogram")
cyto.info = UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

cyto.info = UCSC.HG38.Human.CytoBandIdeogram
cyto.info$Name = NA
cyto.info$Stain = NA
RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)

chr_order = unique(cyto.info$Chromosome)

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

ideo = RCircos.Get.Plot.Ideogram()
ideo$BandColor = 'salmon'
num = which(ideo$Chromosome == 'chrX')
ideo[num, 'BandColor'] = 'chartreuse'

num = which(ideo$Chromosome == 'chrY')
ideo[num, 'BandColor'] = 'purple'


RCircos.Reset.Plot.Ideogram(ideo)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

num = which(ideo$Chromosome == 'chr1')
ideo[num, 'ChrColor'] = 'red'

RCircos.Reset.Plot.Ideogram(ideo)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()


library(biomaRt)

mat = readRDS("mat.rds")
print(mat)

m = useMart('ensembl', dataset='hsapiens_gene_ensembl')

coords = getBM(attributes=c('chromosome_name', 'start_position', 
                            'end_position', 'hgnc_symbol'),
               filters = c('hgnc_symbol'),
               values = list(rownames(mat)),
               mart = m)

write.csv(coords, file ='coords.csv')
coords$chromosome_name = paste0('chr', coords$chromosome_name)
coords$chromosome_name = factor(coords$chromosome_name, levels = chr_order)

num = which(is.na(coords$chromosome_name))
coords = coords[-num, ]

up = which((mat$pval < 0.01) &
             (mat$log2FC > 1))
upmat = mat[up, ]

num = which(coords$hgnc_symbol %in% rownames(upmat))
coords1 = coords[num, ]

RCircos.Gene.Name.Plot(coords1, name.col=4, track.num = 2, side = "in",
                       is.sorted = F)

