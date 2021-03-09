library(dplyr)
library(gplots)



source("./utils/quy_TCR_utils.R")

tcr.res <- read.csv("filtered_contig_annotations.csv", stringsAsFactors = FALSE)

tcr.res.combined <- TCRprofle(tcr.res)

# separate by each sample
tcr.res.combined.freq <- calClonalFreq(tcr.res.combined, cloneType="TCRgene")

# compute TRA TRB gene composition
# separate by each sample
tcr.barcode <- tcr.res.combined[,c("barcode","TCRgene")]
tcr.gene.matrix <- TCRgeneCompMatrix(tcr.barcode)
tcr.gene.matrix <- tcr.gene.matrix*100/nrow(tcr.barcode)


bk = unique(c(seq(0, 20, length=30)))
col = colorRampPalette(c("white","#162955"))(length(bk)-1)
pdf("heatmap_of_TCR_gene_composition.pdf")
heatmap.2(as.matrix(tcr.gene.matrix), Rowv = T, Colv = T, col = col, 
          breaks=bk,symm=F, symbreaks=F, scale="none", cexRow = .2 ,
          trace="none", cexCol = .2,density.info = "none")
dev.off()


## BCR
source("./utils/quy_BCR_utils.R")
bcr.res <- read.csv("./BCR/filtered_contig_annotations.csv", stringsAsFactors = FALSE)

# separate by each sample
bcr.res.combined <- BCRprofle(bcr.res)


# compute heavy and light chain gene composition
# separate by each sample
bcr.res.combined.freq <- calClonalFreq(bcr.res.combined, cloneType="BCRgene")

bcr.barcode <- bcr.res.combined[,c("barcode","BCRgene")]
bcr.gene.matrix <- BCRgeneCompMatrix(bcr.barcode)
bcr.gene.matrix <- bcr.gene.matrix*100/nrow(bcr.barcode)
bk = unique(c(seq(0, 20, length=40)))

col = colorRampPalette(c("white","#162955"))(length(bk)-1)
pdf("heatmap_of_BCR_gene_composition.pdf")
heatmap.2(as.matrix(bcr.gene.matrix), Rowv = T, Colv = T, col = col, 
          breaks=bk,symm=F, symbreaks=F, scale="none", cexRow = .2 ,
          trace="none", cexCol = .2,density.info = "none")
dev.off()



