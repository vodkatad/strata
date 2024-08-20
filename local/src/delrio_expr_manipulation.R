# BiocManager::install("hgu133plus2.db") ### array used by delrio
library(hgu133plus2.db)
library(AnnotationDbi)

delrio1 <- read.table(snakemake@input[['one']], sep='\t', quote="", header=TRUE, stringsAsFactors = FALSE)
delrio2 <- read.table(snakemake@input[['two']], sep='\t', quote="", header=TRUE, stringsAsFactors = FALSE)

mm <- merge(delrio1, delrio2, by=1)

mm$gene <- mapIds(hgu133plus2.db, keys=mm$ID_REF, column="SYMBOL", keytype="PROBEID")
mm <- mm[,c(1,147,2:146)]

write.table(mm, snakemake@output[['symbol']], sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)

mm$ID_REF <- NULL
mm <- mm[!is.na(mm$gene),]

### mean of same genes
genes <- unique(mm$gene)
new <- data.frame(matrix(nrow=length(genes), ncol=ncol(mm)-1))
colnames(new) <- colnames(mm)[2:ncol(mm)]
rownames(new) <- genes

for ( g in seq_along(genes) ) {
  tmp <- mm[mm$gene==genes[g],]
  mean <- data.frame(colMeans(tmp[,2:ncol(tmp)]))
  colnames(mean) <- genes[g] 
  mean <- as.data.frame(t(mean))
  new[genes[g],] <- mean
}

write.table(new, snakemake@output[['mean']], sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)


