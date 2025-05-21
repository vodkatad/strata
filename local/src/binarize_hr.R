muts_f <- snakemake@input[['mut']]
annot_f <- snakemake@input[['annot']]
models_f <- snakemake@input[['all_models']]
outfile <- snakemake@output[['mat']]
thr <- as.numeric(snakemake@params[['thr']])

mdata <- read.table(muts_f, sep="\t" , header=TRUE)
anndata <- read.table(annot_f, sep="\t" , header=TRUE, stringsAsFactors = FALSE)
all_models <- read.table(models_f, sep="\t" , header=FALSE, stringsAsFactors = FALSE)

# reduce longgen to smodel
newcols <- unique(substr(colnames(mdata)[-1], 0,7))
if (ncol(mdata)-1 != length(newcols)) {
    stop('Error, CRC models are repeated')
}

#mut ids conversions and to rownames
# CRC0196 and ATM WTF
anndata[anndata$Ref == '-', 'Ref'] <- ''
anndata[anndata$Alt == '-', 'Alt'] <- ''
rownames(anndata) <- paste0('chr', anndata$Chr, ":", anndata$Start, ":", anndata$Ref, ":", anndata$Alt)
rownames(mdata) <- mdata$ID
mdata$ID <- NULL
colnames(mdata) <- newcols

mdata <- ifelse(mdata>thr, 1, 0)

# now we need to collapse to genes
genes <- unique(anndata$Gene.refGene)
#binary <- data.frame(matrix(rep(0, length(genes)*length(newcols)), nrow=length(genes)), row.names=unique(anndata$Gene.refGene))
#colnames(binary) <- newcols

search_gene_muts <- function(gene, muts, annot) {
   genemuts <- annot[annot$Gene.refGene == gene,]
   genemuts <- genemuts[genemuts$ExonicFunc.refGene != 'synonymous SNV',]
   if (nrow(genemuts) > 0) {
      genedata <- muts[rownames(muts) %in% rownames(genemuts),, drop=FALSE]
      return(colSums(genedata))
   } else {
     print(gene)
     return(rep(0, ncol(muts)))
   }
}

binary <- sapply(genes, search_gene_muts, mdata, anndata)
binary <- ifelse(binary>0, 1, 0)

toadd <- setdiff(rownames(binary), all_models[,1])
if (length(toadd) > 0) {
  stop('TODO implement add samples without any mut!')
}

write.table(binary, outfile, sep= "\t", quote=FALSE)
#save.image('pippo.Rdata')