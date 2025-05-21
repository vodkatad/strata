casi_f <- snakemake@input[["wanted"]]
muts_f <- snakemake@input[['mut']]
annot_f <- snakemake@input[['annot']]
models_f <- snakemake@input[['all_models']]
vaffile <- snakemake@output[['vaf']]
protfile <- snakemake@output[['prot']]
thr <- as.numeric(snakemake@params[['thr']])

mdata <- read.table(muts_f, sep="\t" , header=TRUE)
anndata <- read.table(annot_f, sep="\t" , header=TRUE, stringsAsFactors = FALSE)
all_models <- read.table(models_f, sep="\t" , header=FALSE, stringsAsFactors = FALSE)
wanted <- read.table(casi_f, sep="\t", header=FALSE, stringsAsFactors = FALSE)

save.image('pippo.Rdata') 
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

wanted_crc <- intersect(wanted[,1], colnames(mdata))

# now we need to collapse to genes
genes <- unique(anndata$Gene.refGene)
#binary <- data.frame(matrix(rep(0, length(genes)*length(newcols)), nrow=length(genes)), row.names=unique(anndata$Gene.refGene))
#colnames(binary) <- newcols

search_gene_muts_vaf <- function(gene, muts, annot) {
  genemuts <- annot[annot$Gene.refGene == gene,]
  genemuts <- genemuts[genemuts$ExonicFunc.refGene != 'synonymous SNV',]
  if (nrow(genemuts) > 0) {
    genedata <- muts[rownames(muts) %in% rownames(genemuts),, drop=FALSE]
    vafs <- apply(genedata, 2, function(x){paste0(x[x>thr], collapse=",")})
    return(vafs)
  } else {
    print(gene)
    return(rep('', ncol(muts)))
  }
}

vaf <- sapply(genes, search_gene_muts_vaf, mdata, anndata)
vaf <- vaf[wanted_crc,]

search_gene_muts_prot <- function(gene, muts, annot) {
  genemuts <- annot[annot$Gene.refGene == gene,]
  genemuts <- genemuts[genemuts$ExonicFunc.refGene != 'synonymous SNV',]
  genemuts$prot <- sapply(strsplit(genemuts$AAChange.refGene, ':'), function(x){x[[length(x)]][1]})
  if (nrow(genemuts) > 0) {
    genedata <- muts[rownames(muts) %in% rownames(genemuts),, drop=FALSE]
    # genemuts$prot was not ordered like genedata...need to get the right mut in the right way.
    genemuts <- genemuts[rownames(genemuts) %in% rownames(genedata),,drop=FALSE]
    genemuts <- genemuts[rownames(genedata),, drop=FALSE]
    res <- data.frame(row.names=rownames(genedata))
    for (g_i in seq(1, nrow(genemuts))) {
      for (m_i in seq(1, ncol(genedata))) {
        if (genedata[g_i, m_i] > thr) {
          res[g_i, m_i] <- genemuts[g_i, 'prot']
        } else {
          res[g_i, m_i] <- ''
        }
      } 
    }
    #for (i in seq(1, ncol(genedata))) {
    #  res[, colnames(genedata)[i]] <- ifelse(genedata[,i] > thr, genemuts$prot, '')  
    #}
    collapsed <- apply(res, 2, function(x){paste0(x[x!=''], collapse=",")})
    names(collapsed) <- colnames(genedata)
    return(collapsed)
  } else {
    return(rep('', ncol(muts)))
  }
}

prot <- sapply(genes, search_gene_muts_prot, mdata, anndata)
prot <- prot[wanted_crc,]
#prot[prot==""] <- 'WT'

write.table(vaf, file=vaffile, sep="\t", quote=FALSE)
write.table(prot, file=protfile, sep="\t", quote=FALSE)
