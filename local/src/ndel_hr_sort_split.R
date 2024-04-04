muts_f <- snakemake@input[['mut']]
annot_f <- snakemake@input[['annot']]
models_f <- snakemake@input[['all_models']]
wf_f <- snakemake@input[['wf_models']]
sort_f <- snakemake@input[['sort']]
outfile <- snakemake@output[['mat']]
thr <- as.numeric(snakemake@params[['thr']])

sortdf <- read.table(sort_f, sep="\t", header=TRUE, stringsAsFactors = FALSE)

mdata <- read.table(muts_f, sep="\t" , header=TRUE)
anndata <- read.table(annot_f, sep="\t" , header=TRUE, stringsAsFactors = FALSE)
all_models <- read.table(models_f, sep="\t" , header=FALSE, stringsAsFactors = FALSE)

# reduce longgen to smodel
newcols <- unique(substr(colnames(mdata)[-1], 0,7))
if (ncol(mdata)-1 != length(newcols)) {
    stop('Error, CRC models are repeated')
}

#mut ids conversions and to rownames
rownames(anndata) <- paste0('chr', anndata$Chr, ":", anndata$Start, ":", anndata$Ref, ":", anndata$Alt)
rownames(mdata) <- mdata$ID
mdata$ID <- NULL
colnames(mdata) <- newcols

##### reduce to wf for results numbers
wf <- read.table(wf_f, sep="\t", stringsAsFactors = FALSE)
colnames(wf) <- c('smodel', 'perc')

keep <- intersect(colnames(mdata), wf$smodel)
logdata <- mdata[ , colnames(mdata) %in% keep,]
#####

mdata <- ifelse(mdata>thr, 1, 0)

# now we need to collapse to genes
genes <- unique(anndata$Gene.refGene)
#binary <- data.frame(matrix(rep(0, length(genes)*length(newcols)), nrow=length(genes)), row.names=unique(anndata$Gene.refGene))
#colnames(binary) <- newcols
wanted_predictions <- c('SIFT_pred', 'fathmm.MKL_coding_pred', 'LRT_pred', 'FATHMM_pred', 'PROVEAN_pred')

count_letter <- function(x, letter='D') {
   length(x[x==letter])
}

# if one sample has two muts on the same gene we keep the worst one (more Deleterious predicion)
# -> no here we keep both with a , in between
search_gene_muts_nDel <- function(gene, muts, annot) {
   genemuts <- annot[annot$Gene.refGene == gene,]
   genemuts <- genemuts[genemuts$ExonicFunc.refGene != 'synonymous SNV',]
   ndel <- as.data.frame(apply(genemuts[, wanted_predictions], 1, count_letter))
   colnames(ndel) <- 'n'
   if (nrow(genemuts) > 0) {
      genedata <- muts[rownames(muts) %in% rownames(genemuts),, drop=FALSE]
      res <- genedata
      for (i in seq(1, nrow(genedata))) {
         for (j in seq(1, ncol(genedata))) {
            if (genedata[i,j] >= 1) {
               res[i,j] <- ndel[rownames(ndel)==rownames(genedata)[i], 'n']
            } else {
               res[i,j] <- ''
            }
         }
      }
      mypaste <- function(x) { x <- x[x!=""]; paste(x, collapse=",")}
      return(apply(res, 2, mypaste))
   } else {
     print(gene)
     return(rep('', ncol(muts)))
   }
}

ndeldf <- sapply(genes, search_gene_muts_nDel, mdata, anndata)

toadd <- setdiff(rownames(ndeldf), all_models[,1])
if (length(toadd) > 0) {
  stop('TODO implement add samples without any mut!')
}

ndeldf <- ndeldf[rownames(sortdf), colnames(sortdf)]

write.table(ndeldf, outfile, sep= "\t", quote=FALSE)
save.image('pippo.Rdata')
