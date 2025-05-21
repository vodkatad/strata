muts_f <- snakemake@input[['mut']]
annot_f <- snakemake@input[['annot']]
models_f <- snakemake@input[['all_models']]
wf_f <- snakemake@input[['wf_models']]
outfile <- snakemake@output[['mat']]
rdata <- snakemake@output[['data']]
thr <- as.numeric(snakemake@params[['thr']])
log_f <- snakemake@log[['log']]

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
               if (res[i,j]==0) {
                  res[i,j] = 0.5 # we want to keep muts with 0 deleterious prediction anyway
               }
            }
         }
      }
      return(apply(res, 2, max))
   } else {
     print(gene)
     return(rep(0, ncol(muts)))
   }
}

ndeldf <- sapply(genes, search_gene_muts_nDel, mdata, anndata)

toadd <- setdiff(rownames(ndeldf), all_models[,1])
if (length(toadd) > 0) {
  stop('TODO implement add samples without any mut!')
}

write.table(ndeldf, outfile, sep= "\t", quote=FALSE)
save.image(rdata)

tot <- length(logdata[logdata!=0])
monoallelic <- length(logdata[logdata!=0 & logdata<=0.99])
biallelic <- length(logdata[logdata>0.99])

# access annotations of muts mutated in at least one sample
tot_mut <- rownames(logdata)[apply(logdata, 1, function(x){any(x!=0)})]
log_ann <- anndata[rownames(anndata) %in% tot_mut,]
ndel <- as.data.frame(apply(log_ann[, wanted_predictions], 1, count_letter))
colnames(ndel) <- 'n'
ndel_n <- sum(ndel$n >= 3)
sink(log_f)
print('tot muts counting double same mut in two samples')
print(tot)
print('monoallelic mut')
print(monoallelic)
print((monoallelic/tot)*100)
print('biallelic mut')
print(biallelic)
print((biallelic/tot)*100)
print('tot muts "normal"')
print(length(tot_mut))
print('tot muts with >= 3/5 d')
print(ndel_n)
print((ndel_n/length(tot_mut))*100)
sink()