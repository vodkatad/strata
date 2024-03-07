binary_f <- snakemake@input[['mutmat']]
wf_f <- snakemake@input[['wf']]
wil_f <- snakemake@output[['wilcox_genes']]
wanted <- snakemake@params[['wanted']]

mut <- read.table(binary_f, sep="\t")
wf <- read.table(wf_f, sep="\t", stringsAsFactors = FALSE, header=FALSE)
colnames(wf) <- c('smodel', 'perc')

data_merged <- merge(mut, wf, by.x="row.names", by.y='smodel')
if (snakemake@wildcards[['week']] == "w3") {
  stopifnot(nrow(data_merged)==wanted)
}

save.image('pi.Rdata')

fisher <- function(gene, mydata) {
  mydata$class <- ifelse(mydata$perc < -50, 'OR', ifelse(mydata$perc > 35, 'PD', 'SD'))
  ct <- table(mydata$class, mydata[,gene] != 0)
  if (ncol(ct) == 2) {
    ft <- fisher.test(ct)
    return(ft$p.value)
  } else {
    return(NA)
  }
}

fish <- as.data.frame(sapply(colnames(mut), fisher, data_merged))
colnames(fish) <- c('pvalue')
fish <- fish[!is.na(fish$pvalue),, drop=FALSE]
fish$padj <- p.adjust(fish$pvalue, method="BH")
fish <- fish[order(fish$pvalue),]

write.table(fish, file=wil_f, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
