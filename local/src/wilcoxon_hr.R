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

wilcoxon_perc <- function(gene, mydata) {
  percwt <- mydata[mydata[,gene] == 0, 'perc']
  percmut <- mydata[mydata[,gene] != 0, 'perc']
  if (length(percwt) != 0 && length(percmut) != 0) {
    wt <- wilcox.test(percwt, percmut)
    return(c(wt$p.value, length(percwt), length(percmut)))
  } else {
    return(c(NA, NA, NA))
  }
}

wil <- as.data.frame(t(sapply(colnames(mut), wilcoxon_perc, data_merged)))
colnames(wil) <- c('pvalue', 'nwt', 'nmut')
wil$padj <- p.adjust(wil$pvalue, method="BH")
wil <- wil[order(wil$pvalue),]

write.table(wil, file=wil_f, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
