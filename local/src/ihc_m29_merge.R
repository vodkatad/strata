library(stringr)

rad51_basal_f <- snakemake@input[['RAD51_NT']]
rad51_irino_f <- snakemake@input[['RAD51_irino']]
h2ax_basal_f <- snakemake@input[['H2AX_NT']]
h2ax_irino_f <- snakemake@input[['H2AX_irino']]
pold_f <- snakemake@input[['POLD']]
hrd_f <- snakemake@input[['HRD']]
wf_f <- snakemake@input[['w3']]

rad51_basal <- read.table(rad51_basal_f, sep="\t")
rad51_irino <- read.table(rad51_irino_f, sep="\t")

pold1_basal <- read.table(pold_f, sep="\t")

w3 <- read.table(wf_f, sep="\t", row.names=1)
colnames(w3) <- 'dvw3'

H2AX_basal <- read.table(h2ax_basal_f, sep="\t")
H2AX_irino <- read.table(h2ax_irino_f, sep="\t")

hrd <- read.table(hrd_f, sep="\t", header=TRUE)

m <- merge(rad51_basal, pold1_basal, by="row.names")
rownames(m) <- m$Row.names
m$Row.names <- NULL
m2 <- merge(m, w3, by="row.names")
rownames(m2) <- m2$Row.names
m2$Row.names <- NULL

m2$recist <- ifelse(m2$dvw3 < -50, 'PR', ifelse(m2$dvw3 > 35, 'PD', 'SD'))
m2$trimmeddw3 <- ifelse(m2$dvw3 > 150, 150, m2$dvw3)

m3 <- merge(m2, rad51_irino, by='row.names')
rownames(m3) <- m3$Row.names
m3$Row.names <- NULL

m4 <- merge(m3, H2AX_basal, by="row.names")
rownames(m4) <- m4$Row.names
m4$Row.names <- NULL
m5 <- merge(m4, H2AX_irino, by="row.names")
rownames(m5) <- m5$Row.names
m5$Row.names <- NULL

m5$H2AX_induction <- m5$H2AX_irino / (m5$H2AX_NT+0.1)

#if (!all(hrd$CASI == as.numeric(unlist(str_extract_all(hrd$smodel, '\\d+'))))) {
#  stop('Issue in HRD scores models ids!')
#}

hrd$CASI <- NULL

m6 <- merge(m5, hrd, by.x="row.names", by.y="smodel")
rownames(m6) <- m6$Row.names
m6$Row.names <- NULL

write.table(m6, file=snakemake@output[['out']], sep="\t", quote=F, row.names=TRUE)
