loweraf <- snakemake@params[["loweraf"]]
data <- snakemake@input[["afmatrix"]]
afcolumn <- snakemake@params[["afcolumn"]]
fit <- snakemake@output[["fit"]]
histo <- snakemake@output[["hist"]]
debug <- snakemake@params[["debug"]]

if (debug == "yes") {
  save.image(file=paste0(fit,'.debug','.RData'))
}

data <- read.table(gzfile(data), sep="\t", header=TRUE)
af <- data[,afcolumn, drop=FALSE]

exsubcl <- af[af<loweraf & af != 0]
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
invf <- 1/exsubcl

pdf(fit)
plot(invf, excum)
graphics.off()
pdf(histo)
hist(exsubcl, breaks=50)
graphics.off()

if (debug == "yes") {
  save.image(file=paste0(fit,'.debug','.RData'))
}

