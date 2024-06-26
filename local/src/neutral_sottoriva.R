loweraf <- snakemake@wildcards[["loweraf"]]
higheraf <- snakemake@wildcards[["higheraf"]]
data <- snakemake@input[["afmatrix"]]
afcolumn <- snakemake@params[["afcolumn"]]
fit <- snakemake@output[["fit"]]
histo <- snakemake@output[["hist"]]
debug <- snakemake@params[["debug"]]
r2 <- snakemake@output[["r2"]]

if (debug == "yes") {
  save.image(file=paste0(fit,'.debug','.RData'))
}

data <- read.table(gzfile(data), sep="\t", header=TRUE)
af <- data[,afcolumn, drop=FALSE]

exsubcl <- af[af<higheraf & af>loweraf]
excum <- sapply(1:length(exsubcl),function(i)sum(exsubcl[i]<=exsubcl[1:length(exsubcl)]))
invf <- 1/exsubcl
maxi <- length(invf)
labels <-  c(1, floor(maxi/5), floor(maxi/2), floor(maxi/0.5), maxi)
model <- lm(excum~invf)
sfit <- summary(model)
dr2 <- data.frame(r=sfit$r.squared)
write.table(dr2, file=r2, sep="\t", quote=F)
pdf(fit)
plot(invf, excum, xaxt='n', main=paste0("R2=", round(sfit$r.squared, digits=3)))
oi <- invf[order(invf)]
oex <- exsubcl[order(-exsubcl)]
axis(1, at=oi[labels],labels=paste0("1/",oex[labels]), las=2)
abline(model)
graphics.off()
pdf(histo)
hist(exsubcl, breaks=50)
graphics.off()

if (debug == "yes") {
  save.image(file=paste0(fit,'.debug','.RData'))
}

