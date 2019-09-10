af_file <- snakemake@input[["af"]]
out_file <- snakemake@output[["minaf"]]


d <- read.table(gzfile(af_file), sep="\t", header=F)
all <- unique(d$V1)
mins <- sapply(all, function(x) {min(d[d$V1==x,]$V2)})
df <- data.frame(id=all, min=mins)
write.table(df,gzfile(out_file), col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")

