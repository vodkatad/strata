xeno_df <- read.table(snakemake@input[['cn']], quote="", sep="\t", header=TRUE, row.names = 1)


xeno_df$start <- NULL
xeno_df$end <- NULL
xeno_df$chromosome <- NULL
#pc <- min(xeno_df_bck[xeno_df_bck>0], pdo_df_bck[pdo_df_bck>0])
#pc <- 0.01
pc <- 0.01
xeno_df <- log2(xeno_df+pc)

write.table(xeno_df, file=gzfile(snakemake@output[['cn']]), quote=FALSE, sep="\t")

