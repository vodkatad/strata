data_f <- snakemake@input[["n"]]
stat_f <- snakemake@output[["stat"]]

save.image('pippo.Rdata')
### chisq 
set.seed(42)
d <- read.table(data_f, sep='\t', header=T, row.names=1)

sink(stat_f)
chisq.test(d, simulate.p.value = TRUE)

#chisq.test(d, simulate.p.value = FALSE)

sink()