library(ggplot2)
d2 <- read.table('/mnt/cold1/snaketree/prj/snakegatk_real/dataset/strata_shallow/all_aligned_dedup.tsv', '\t', header=F)
d1 <- read.table('/mnt/trcanmed/snaketree/prj/snakegatk/dataset/biobanca_shallowseq_pdx/all_aligned_dedup.tsv', '\t', header=F)
colnames(d1) <- c('id', 'aligned_dedup_reads', 'mean_coverage')
colnames(d2) <- c('id', 'aligned_dedup_reads', 'mean_coverage')

ggplot(data=d1, aes(x=mean_coverage))+geom_histogram(bins=20)+theme_bw()+theme(text=element_text(size=20))+coord_cartesian(xlim=c(0.2, 0.8))
ggplot(data=d2, aes(x=mean_coverage))+geom_histogram(bins=20)+theme_bw()+theme(text=element_text(size=20))+coord_cartesian(xlim=c(0.2, 0.8))

