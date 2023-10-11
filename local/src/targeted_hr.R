library(ggplot2)
library(ggsignif)
d <- read.table('~/comet.txt', sep="\t" , header=T)
d2 <- read.table('~/comet_CRC0096.tsv', sep="\t" , header=T)
d3 <- read.table('~/comet_CRC0196.tsv', sep="\t" , header=T)


ggplot(data=d, aes(x=sample, y=circularity))+geom_violin()+geom_boxplot(width=0.1, outlier.shape=NA)+
geom_jitter()+theme_bw()+
geom_signif(comparisons = list(t1=c("322 plvx nt", "322 plvx sn38"), 
                               t2=c('322 prad51 nt', '322 prad51 sn38'),
                               t3=c("542 plvx nt", "542 plvx sn38"), 
                               t4=c('542 prad51 nt', '542 prad51 sn38')), test.args=list(alternative = c("greater")), test="wilcox.test" )

wilcox.test(d[d$sample== "322 plvx nt",'circularity'], d[d$sample== "322 prad51 nt",'circularity'], alternative ="less")

d <- rbind(d2, d3)
ggplot(data=d, aes(x=sample, y=circularity))+geom_violin()+geom_boxplot(width=0.1, outlier.shape=NA)+
  geom_jitter()+theme_bw()+
  geom_signif(comparisons = list(t1=c("96 plvx nt", "96 plvx sn38"), 
                                 t2=c('96 prad51 nt', '96 prad51 sn38'),
                                 t3=c("196 plvx nt", "196 plvx sn38"), 
                                 t4=c('196 rad51 nt', '196 rad51 sn38')), test.args=list(alternative = c("greater")), test="wilcox.test" )


ggplot(data=d, aes(x=sample, y=circularity))+geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = list(t1=c("96 plvx nt", "96 plvx sn38"), 
                                 t2=c('96 prad51 nt', '96 prad51 sn38'),
                                 t3=c("196 plvx nt", "196 plvx sn38"), 
                                 t4=c('196 rad51 nt', '196 rad51 sn38')), test.args=list(alternative = c("greater")), test="wilcox.test" )


ddef <- rbind(d, d2, d3)
ddef <- ddef[!grepl('542', ddef$sample),]
# split sample to have three columns
# iterate on models to get single plots for one model + do stats
# assemble with ggpubr

# poi statistica 4 vs 4 delta mediane/delta medie