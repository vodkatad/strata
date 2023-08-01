d <- read.table('~/comet.txt', sep="\t" , header=T)

ggplot(data=d, aes(x=sample, y=circularity))+geom_violin()+geom_boxplot(width=0.1, outlier.shape=NA)+
geom_jitter()+theme_bw()+
geom_signif(comparisons = list(t1=c("322 plvx nt", "322 plvx sn38"), 
                               t2=c('322 prad51 nt', '322 prad51 sn38'),
                               t3=c("542 plvx nt", "542 plvx sn38"), 
                               t4=c('542 prad51 nt', '542 prad51 sn38')), test.args=list(alternative = c("greater")), test="wilcox.test" )

wilcox.test(d[d$sample== "322 plvx nt",'circularity'], d[d$sample== "322 prad51 nt",'circularity'], alternative ="less")
