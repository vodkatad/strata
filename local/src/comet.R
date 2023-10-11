library(ggplot2)
library(ggsignif)
library(ggpubr)


data_f <- snakemake@input[["comet"]]
plot_f <- snakemake@output[["boxplot"]]
log_f <- snakemake@log[['log']]


save.image('pippo.Rdata')
d <- read.table(data_f, sep="\t" , header=TRUE, stringsAsFactors = FALSE)

# p  <- ggplot(data=d, aes(x=sample, y=circularity))+geom_boxplot()+
#   theme_bw()+
#   geom_signif(comparisons = list(t1=c("96 plvx nt", "96 plvx sn38"), 
#                                  t2=c('96 prad51 nt', '96 prad51 sn38'),
#                                  t3=c("196 plvx nt", "196 plvx sn38"), 
#                                  t4=c('196 rad51 nt', '196 rad51 sn38'),
#                                  t5=c('729 plvx nt', '729 plvx sn38'),
#                                  t6=c('729 prad51 nt', '729 prad51 sn38')), test.args=list(alternative = c("greater")), test="wilcox.test" )

spl <- strsplit(d$sample, ' ')
d$model <- as.character(sapply(spl, '[[', 1 ))
d$vector <- sapply(spl, '[[', 2 )
d$treatment <- sapply(spl, '[[', 3 )
pad <- 4 - nchar(d$model)
for (i in seq(1, length(pad))) {
  d$model[i] <- paste0('CRC',  paste0(rep('0', pad[i]), collapse=''), d$model[i])
}

d$vector <- gsub('plvx', 'pLVX', d$vector)
d$vector <- gsub('prad51', 'pRAD51', d$vector)
d$vector <- gsub('rad51', 'pRAD51', d$vector)
d$treatment <- gsub('nt', 'Untreated', d$treatment)
d$treatment <- gsub('sn38', 'Treated', d$treatment)

d <- d[d$model != "CRC0542",]
oneplot <- function(mod, data) {
  pd <- d[d$model==mod,]
  #relevel factor SN38 Untreated
  pd$treatment <- factor(pd$treatment, levels=c('Untreated', 'Treated'))
  ggplot(data=pd, aes(x=treatment, y=circularity))+geom_boxplot()+facet_wrap(~vector,strip.position = 'right')+theme_bw()+xlab('')
  # statistics
}


plots <- lapply(unique(d$model), oneplot, d)

p <- ggarrange(plotlist=plots, nrow=length(unique(d$model)), labels=unique(d$model), label.x=0.42, label.y=1.2)
ggsave(file=plot_f, p, height=89, width=56, units="mm")

#colori da figura Marco A
# 542 togliere
# split sample to have three columns
# iterate on models to get single plots for one model + do stats
# assemble with ggpubr

# poi statistica 4 vs 4 delta mediane/delta medie