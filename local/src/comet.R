library(ggplot2)
library(ggsignif)
library(ggpubr)


textSize <- 8

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
unmute_theme <- theme(
  text = element_text(size = textSize, family='sans'),
  axis.title = element_text(size = textSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = textSize, hjust = 0.5),
  legend.title = element_text(size=textSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(color = "black"),
  panel.background = element_blank()
)


mute_theme <- theme_bw() +
  theme(
    text = element_blank()
  )


# function that given values to be plotted on an axis will return:
# vector of breaks, trying to guess which max will be the best one
# this will be used as scale_y_continuous(breaks=  and as ylim(min, max) to have the - also limits-c()
# last tick at the extremity of the axis.
# other parameter is n. of ticks
guess_ticks <- function(values, nticks=5, fixed_max=NULL, fixed_min=0) {
  vmax <- max(values)
  if (is.null(fixed_max)) { 
    round_max <- ceiling(vmax)
  } else {
    round_max <- fixed_max
  }
  my_breaks <- seq(fixed_min, round_max, length.out=nticks)
  return(my_breaks)
}


data_f <- snakemake@input[["comet"]]
plot_f <- snakemake@output[["boxplot"]]
log_f <- snakemake@log[['log']]

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
  ggplot(data=pd, aes(x=treatment, y=circularity))+geom_boxplot(outlier.size = 0.15, size=0.2)+
    facet_wrap(~vector,strip.position = 'right')+xlab('')+
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1),expand = c(0, 0))+
    unmute_theme#+
    #scale_color_manual(values=c('black', '#3d752fff'))+theme(legend.position="none")
}


plots <- lapply(unique(d$model), oneplot, d)

p <- ggarrange(plotlist=plots, nrow=length(unique(d$model)), 
               labels=unique(d$model), label.x=0.44, label.y=1.05, font.label = list(size = 8))
#ggsave(file=plot_f, p, height=80, width=78, units="mm", device = cairo_ps)
ggsave(file=plot_f, p, height=80, width=78, units="mm", device = cairo_ps)

#colori da figura Marco A
#fix first label
# [542 togliere
# [split sample to have three columns
# [iterate on models to get single plots for one model + do stats
# [assemble with ggpubr

# poi statistica 4 vs 4 delta mediane/delta medie


onestat <- function(mod, data) {
  pd <- data[data$model==mod,]
  #relevel factor SN38 Untreated
  plvx <- pd[pd$vector=="pLVX",]
  prad51 <- pd[pd$vector=="pRAD51",]
  wmock <- wilcox.test(x=plvx[plvx$treatment=='Untreated', 'circularity'], plvx[plvx$treatment=='Treated', 'circularity'], alternative="greater")
  wrad <- wilcox.test(x=prad51[prad51$treatment=='Untreated', 'circularity'], prad51[prad51$treatment=='Treated', 'circularity'], alternative="greater")
  c(wmock$statistic, wmock$p.value, wrad$statistic, wrad$p.value)
}


stats <- sapply(unique(d$model), onestat, d)

sink(log_f)
print(stats)
sink()


deltaavg <- function(mod, data) {
  pd <- data[data$model==mod,]
  #relevel factor SN38 Untreated
  plvx <- pd[pd$vector=="pLVX",]
  prad51 <- pd[pd$vector=="pRAD51",]
  dmock <- mean(plvx[plvx$treatment=='Untreated', 'circularity']) - mean(plvx[plvx$treatment=='Treated', 'circularity'])
  drad <- mean(x=prad51[prad51$treatment=='Untreated', 'circularity']) - mean(prad51[prad51$treatment=='Treated', 'circularity'])
  c(dmock,drad)
}


deltas <- sapply(unique(d$model), deltaavg, d)

sink(log_f, append=TRUE)
wilcox.test(deltas[1,], deltas[2,], alternative='greater')
sink()
