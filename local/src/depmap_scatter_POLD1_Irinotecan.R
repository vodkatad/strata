library(tidyverse)
library(ggplot2)
library(openxlsx)

data_f <- snakemake@input[["POLD_Irino"]]
plot_cor_png <- snakemake@output[["png"]]
plot_cor_eps <- snakemake@output[["eps"]]
excel <- snakemake@output[["excel_cor"]]
log_f <- snakemake@log[['log']]

#data_f <- '/mnt/trcanmed/snaketree/prj/strata/local/share/data/POLD1_logTPMpc1_23Q2_Irinotecan_IC50_GDSC2_010823.csv'
data <- read.table(data_f, sep = ",", header = TRUE, stringsAsFactors = FALSE)

#load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_5.Rdata')
size <- 8

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

# Da Marti e https://www.christophenicault.com/post/understand_size_dimension_ggplot2/
#showtext_opts(dpi = 300) 
# since we are not changing fonts in the end cause myriad end up not being text object I'm not sure it's needed
#showtext_auto(enable = TRUE)

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
death_conversion_dpi96 = 96/72

textSize <- size * death_conversion_dpi96
largerSize <- (size) * death_conversion_dpi96

unmute_theme <- theme(
  text = element_text(size = textSize, family='Arial'),
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black", size=0.508/0.564), # origin of this ratio is honestly not known, empirical
  axis.ticks = element_line(color = "black", size=0.508/0.564),
  axis.ticks.length= unit(1.905*death_conversion_dpi96, "mm"),
  panel.background = element_blank()
)
#axis.ticks.length= unit(1.905*death_conversion_dpi96, "mm"),


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



colnames(data)[2] <- 'IC50_Irinotecan'
colnames(data)[3] <- 'POLD1'

sink(log_f, append=TRUE)
cor.test(data$IC50_Irinotecan, data$POLD1)
sink()

y_breaks <- guess_ticks(data$IC50_Irinotecan, fixed_min=-3, fixed_max=4)
x_breaks <- guess_ticks(data$POLD1, fixed_min=2, fixed_max=8)

p <- ggplot(data, aes(x=POLD1, y=IC50_Irinotecan)) +  geom_point(size=1) + geom_smooth(method='lm', size=1)+ # ratio is ... 1 becomes 0.939
  unmute_theme+theme(legend.position="none") + xlab('POLD1 log2(TPM+1)')+ylab('Irinotecan IC50')+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))

ggsave(p, file=plot_cor_png, width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")

p <- ggplot(data, aes(x=POLD1, y=IC50_Irinotecan)) +  geom_point(size=1) + geom_smooth(method='lm', size=1)+ # ratio is ... 1 becomes 0.939
  unmute_theme+theme(legend.position="none")+theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))

ggsave(p, file=plot_cor_eps, width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")

write.xlsx(data,file = excel)

