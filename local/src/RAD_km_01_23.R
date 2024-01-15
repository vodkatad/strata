library(survival)
library(survminer)

data_f <- snakemake@input[["surviv"]]
plot <- snakemake@output[['plot']]
log_f <- snakemake@log[['log']]
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

##
d <- read.table(data_f, sep="\t",header=TRUE, stringsAsFactors = FALSE)
d$hospital <- sapply(strsplit(d$ID, " "), '[[', 1)
table(d$hospital, d$quartile)

folfiri_nona <- d
folfiri_nona$status <- ifelse(d$X == 1, 1,0)
#lab meeting: colonna AA se == 1 evento altrimenti perso al followup.
# 1 progression 2 toxicity 3 other

folfiri_nona$quartilimerged <- ifelse(folfiri_nona$quartile %in% c('0','1+'), '0, 1+','2+, 3+')
sfit <- survfit(Surv(PFS, status) ~ quartilimerged, data=folfiri_nona)
#ggsurvplot(sfit, data = folfiri_nona, pval = TRUE) # pval.method=Log-rank
#surv_pvalue(sfit)

#RESISTANT: ROSSO (R 165 G 0 B 25)
#SENSITIVE: AZZURRO (R 30 G 85 B 130)

#th <- theme_minimal() + scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1),expand = c(0, 0))
pl <- ggsurvplot(sfit, data = folfiri_nona, pval = FALSE, palette=c('#1E5582', '#A50019'), risk.table=FALSE)

p1 <- pl$plot + scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1),expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0), limits=c(0, 3000), breaks=c(0, 500, 1000, 1500, 2000))

p1 <- p1+unmute_theme
ggsave(file = plot, plot = p1, width=88*death_conversion_dpi96, height=64*death_conversion_dpi96, units="mm")

sink(log_f)
surv_pvalue(sfit)
sink()

q('no')

folfiri_nona$hospital <- as.factor(folfiri_nona$hospital)
sfit2 <- survfit(Surv(PFS, status) ~ quartilimerged+hospital, data=folfiri_nona)
table(folfiri_nona$hospital, folfiri_nona$X.1)
table(folfiri_nona$hospital, folfiri_nona$status)
plot(table(folfiri_nona$hospital, folfiri_nona$status))
ggplot(data=folfiri_nona, aes(x=hospital, y=PFS, fill=quartilimerged))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_minimal()


ggsurvplot(sfit2, data = folfiri_nona, pval = TRUE, risk.table=FALSE)

folfiri_nonin <- folfiri_nona[folfiri_nona$hospital != "NIG",]
sfit <- survfit(Surv(PFS, status) ~ quartilimerged, data=folfiri_nonin)
pl <- ggsurvplot(sfit, data = folfiri_nonin, pval = FALSE, palette=c('#1E5582', '#A50019'), risk.table=TRUE)

surv_pvalue(sfit)
