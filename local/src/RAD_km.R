library(survival)
library(survminer)

data_f <- snakemake@input[["surviv"]]
stat_f <- snakemake@output[["stat"]]
##
d <- read.table('/mnt/trcanmed/snaketree/prj/strata/local/share/data/rad51_km.txt', sep="\t",header=TRUE, stringsAsFactors = FALSE)
d$hospital <- sapply(strsplit(d$ID, " "), '[[', 1)
table(d$hospital, d$quartile)

folfiri_nona <- d
folfiri_nona$status <- ifelse(d$X == 1, 1,0)

folfiri_nona$quartilimerged <- ifelse(folfiri_nona$quartile %in% c('0','1+','2+'), 'low','high')
sfit <- survfit(Surv(PFS, status) ~ quartilimerged, data=folfiri_nona)
#ggsurvplot(sfit, data = folfiri_nona, pval = TRUE) # pval.method=Log-rank
#surv_pvalue(sfit)

#RESISTANT: ROSSO (R 165 G 0 B 25)
#SENSITIVE: AZZURRO (R 30 G 85 B 130)

#th <- theme_minimal() + scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1),expand = c(0, 0))
pl <- ggsurvplot(sfit, data = folfiri_nona, pval = TRUE, palette=c('#A50019','#1E5582'), risk.table='absolute')

pl$plot + scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits=c(0, 1),expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
