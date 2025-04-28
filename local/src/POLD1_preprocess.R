library(ggplot2)
library(stringr)

fra_f <- snakemake@input[["fra"]]
rad51_f <- snakemake@input[["orig_pfs"]]
plot <- snakemake@output[['outplot']]
outtsv <- snakemake@output[['outtsv']]
outtsvboth <- snakemake@output[['outtsvboth']]
log_f <- snakemake@log[['log']]
outchisq <- snakemake@output[['outchisq']]

save.image('p.Rdata')

pold1 <- read.table(fra_f, sep="\t", stringsAsFactors = FALSE, header=TRUE)
rad51 <- read.table(rad51_f, sep="\t", stringsAsFactors = FALSE, header=TRUE)

pold1$response <- factor(x=pold1$response, levels=rev(c('CR', 'PR', 'SD', 'PD')))
# ggplot(data=pold1, aes(x=reorder(pid_new, POLD1), y=POLD1, fill=response))+geom_col()+
#   theme_bw(base_size = 15)+ theme( axis.text.x=element_blank())+scale_fill_brewer(palette = "Spectral")

p <- ggplot(data=pold1, aes(x=reorder(pid_new, POLD1), y=POLD1, fill=POLD1_class))+geom_col()+
  theme_bw(base_size = 15)+ theme( axis.text.x=element_blank())+scale_fill_brewer(palette="YlOrRd")
ggsave(plot, plot=p)

id_numeric <- unlist(str_extract_all(pold1$pid_new, '\\d+'))
id_hospital <- toupper(unlist(str_extract_all(pold1$pid_new, '[:alpha:]+')))
pold1$ID <- paste0(id_hospital, ' ', id_numeric)

rad51[rad51$ID=='HMAR 56', 'ID'] = 'HMAR 42'
rad51[rad51$ID=='HMAR 55', 'ID'] = 'HMAR 46'

n1 <- length(intersect(pold1$ID, rad51$ID))
n2 <- nrow(rad51)
n3 <- nrow(pold1)

stopifnot('problems pids 1-2'= n1==n2)
stopifnot('problems pids 2-3'= n2==n3)

m <- merge(pold1, rad51, by='ID')

sink(log_f)
print('We want 82:')
print(nrow(m))
print(all(m$PFS.x == m$PFS.y))
sink()
# m$p <- as.numeric(unlist(str_extract_all(m$POLD1_class, '\\d+')))
# m$r <- as.numeric(unlist(str_extract_all(m$quartile, '\\d+')))
# table(m$p, m$r)
mboth <- m
m$quartile <- m$POLD1_class

stopifnot('PFS'= all(m$PFS.x == m$PFS.y))
m$PFS <- m$PFS.x
m$PFS.x <- NULL
m$PFS.y <- NULL
write.table(m, file=outtsv, sep="\t", quote=FALSE, row.names=FALSE)

# setwd('/mnt/trcanmed/snaketree/prj/strata/dataset/figures')
# load('p.Rdata')
forchisq <- as.matrix(table(m$POLD1_class, m$response))
write.table(forchisq, file=outchisq, sep="\t", quote=FALSE, row.names=TRUE)

colnames(mboth)[colnames(mboth)=='quartile'] <- 'RAD51_class'
mboth$POLD1 <- NULL
#mboth$response <- NULL
mboth$PFS <- mboth$PFS.x
mboth$PFS.x <- NULL
mboth$PFS.y <- NULL
write.table(mboth, file=outtsvboth, sep="\t", quote=FALSE, row.names=FALSE)

q('no')
ggplot(data=m, aes(x=response, y=POLD1_class))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size = 13)

