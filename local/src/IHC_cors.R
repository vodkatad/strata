library(ggplot2)
library(RColorBrewer)
setwd('/mnt/trcanmed/snaketree/prj/strata/dataset/figures')

rad51_basal <- read.table('avgall_RAD51_NT.tsv', sep="\t")
rad51_irino <- read.table('avgall_RAD51_irino.tsv', sep="\t")

pold1_basal <- read.table('avgall_POLD1.tsv', sep="\t")

w3 <- read.table('w3_waterfall.tsv', sep="\t", row.names=1)
colnames(w3) <- 'dvw3'


m <- merge(rad51_basal, pold1_basal, by="row.names")
rownames(m) <- m$Row.names
m$Row.names <- NULL
m2 <- merge(m, w3, by="row.names")
rownames(m2) <- m2$Row.names
m2$Row.names <- NULL

m2$recist <- ifelse(m2$dvw3 < -50, 'PR', ifelse(m2$dvw3 > 35, 'PD', 'SD'))
m2$trimmeddw3 <- ifelse(m2$dvw3 > 150, 150, m2$dvw3)
# add Eugy list of treated with 5FU too
ggplot(data=m2, aes(x=POLD1, y=RAD51_NT, color=trimmeddw3))+geom_point()+theme_bw() +scale_colour_gradient2(midpoint = 35, low = "blue",
                                                                                                      mid = "white",
                                                                                                      high = "red")

ggplot(data=m2, aes(x=POLD1, y=RAD51_NT, color=recist))+geom_point()+theme_bw()
ggplot(data=m2, aes(x=RAD51_NT, y=dvw3))+geom_point()+theme_bw()
ggplot(data=m2, aes(x=POLD1, y=dvw3))+geom_point()+theme_bw()


luraghi <- read.table('/mnt/trcanmed/snaketree/prj/strata/local/share/data/luraghi', sep="\t")
w3$luraghi <- ifelse(rownames(w3) %in% luraghi$V1, 'folfiri', 'irino')

# non sembra un gran problema.
w3$model <- rownames(w3)
ggplot(data=w3, aes(order(-dvw3, model), y=dvw3, fill=luraghi))+geom_col()+
  theme_bw()+geom_hline(yintercept=-50)+geom_hline(yintercept=35)


## RAD51 vs response
m3 <- merge(m2, rad51_irino, by='row.names')
rownames(m3) <- m3$Row.names
m3$Row.names <- NULL
