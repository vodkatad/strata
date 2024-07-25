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
w3$luraghi <- ifelse(rownames(w3) %in% luraghi$V1, 'irino', 'folfiri')

# non sembra un gran problema.
w3$model <- rownames(w3)
ggplot(data=w3, aes(order(-dvw3, model), y=dvw3, fill=luraghi))+geom_col()+
  theme_bw()+geom_hline(yintercept=-50)+geom_hline(yintercept=35)


## RAD51 vs response
m3 <- merge(m2, rad51_irino, by='row.names')
rownames(m3) <- m3$Row.names
m3$Row.names <- NULL

plot_lmsmooth_info <- function(xvar, yvar, md) {
  pear <- cor.test(md[,xvar], md[,yvar])
  title <- paste0("pearson=", signif(pear$estimate, 3), " P=", signif(pear$p.value, 3))
  print(ggplot(data=md, aes_string(x=xvar, y=yvar))+geom_point()+theme_bw()+geom_smooth(method="lm")+ggtitle(title))
}


plot_lmsmooth_info('RAD51_NT', 'dvw3', m3)
plot_lmsmooth_info('RAD51_irino', 'dvw3', m3)
m3$RAD51_induction <- m3$RAD51_irino / m3$RAD51_NT
plot_lmsmooth_info('RAD51_induction', 'dvw3', m3)

# H2ax

H2AX_basal <- read.table('avgall_H2AX_NT.tsv', sep="\t")
H2AX_irino <- read.table('avgall_H2AX_irino.tsv', sep="\t")


m4 <- merge(m3, H2AX_basal, by="row.names")
rownames(m4) <- m4$Row.names
m4$Row.names <- NULL
m5 <- merge(m4, H2AX_irino, by="row.names")
rownames(m5) <- m5$Row.names
m5$Row.names <- NULL

m5$H2AX_induction <- m5$H2AX_irino / (m5$H2AX_NT+0.1)
plot_lmsmooth_info('H2AX_NT', 'dvw3', m5)
plot_lmsmooth_info('H2AX_irino', 'dvw3', m5)
plot_lmsmooth_info('H2AX_induction', 'dvw3', m5)


plot_lmsmooth_info('H2AX_NT', 'RAD51_NT', m5)
plot_lmsmooth_info('H2AX_irino', 'RAD51_irino', m5)
plot_lmsmooth_info('H2AX_induction', 'RAD51_induction', m5)

plot_lmsmooth_info('H2AX_induction', 'H2AX_irino', m5)

m5$luraghi <- ifelse(rownames(m5) %in% luraghi$V1, 'irino', 'folfiri')
ggplot(data=m5, aes(x=H2AX_induction, y=RAD51_induction, color=luraghi))+geom_point()+theme_bw() 


ggplot(data=m5, aes(x=POLD1, y=RAD51_NT, color=trimmeddw3, shape=luraghi))+geom_point()+theme_bw() +scale_colour_gradient2(midpoint = 35, low = "blue",
                                                                                                            mid = "white",
                                                                                                            high = "red")




ggplot(data=m5, aes(x=POLD1, y=RAD51_irino, color=trimmeddw3, shape=luraghi))+geom_point()+theme_bw() +scale_colour_gradient2(midpoint = 35, low = "blue",
                                                                                                                           mid = "white",
                                                                                                                           high = "red")




