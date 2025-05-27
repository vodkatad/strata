library(ggplot2)
library(RColorBrewer)

size <- 8

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
death_conversion_dpi96 = 96/72

textSize <- size * death_conversion_dpi96
largerSize <- (size) * death_conversion_dpi96

unmute_theme <- theme(
  text = element_text(size = textSize),#, family='Arial'),
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

## plot figo supplementary 11 
# Togliamo dal deg/dem sui m29 i due gruppi:
#   - RAD51 > 15 & POLD1 < 35
# - RAD51 < 10 & POLD1 > 40

ggplot(data=m5, aes(x=POLD1, y=RAD51_irino, color=trimmeddw3))+geom_point()+theme_bw() +
  scale_colour_gradient2(midpoint = 35, low = "blue", high = "red")+
  geom_hline(yintercept=c(15,10))+geom_vline(xintercept = c(35,40))


remove <- m5
for (i in rownames(remove)) {
  if (remove[i, "POLD1"] < 35 && remove[i, "RAD51_irino"] > 15) {
    remove[i, "remove"] <- "yes"
  } else if (remove[i, "POLD1"] > 40 && remove[i, "RAD51_irino"] < 10) {
    remove[i, "remove"] <- "yes"
  } else {
    remove[i, "remove"] <- "no"
  }
}

sd <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/samples_data"
sd <- read.table(sd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sd <- sd[!duplicated(sd$model),]

remove$model <- rownames(remove)

remove <- merge(remove, sd, by="model")
remove$type <- as.factor(remove$type)

x_breaks <- guess_ticks(remove$POLD1, fixed_min=0, fixed_max=70)
y_breaks <- guess_ticks(remove$RAD51_irino, fixed_min=-10, fixed_max=40)

ggplot(data = remove, aes(x = POLD1, y = RAD51_irino)) +
  geom_point(data = subset(remove, remove == "no"), aes(color = type), fill = "grey", shape = 21, size=4, stroke = 1) +
  geom_point(data = subset(remove, remove == "yes"), aes(color = type, fill=type), shape = 21, size = 4, stroke = 1) +
  scale_colour_manual(values = c("resistant" = rgb(red=165,green=0,blue=25, max = 255), "sensitive" = rgb(30, 85, 130, max = 255))) +
  scale_fill_manual(values = c("resistant" = rgb(red=165,green=0,blue=25, max = 255), "sensitive" = rgb(30, 85, 130, max = 255))) +
  geom_hline(yintercept = c(15, 10), linetype = "dashed") +
  geom_vline(xintercept = c(35, 40), linetype = "dashed") +
  unmute_theme+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))

remove2 <- remove
remove2 <- remove2 %>% filter(remove == "yes")

remove3 <- remove
remove3 <- remove3 %>% filter(remove =="no")

write.table(remove2, file="/scratch/trcanmed/DE_RNASeq/dataset/new_chemio_groups/removefromDEG.tsv", quote = FALSE,
            sep = "\t", col.names = TRUE, row.names = TRUE)
