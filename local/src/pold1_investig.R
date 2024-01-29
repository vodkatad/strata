library(reshape)
library(ggplot2)
targeted <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_thresholded.by_genes.txt"
targeted <- read.table(targeted, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  targeted_p <- targeted %>% filter(Gene.Symbol == "POLD1")
  targeted_p$Locus.ID <- NULL
  targeted_p$Cytoband <- NULL
  rownames(targeted_p) <- targeted_p$Gene.Symbol
  targeted_p$Gene.Symbol <- NULL
  targeted_p <- as.data.frame(t(targeted_p))
  #targeted_p <- targeted_p %>% filter(!POLD1 == "0")
  targeted_p$model <- substr(rownames(targeted_p), 1, 7)
  
  casi_f <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/w3_waterfall.tsv"
  casi <- read.table(casi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(casi) <- c('model', 'perc')
casi$folfiri.recist_3w <- ifelse(casi$perc < -50, 'PR', ifelse(casi$perc > 35, 'PD', 'SD'))

pold1 <- merge(casi, targeted_p, by = "model")
wilcox.test(pold1[pold1$POLD1==-1, 'perc'], pold1[pold1$POLD1!=-1, 'perc'])

pold1_ct <- as.matrix(table(pold1$folfiri.recist_3w, pold1$POLD1))

#del vs ndel
ct <- data.frame(pold1_ct[,1], pold1_ct[,2]+pold1_ct[,3])
ct2 <- data.frame(del=c(ct[1,1], ct[2,1]+ct[3,1]), ndel=c(ct[1,2], ct[2,2]+ct[3,2]))
rownames(ct2) <- c('PD', 'PR_SD')
fisher.test(ct2)

ct3 <- data.frame(del=c(ct[2,1], ct[1,1]+ct[3,1]), ndel=c(ct[2,2], ct[1,2]+ct[3,2]))
rownames(ct3) <- c('PR', 'PD_SD')
fisher.test(ct3)

ct4 <- data.frame(del=c(ct[2,1], ct[1,1]), ndel=c(ct[2,2], ct[1,2]))
rownames(ct4) <- c('PR', 'PD')
fisher.test(ct4)


tdata <- t(pold1_ct)
dataf <- t(t(tdata)/colSums(tdata))
colSums(dataf)# check all 1
pd <- melt(dataf)
colnames(pd) <-c('POLD1','Irino', 'frac')

pdm <- melt(pold1_ct)
colnames(pdm) <- c('Irino', 'POLD1', 'n')

pd$id <- paste0(pd$Irino, pd$POLD1)
pdm$id <- paste0(pdm$Irino, pdm$POLD1)
pdm$Irino <- NULL
pdm$POLD1 <- NULL
pdd <- merge(pd, pdm, by="id")
pdd$recist <- factor(pdd$Irino, levels=c('PD', 'SD', 'PR'))

ggplot(data=pdd, aes(y=frac, x=factor(1),fill=as.factor(POLD1)))+geom_col()+theme_bw()+facet_wrap(~Irino)+
  coord_polar(theta = "y")+scale_fill_manual(values=c('blue', 'darkgrey', 'red'))+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))


tdata <- t(ct2)
dataf <- t(t(tdata)/colSums(tdata))
colSums(dataf)# check all 1
pd <- melt(dataf)
colnames(pd) <-c('POLD1','Irino', 'frac')

ct2$id <- rownames(ct2)
pdm <- melt(ct2)
colnames(pdm) <- c('Irino', 'POLD1', 'n')

pd$id <- paste0(pd$Irino, pd$POLD1)
pdm$id <- paste0(pdm$Irino, pdm$POLD1)
pdm$Irino <- NULL
pdm$POLD1 <- NULL
pdd <- merge(pd, pdm, by="id")
pdd$recist <- factor(pdd$Irino, levels=c('PD', 'PR_SD'))

ggplot(data=pdd, aes(y=frac, x=factor(1),fill=as.factor(POLD1)))+geom_col()+theme_bw()+facet_wrap(~Irino)+
  coord_polar(theta = "y")+scale_fill_manual(values=c('blue', 'darkgrey'))+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))
########################################
  casi_f <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/w3_waterfall.tsv"
  casi <- read.table(casi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(casi) <- c('model', 'perc')
casi$folfiri.recist_3w <- ifelse(casi$perc < -50, 'PR', ifelse(casi$perc > 35, 'PD', 'SD'))

targeted_region <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_lesions.conf_90.txt"
#targeted_region <- "/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/gistic/gistic_xeno/all_lesions.conf_90.txt"

targeted_region <- read.table(targeted_region, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
targeted_region$X <- NULL #?
targeted_31 <- targeted_region %>% filter(Descriptor == "19q13.31")
targeted_31 <- targeted_31[1,]
rownames(targeted_31) <- targeted_31$Descriptor
targeted_31$Descriptor <- NULL
targeted_31$Unique.Name <- NULL
targeted_31 <- targeted_31[,-c(1:7)]
targeted_31 <- as.data.frame(t(targeted_31))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_31$model <- substr(rownames(targeted_31), 1, 7)

pold1 <- merge(casi, targeted_31, by = "model")
wilcox.test(pold1[pold1$`19q13.31`==1, 'perc'], pold1[pold1$`19q13.31`==0, 'perc'])

pold1_ct <- as.matrix(table(pold1$folfiri.recist_3w, pold1$`19q13.31`))

ct <- pold1_ct
ct2 <- data.frame(del=c(ct[1,2], ct[2,2]+ct[3,2]), ndel=c(ct[1,1], ct[2,1]+ct[3,1]))
rownames(ct2) <- c('PD', 'PR_SD')
fisher.test(ct2)

ct3 <- data.frame(del=c(ct[2,2], ct[1,2]+ct[3,2]), ndel=c(ct[2,1], ct[1,1]+ct[3,1]))
rownames(ct3) <- c('PR', 'PD_SD')
fisher.test(ct3)

ct4 <- data.frame(del=c(ct[2,2], ct[1,2]), ndel=c(ct[2,1], ct[1,1]))
rownames(ct4) <- c('PR', 'PD')
fisher.test(ct4)


tdata <- t(ct2)
dataf <- t(t(tdata)/colSums(tdata))
colSums(dataf)# check all 1
pd <- melt(dataf)
colnames(pd) <-c('19q13.31','Irino', 'frac')

ct2$id <- rownames(ct2)
pdm <- melt(ct2)
colnames(pdm) <- c('Irino', '19q13.31', 'n')

pd$id <- paste0(pd$Irino, pd$`19q13.31`)
pdm$id <- paste0(pdm$Irino, pdm$`19q13.31`)
pdm$Irino <- NULL
pdm$`19q13.31` <- NULL
pdd <- merge(pd, pdm, by="id")
pdd$recist <- factor(pdd$Irino, levels=c('PD', 'PR_SD'))

ggplot(data=pdd, aes(y=frac, x=factor(1),fill=as.factor(`19q13.31`)))+geom_col()+theme_bw()+facet_wrap(~Irino)+
coord_polar(theta = "y")+scale_fill_manual(values=c('blue', 'darkgrey'))+
geom_text(aes(label = n), position = position_stack(vjust = 0.5))+theme(legend.title=element_blank(), axis.title.y=element_blank())


### other calling methods/log2FC
call <- read.table('/mnt/trcanmed/snaketree/prj/strata/dataset/shallow/pold1_calls.tsv', sep="\t", header=TRUE)
casi_f <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/w3_waterfall.tsv"
casi <- read.table(casi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(casi) <- c('model', 'perc')
casi$folfiri.recist_3w <- ifelse(casi$perc < -50, 'PR', ifelse(casi$perc > 35, 'PD', 'SD'))

call$chromosome <- NULL; call$start <- NULL; call$end <- NULL
call <- t(call)

table(apply(call, 1, function(x){length(unique(x))==1}))# check all ==

call_p <- data.frame(model=substr(rownames(call), 0, 7), pold1=call[,1])
m <- merge(call_p, casi, by="model")
stopifnot(nrow(m) == 84)

t1 <- as.matrix(table(m$folfiri.recist_3w,m$pold1))
ct <- data.frame(t1[, 1], t1[,2]+t1[,3])

ct2 <- data.frame(del=c(ct[1,1], ct[1,2]+ct[1,3]), ndel=c(ct[1,2], ct[2,2]+ct[3,2]))
rownames(ct2) <- c('PD', 'PR_SD')
fisher.test(ct2)

call <- read.table('/mnt/trcanmed/snaketree/prj/strata/dataset/shallow/pold1_log2.tsv', sep="\t", header=TRUE)
casi_f <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/w3_waterfall.tsv"
casi <- read.table(casi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(casi) <- c('model', 'perc')
casi$folfiri.recist_3w <- ifelse(casi$perc < -50, 'PR', ifelse(casi$perc > 35, 'PD', 'SD'))

call$chromosome <- NULL; call$start <- NULL; call$end <- NULL
call <- t(call)

call <- rowMeans(call)

call_p <- data.frame(model=substr(names(call), 0, 7), pold1=call)
m <- merge(call_p, casi, by="model")
stopifnot(nrow(m) == 84)
m$folfiri.recist_3w <- factor(m$folfiri.recist_3w, levels=c('PR', 'SD', 'PD'))
ggplot(data=m, aes(x=folfiri.recist_3w, y=pold1))+geom_boxplot(outlier.shape = NA)+geom_jitter(height=NULL)+theme_bw()
ggplot(data=m, aes(x=perc, y=pold1))+geom_point()+theme_bw()


### pold1 vs locus
targeted <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_thresholded.by_genes.txt"
targeted <- read.table(targeted, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

targeted_p <- targeted %>% filter(Gene.Symbol == "POLD1")
targeted_p$Locus.ID <- NULL
targeted_p$Cytoband <- NULL
rownames(targeted_p) <- targeted_p$Gene.Symbol
targeted_p$Gene.Symbol <- NULL
targeted_p <- as.data.frame(t(targeted_p))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_p$model <- substr(rownames(targeted_p), 1, 7)
targeted_region <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_lesions.conf_90.txt"
#targeted_region <- "/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/gistic/gistic_xeno/all_lesions.conf_90.txt"

targeted_region <- read.table(targeted_region, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
targeted_region$X <- NULL #?
targeted_31 <- targeted_region %>% filter(Descriptor == "19q13.31")
targeted_31 <- targeted_31[1,]
rownames(targeted_31) <- targeted_31$Descriptor
targeted_31$Descriptor <- NULL
targeted_31$Unique.Name <- NULL
targeted_31 <- targeted_31[,-c(1:7)]
targeted_31 <- as.data.frame(t(targeted_31))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_31$model <- substr(rownames(targeted_31), 1, 7)

m <- merge(targeted_p, targeted_31, by="model")
stopifnot(nrow(m) == 84)
table(m$POLD1, m$`19q13.31`)



# expr
ex <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMX_BASALE-POLD1_fpkm_ave.tsv', sep="\t", header=TRUE, stringsAsFactors = F)

targeted <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_thresholded.by_genes.txt"
targeted <- read.table(targeted, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

targeted_p <- targeted %>% filter(Gene.Symbol == "POLD1")
targeted_p$Locus.ID <- NULL
targeted_p$Cytoband <- NULL
rownames(targeted_p) <- targeted_p$Gene.Symbol
targeted_p$Gene.Symbol <- NULL
targeted_p <- as.data.frame(t(targeted_p))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_p$model <- substr(rownames(targeted_p), 1, 7)

m <- merge(targeted_p, ex, by="model")

ggplot(data=m, aes(x=as.factor(POLD1), y=log(expr+1)))+geom_boxplot(outlier.shape = NA)+geom_jitter(height=NULL)+theme_bw()

ex <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMX_BASALE-POLD1_fpkm_ave.tsv', sep="\t", header=TRUE, stringsAsFactors = F)

targeted_region <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_lesions.conf_90.txt"
#targeted_region <- "/mnt/trcanmed/snaketree/prj/biobanca/local/share/data/shallowseq/gistic/gistic_xeno/all_lesions.conf_90.txt"

targeted_region <- read.table(targeted_region, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
targeted_region$X <- NULL #?
targeted_31 <- targeted_region %>% filter(Descriptor == "19q13.31")
targeted_31 <- targeted_31[1,]
rownames(targeted_31) <- targeted_31$Descriptor
targeted_31$Descriptor <- NULL
targeted_31$Unique.Name <- NULL
targeted_31 <- targeted_31[,-c(1:7)]
targeted_31 <- as.data.frame(t(targeted_31))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_31$model <- substr(rownames(targeted_31), 1, 7)

m <- merge(targeted_31, ex, by="model")

ggplot(data=m, aes(x=as.factor(`19q13.31`), y=log(expr+1)))+geom_boxplot(outlier.shape = NA)+geom_jitter(height=NULL)+theme_bw()


## 2 metodi di chiamata per pold1 gistic cghall
targeted <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_thresholded.by_genes.txt"
targeted <- read.table(targeted, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

targeted_p <- targeted %>% filter(Gene.Symbol == "POLD1")
targeted_p$Locus.ID <- NULL
targeted_p$Cytoband <- NULL
rownames(targeted_p) <- targeted_p$Gene.Symbol
targeted_p$Gene.Symbol <- NULL
targeted_p <- as.data.frame(t(targeted_p))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_p$model <- substr(rownames(targeted_p), 1, 7)

call <- read.table('/mnt/trcanmed/snaketree/prj/strata/dataset/shallow/pold1_calls.tsv', sep="\t", header=TRUE)
call$chromosome <- NULL; call$start <- NULL; call$end <- NULL
call <- t(call)

table(apply(call, 1, function(x){length(unique(x))==1}))# check all ==

call_p <- data.frame(model=substr(rownames(call), 0, 7), pold1=call[,1])

m <- merge(targeted_p, call_p, by="model")
table(m$POLD1, m$pold1)


# vs expr + casi
ex <- read.table('/mnt/trcanmed/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/LMX_BASALE-POLD1_fpkm_ave.tsv', sep="\t", header=TRUE, stringsAsFactors = F)

targeted <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_thresholded.by_genes.txt"
targeted <- read.table(targeted, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

targeted_p <- targeted %>% filter(Gene.Symbol == "POLD1")
targeted_p$Locus.ID <- NULL
targeted_p$Cytoband <- NULL
rownames(targeted_p) <- targeted_p$Gene.Symbol
targeted_p$Gene.Symbol <- NULL
targeted_p <- as.data.frame(t(targeted_p))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_p$model <- substr(rownames(targeted_p), 1, 7)

m <- merge(targeted_p, ex, by="model")


casi_f <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/w3_waterfall.tsv"
casi <- read.table(casi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(casi) <- c('model', 'perc')
casi$folfiri.recist_3w <- ifelse(casi$perc < -50, 'PR', ifelse(casi$perc > 35, 'PD', 'SD'))

m2 <- merge(m, casi, by= "model")
m2$recist <- factor(m2$folfiri.recist_3w, levels=c('PD', 'SD', 'PR'))

ggplot(data=m2, aes(x=as.factor(POLD1), y=log(expr+1)))+geom_boxplot(outlier.shape = NA)+geom_jitter(aes(color=recist),height=NULL)+theme_bw()+scale_color_manual(values=c('red', 'darkgrey', 'green'))

wilcox.test(m2[m2$POLD1==-1, 'expr'], m2[m2$POLD1==0, 'expr'])
wilcox.test(m2[m2$POLD1==-1, 'expr'], m2[m2$POLD1==1, 'expr'])

################################################################################# inverted pies
library(reshape)
library(ggplot2)
targeted <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_thresholded.by_genes.txt"
targeted <- read.table(targeted, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

targeted_p <- targeted %>% filter(Gene.Symbol == "POLD1")
targeted_p$Locus.ID <- NULL
targeted_p$Cytoband <- NULL
rownames(targeted_p) <- targeted_p$Gene.Symbol
targeted_p$Gene.Symbol <- NULL
targeted_p <- as.data.frame(t(targeted_p))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_p$model <- substr(rownames(targeted_p), 1, 7)

casi_f <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/w3_waterfall.tsv"
casi <- read.table(casi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(casi) <- c('model', 'perc')
casi$folfiri.recist_3w <- ifelse(casi$perc < -50, 'PR', ifelse(casi$perc > 35, 'PD', 'SD'))

pold1 <- merge(casi, targeted_p, by = "model")
wilcox.test(pold1[pold1$POLD1==-1, 'perc'], pold1[pold1$POLD1!=-1, 'perc'])

pold1_ct <- as.matrix(table(pold1$folfiri.recist_3w, pold1$POLD1))
ct <- data.frame(pold1_ct[, 1], pold1_ct[,2]+pold1_ct[,3])


dataf <- t(t(ct)/colSums(ct))
colSums(dataf)# check all 1
colnames(dataf) <- c('Single copy loss', 'Not deleted')
pd <- melt(dataf)
colnames(pd) <-c('Irino', 'POLD1','frac')

colnames(ct) <- c('Single copy loss', 'Not deleted')
ct$id <- rownames(ct)
pdm <- melt(ct)
colnames(pdm) <- c('Irino', 'POLD1', 'n')

pd$id <- paste0(pd$Irino, pd$POLD1)
pdm$id <- paste0(pdm$Irino, pdm$POLD1)
pdm$Irino <- NULL
pdm$POLD1 <- NULL
pdd <- merge(pd, pdm, by="id")
pdd$recist <- factor(pdd$Irino, levels=c('PD', 'SD', 'PR'))

ggplot(data=pdd, aes(y=frac, x=factor(1),fill=as.factor(recist)))+geom_col()+theme_bw()+facet_wrap(~POLD1)+
  coord_polar(theta = "y")+scale_fill_manual(values=c('red', 'darkgrey', 'blue'))+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))
################################################################################# inverted pies
# 6 weeks
library(reshape)
library(ggplot2)
targeted <- "/mnt/trcanmed/snaketree/prj/strata/local/share/data/gistic/all/558381/all_thresholded.by_genes.txt"
targeted <- read.table(targeted, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

targeted_p <- targeted %>% filter(Gene.Symbol == "POLD1")
targeted_p$Locus.ID <- NULL
targeted_p$Cytoband <- NULL
rownames(targeted_p) <- targeted_p$Gene.Symbol
targeted_p$Gene.Symbol <- NULL
targeted_p <- as.data.frame(t(targeted_p))
#targeted_p <- targeted_p %>% filter(!POLD1 == "0")
targeted_p$model <- substr(rownames(targeted_p), 1, 7)

casi_f <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/w6_waterfall.tsv"
casi <- read.table(casi_f, quote = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

colnames(casi) <- c('model', 'perc')
casi$folfiri.recist_3w <- ifelse(casi$perc < -50, 'PR', ifelse(casi$perc > 35, 'PD', 'SD'))

pold1 <- merge(casi, targeted_p, by = "model")
wilcox.test(pold1[pold1$POLD1==-1, 'perc'], pold1[pold1$POLD1!=-1, 'perc'])

pold1_ct <- as.matrix(table(pold1$folfiri.recist_3w, pold1$POLD1))

#del vs ndel
ct <- data.frame(pold1_ct[,1], pold1_ct[,2]+pold1_ct[,3])
ct2 <- data.frame(del=c(ct[1,1], ct[2,1]+ct[3,1]), ndel=c(ct[1,2], ct[2,2]+ct[3,2]))
rownames(ct2) <- c('PD', 'PR_SD')
fisher.test(ct2)

ct3 <- data.frame(del=c(ct[2,1], ct[1,1]+ct[3,1]), ndel=c(ct[2,2], ct[1,2]+ct[3,2]))
rownames(ct3) <- c('PR', 'PD_SD')
fisher.test(ct3)

ct4 <- data.frame(del=c(ct[2,1], ct[1,1]), ndel=c(ct[2,2], ct[1,2]))
rownames(ct4) <- c('PR', 'PD')
fisher.test(ct4)


tdata <- t(pold1_ct)
dataf <- t(t(tdata)/colSums(tdata))
colSums(dataf)# check all 1
pd <- melt(dataf)
colnames(pd) <-c('POLD1','Irino', 'frac')

pdm <- melt(pold1_ct)
colnames(pdm) <- c('Irino', 'POLD1', 'n')

pd$id <- paste0(pd$Irino, pd$POLD1)
pdm$id <- paste0(pdm$Irino, pdm$POLD1)
pdm$Irino <- NULL
pdm$POLD1 <- NULL
pdd <- merge(pd, pdm, by="id")
pdd$recist <- factor(pdd$Irino, levels=c('PD', 'SD', 'PR'))

ggplot(data=pdd, aes(y=frac, x=factor(1),fill=as.factor(POLD1)))+geom_col()+theme_bw()+facet_wrap(~Irino)+
  coord_polar(theta = "y")+scale_fill_manual(values=c('blue', 'darkgrey', 'red'))+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))


tdata <- t(ct2)
dataf <- t(t(tdata)/colSums(tdata))
colSums(dataf)# check all 1
pd <- melt(dataf)
colnames(pd) <-c('POLD1','Irino', 'frac')

ct2$id <- rownames(ct2)
pdm <- melt(ct2)
colnames(pdm) <- c('Irino', 'POLD1', 'n')

pd$id <- paste0(pd$Irino, pd$POLD1)
pdm$id <- paste0(pdm$Irino, pdm$POLD1)
pdm$Irino <- NULL
pdm$POLD1 <- NULL
pdd <- merge(pd, pdm, by="id")
pdd$recist <- factor(pdd$Irino, levels=c('PD', 'PR_SD'))

ggplot(data=pdd, aes(y=frac, x=factor(1),fill=as.factor(POLD1)))+geom_col()+theme_bw()+facet_wrap(~Irino)+
  coord_polar(theta = "y")+scale_fill_manual(values=c('blue', 'darkgrey'))+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))