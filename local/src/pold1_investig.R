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
####
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


tdata <- t(pold1_ct)
dataf <- t(t(tdata)/colSums(tdata))
colSums(dataf)# check all 1
pd <- melt(dataf)
colnames(pd) <-c('19q13.31','Irino', 'frac')

pdm <- melt(pold1_ct)
colnames(pdm) <- c('Irino', '19q13.31', 'n')

pd$id <- paste0(pd$Irino, pd$`19q13.31`)
pdm$id <- paste0(pdm$Irino, pdm$`19q13.31`)
pdm$Irino <- NULL
pdm$`19q13.31` <- NULL
pdd <- merge(pd, pdm, by="id")
pdd$recist <- factor(pdd$Irino, levels=c('PD', 'SD', 'PR'))

ggplot(data=pdd, aes(y=frac, x=factor(1),fill=as.factor(`19q13.31`)))+geom_col()+theme_bw()+facet_wrap(~Irino)+
  coord_polar(theta = "y")+scale_fill_manual(values=c('darkgrey', 'blue'))+
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))
