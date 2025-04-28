d <- read.table('/mnt/cold1/snaketree/prj/strata/local/share/data/qPCR_data_dCT.txt', sep="\t", header=T)
rownames(d) <- d$id
d$id <- NULL
info <- d[, c('response'), drop=FALSE]
d$response <- NULL


basali <- d[grepl('NT', rownames(d)),]
info_b <- info[grepl('NT', rownames(info)),, drop=F]

aveR <- colMeans(basali[info_b$response=="RESISTANT",])
aveS <- colMeans(basali[info_b$response=="SENSITIVE",])


ddCT <- aveR-aveS
library(pheatmap)
pheatmap(basali, annotation_row = info_b)

g <- read.table('/mnt/cold1/snaketree/prj/strata/local/share/data//m29_LFC_qPCR_selection.tsv', sep="\t", header=T)
rownames(g) <- g$genes
g$genes <- NULL
g$ann <- ifelse(g$log2FoldChange_m29 > 0, 'up', 'down')



pheatmap(basali, annotation_row = info_b, annotation_col=g)

basali <- basali[order(info_b$response),]
pheatmap(basali, annotation_row = info_b, annotation_col=g, cluster_rows = F)


dd <- as.data.frame(ddCT)
m <- merge(dd, g,  by="row.names")
ggplot(data=m, aes(x=ddCT, y=log2FoldChange_m29))+geom_point()+geom_smooth(method="lm")
cor.test(m$ddCT, m$log2FoldChange_m29)
