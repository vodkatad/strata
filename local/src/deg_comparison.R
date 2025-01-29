library(pheatmap)
library(ggplot2)
grigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_H_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T)
m29 <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_H_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T)
nongrigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/complementare_m29_grigi/GSEA_results_H_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T)
colnames(grigi) <- paste0(colnames(grigi), '_grigi')
colnames(nongrigi) <- paste0(colnames(nongrigi), '_non_grigi')
colnames(m29) <- paste0(colnames(m29), '_all_29')

m <- merge(grigi, m29, by="row.names")
dim(m)
cor.test(m$NES_grigi, m$NES_all_29)
plot(m$NES_grigi, m$NES_all_29)
m[m$NES.x > 0 & m$NES.y < 0,]
m[m$NES.x > 0 & m$NES.y < 0, c('NES.x', 'NES.x', 'qvalue.x', 'qvalue.y')]
m[m$NES.x > 0 & m$NES.y < 0, c('NES.x', 'NES.x', 'qvalues.x', 'qvalues.y')]
m[m$NES.x > 0 & m$NES.y < 0, c('ID.x', 'NES.x', 'NES.x', 'qvalues.x', 'qvalues.y')]
history()
m[m$NES.x > 0 & m$NES.y < 0, c('ID.x', 'NES.x', 'NES.y', 'qvalues.x', 'qvalues.y')]


ggplot(data=m, aes(x=NES.x, y=NES.y, color=-log10(p.adjust.x)))+geom_point()+theme_bw()
ggplot(data=m, aes(x=NES.x, y=NES.y, color=-log10(p.adjust.y)))+geom_point()+theme_bw()

rownames(m) <- m$Row.names
m$Row.names <- NULL
m3 <- merge(m, nongrigi, by="row.names")


top <- m3[, c('NES_all_29', 'NES_grigi', 'NES_non_grigi')]
rownames(top) <- m3$Row.names
rownames(top) <- gsub('HALLMARK_', '', rownames(top))

minv <- min(top)
maxv <- max(top)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))

pheatmap(top,
         breaks = bk, color=my_palette)

pp <- m3[, c('p.adjust_all_29', 'p.adjust_grigi', 'p.adjust_non_grigi')]
sel <- apply(pp, 1, function(x) any(x < 0.1))


top2 <- top[sel,]
pheatmap(top2,
         breaks = bk, color=my_palette)

sel <- apply(pp, 1, function(x) all(x < 0.3))


top2 <- top[sel,]
pheatmap(top2,
         breaks = bk, color=my_palette)
###
grigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T)
m29 <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T)
m <- merge(grigi, m29, by="row.names")
dim(m)
dim(m29)
cor.test(m$NES.x, m$NES.y)
plot(m$NES.x, m$NES.y)

ggplot(data=m, aes(x=NES.x, y=NES.y, color=-log10(p.adjust.x)))+geom_point()+theme_bw()
ggplot(data=m, aes(x=NES.x, y=NES.y, color=-log10(p.adjust.y)))+geom_point()+theme_bw()


grigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_C6_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T)
m29 <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C6_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T)
m <- merge(grigi, m29, by="row.names")
dim(m)
dim(m29)
cor.test(m$NES.x, m$NES.y)
plot(m$NES.x, m$NES.y)
ggplot(data=m, aes(x=NES.x, y=NES.y, color=-log10(p.adjust.x)))+geom_point()+theme_bw()
ggplot(data=m, aes(x=NES.x, y=NES.y, color=-log10(p.adjust.y)))+geom_point()+theme_bw()

### C2 reactome
grigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
m29 <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
nongrigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/complementare_m29_grigi/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
colnames(grigi) <- paste0(colnames(grigi), '_grigi')
colnames(nongrigi) <- paste0(colnames(nongrigi), '_non_grigi')
colnames(m29) <- paste0(colnames(m29), '_all_29')

m <- merge(grigi, m29, by="row.names")
dim(m)
cor.test(m$NES_grigi, m$NES_all_29)

rownames(m) <- m$Row.names
m$Row.names <- NULL
m3 <- merge(m, nongrigi, by="row.names")


top <- m3[, c('NES_all_29', 'NES_grigi', 'NES_non_grigi'),]
pp <- m3[, c('p.adjust_all_29', 'p.adjust_grigi', 'p.adjust_non_grigi')]

rownames(top) <- m3$Row.names
rownames(pp) <- m3$Row.names

top <- top[grepl('REACTOME', rownames(top)),]
pp <- pp[grepl('REACTOME', rownames(pp)),]

rownames(top) <- gsub('REACTOME', '', rownames(top))

sel <- apply(pp, 1, function(x) any(x < 0.03))


top2 <- top[sel,]
minv <- min(top2)
maxv <- max(top2)
neutral_value <- 0
bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))

pheatmap(top2,
         breaks = bk, color=my_palette)

sel <- apply(pp, 1, function(x) all(x < 0.2))


top3 <- top2[top2$NES_grigi > 0 & top2$NES_all_29 < 0,]
pheatmap(top3,
         breaks = bk, color=my_palette)

sel <- apply(pp, 1, function(x) all(x < 0.2))
top4 <- top[sel,]
pheatmap(top4,
         breaks = bk, color=my_palette)

# c2 all
grigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
m29 <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
nongrigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/complementare_m29_grigi/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
colnames(grigi) <- paste0(colnames(grigi), '_grigi')
colnames(nongrigi) <- paste0(colnames(nongrigi), '_non_grigi')
colnames(m29) <- paste0(colnames(m29), '_all_29')

m <- merge(grigi, m29, by="row.names")
dim(m)
cor.test(m$NES_grigi, m$NES_all_29)

rownames(m) <- m$Row.names
m$Row.names <- NULL
m3 <- merge(m, nongrigi, by="row.names")


top <- m3[, c('NES_all_29', 'NES_grigi', 'NES_non_grigi'),]
pp <- m3[, c('p.adjust_all_29', 'p.adjust_grigi', 'p.adjust_non_grigi')]

rownames(top) <- m3$Row.names
rownames(pp) <- m3$Row.names



sel <- apply(pp, 1, function(x) all(x < 0.2))
top4 <- top[sel,]
top5 <- top4[top4$NES_grigi > 0 & top4$NES_all_29 < 0,]

pheatmap(top5,
         breaks = bk, color=my_palette)


## provare GO?
