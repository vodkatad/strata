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


## provare GO - C5
grigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/GSEA_results_C5_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
m29 <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C5_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
nongrigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/complementare_m29_grigi/GSEA_results_C5_type_cutoff0.05-resistant.vs.sensitive.tsv', sep="\t", header=T, stringsAsFactors = F)
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



sel <- apply(pp, 1, function(x) any(x < 0.08))
top4 <- top[sel,]

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
top5 <- top4[top4$NES_grigi > 0 & top4$NES_all_29 < 0,]

pheatmap(top5,
         breaks = bk, color=my_palette)


sel <- apply(pp, 1, function(x) all(x < 0.3))
top4 <- top[sel,]

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

pheatmap(top4,
         breaks = bk, color=my_palette)

### selected pathways

grigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/m29_new_chemio_groups/risultati_gsea_totali.tsv', sep="\t", header=T)
m29 <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/risultati_gsea_totali.tsv', sep="\t", header=T)
nongrigi <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/complementare_m29_grigi/risultati_gsea_totali.tsv', sep="\t", header=T)

blu <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", 
         "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING") ##c("#74A9CF", "#3690C0", "#0570B0"), 

bordeaux <- c("HALLMARK_FATTY_ACID_METABOLISM", "REACTOME_FATTY_ACID_METABOLISM") # c("#A50026", "#D73027")

grigio <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") # black

verdi1 <- c("LIN_APC_TARGETS", "KIM_MYC_AMPLIFICATION_TARGETS_DN",
            "SANSOM_APC_TARGETS_DN", "LEF1_UP.V1_UP", "LEF1_UP.V1_DN", "BCAT.100_UP.V1_UP") #"#41AB5D", "#238B45", "#006D2C"
#c("#74C476", "#41AB5D", "#238B45")

wpath <- data.frame(name=c(blu, bordeaux, grigio, verdi1), 
                    color=c(rep('#74A9CF', length(blu)), rep("#A50026", length(bordeaux)), 'black', rep('#41AB5D', length(verdi1))),
                    order=c(rep('l1', length(blu)), rep('l2', length(bordeaux)), 'l3', rep('l4', length(verdi1))))

colnames(grigi) <- paste0(colnames(grigi), '_grigi')
colnames(nongrigi) <- paste0(colnames(nongrigi), '_non_grigi')
colnames(m29) <- paste0(colnames(m29), '_all_29')
m <- merge(grigi, m29, by='row.names')
rownames(m) <- m$Row.names
m$Row.names <- NULL
m2 <- merge(m, nongrigi, by='row.names')
rownames(m2) <- m2$Row.names
m2$Row.names <- NULL
mw <- m2[rownames(m2) %in% wpath$name,]

top <- mw[, grepl('NES', colnames(mw))]
#pval <- mw[, grepl('adjust', colnames(mw))]
pval <- mw[, grepl('pvalue', colnames(mw))]

#ns P > 0.05

#* P ≤ 0.05

#** P ≤ 0.01

#*** P ≤ 0.001

#**** P ≤ 0.0001 (For the last two choices only)
ast <- ifelse(pval < 0.05, '*', '')

minv <- min(top)
maxv <- max(top)
neutral_value <- 0
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))

rownames(wpath) <- wpath$name
wpath$name <- NULL

top <- merge(top, wpath, by="row.names")
rownames(top) <- top$Row.names
top$Row.names <- NULL
top <- top[order(top$order),]
top$color <- NULL
top$order <- NULL

ast <- merge(ast, wpath, by="row.names")
rownames(ast) <- ast$Row.names
ast$Row.names <- NULL
ast <- ast[order(ast$order),]
ast$color <- NULL
ast$order <- NULL


wpath$color <- NULL

wpath$order <- factor(wpath$order, levels = c('l1','l2','l3','l4'))

annot_color_list <- list(order=c(l1='#74A9CF', l2="#A50026", l3='black', l4='#41AB5D'))
pheatmap(top,
         breaks = bk, color=my_palette, cluster_rows = F, cluster_cols = F, annotation_row=wpath,
         annotation_colors = annot_color_list,
         display_numbers=ast, number_color='black', fontsize_number=10)


##
m29 <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/risultati_gsea_totali.tsv', sep="\t", header=T)
all <- read.table('/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/risultati_gsea_totali.tsv', sep="\t", header=T)

blu <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", 
         "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING") ##c("#74A9CF", "#3690C0", "#0570B0"), 

bordeaux <- c("HALLMARK_FATTY_ACID_METABOLISM", "REACTOME_FATTY_ACID_METABOLISM") # c("#A50026", "#D73027")

grigio <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") # black

verdi1 <- c("LIN_APC_TARGETS", "KIM_MYC_AMPLIFICATION_TARGETS_DN",
            "SANSOM_APC_TARGETS_DN", "LEF1_UP.V1_UP", "LEF1_UP.V1_DN", "BCAT.100_UP.V1_UP") #"#41AB5D", "#238B45", "#006D2C"
#c("#74C476", "#41AB5D", "#238B45")

wpath <- data.frame(name=c(blu, bordeaux, grigio, verdi1), 
                    color=c(rep('#74A9CF', length(blu)), rep("#A50026", length(bordeaux)), 'black', rep('#41AB5D', length(verdi1))),
                    order=c(rep('l1', length(blu)), rep('l2', length(bordeaux)), 'l3', rep('l4', length(verdi1))))

colnames(m29) <- paste0(colnames(m29), '_m29')
colnames(all) <- paste0(colnames(all), '_all')

m2 <- merge(m29, all, by='row.names')
rownames(m2) <- m2$Row.names
m2$Row.names <- NULL
mw <- m2[rownames(m2) %in% wpath$name,]

top <- mw[, grepl('NES', colnames(mw))]
#pval <- mw[, grepl('adjust', colnames(mw))]
pval <- mw[, grepl('pvalue', colnames(mw))]


ast <- ifelse(pval < 0.05, '*', '')


minv <- min(top)
maxv <- max(top)
neutral_value <- 0
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))

rownames(wpath) <- wpath$name
wpath$name <- NULL

top <- merge(top, wpath, by="row.names")
rownames(top) <- top$Row.names
top$Row.names <- NULL
top <- top[order(top$order),]t
top$color <- NULL
top$order <- NULL

ast <- merge(ast, wpath, by="row.names")
rownames(ast) <- ast$Row.names
ast$Row.names <- NULL
ast <- ast[order(ast$order),]
ast$color <- NULL
ast$order <- NULL


wpath$color <- NULL

wpath$order <- factor(wpath$order, levels = c('l1','l2','l3','l4'))

annot_color_list <- list(order=c(l1='#74A9CF', l2="#A50026", l3='black', l4='#41AB5D'))
pheatmap(top,
         breaks = bk, color=my_palette, cluster_rows = F, cluster_cols = F, annotation_row=wpath,
         annotation_colors = annot_color_list,
         display_numbers=ast, number_color='black', fontsize_number=10)
