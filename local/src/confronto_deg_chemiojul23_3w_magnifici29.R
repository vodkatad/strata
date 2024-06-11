### confronto deg chemio 3w e magnifici29

library(tidyverse)
library(ggplot2)
library(ggvenn)

w3_f <- snakemake@input[["deg_3"]]
m29_f <- snakemake@input[["deg_select"]]
#h_w3_f <- snakemake@input[["gsea_h_w3"]]
#h_m29_f <- snakemake@input[["gsea_h_m29"]]
#c2_w3_f <- snakemake@input[["gsea_c2_w3"]]
#c2_m29_f <- snakemake@input[["gsea_c2_m29"]]
#c6_w3_f <- snakemake@input[["gsea_c6_w3"]]
#c6_m29_f <- snakemake@input[["gsea_c6_m29"]]
plot_cor_png <- snakemake@output[["png"]]
#plot_cor_eps <- snakemake@output[["eps"]]
#cor_gsea_h <- snakemake@output[["cor_h_tsv"]]
#cor_gsea_c2 <- snakemake@output[["cor_c2_tsv"]]
#cor_gsea_c6 <- snakemake@output[["cor_c6_tsv"]]
log_f <- snakemake@log[['log']]

#load('/scratch/trcanmed/AF_spectra/dataset_Figures_Tables/theme_5.Rdata')
size <- 8

#font_add(family = "myriad", regular = snakemake@input[['myriad']])
#showtext_auto()

# Da Marti e https://www.christophenicault.com/post/understand_size_dimension_ggplot2/
# showtext_opts(dpi = 300) 
# since we are not changing fonts in the end cause myriad end up not being text object I'm not sure it's needed
# showtext_auto(enable = TRUE)

#textSize <- textSize * (96/72) # these conversion were needed because the default dpi for text was 96?
# in the svg the number passed to theme was reported as size = ..px.. rather than pt (?)
#largerSize <- largerSize * (96/72) 
death_conversion_dpi96 = 96/72

textSize <- size * death_conversion_dpi96
largerSize <- size* death_conversion_dpi96

unmute_theme <- theme(
  text = element_text(size = textSize, family='Arial'),
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
#axis.ticks.length= unit(1.905*death_conversion_dpi96, "mm"),


# function that given values to be plotted on an axis will return:
# vector of breaks, trying to guess which max will be the best one
# this will be used as scale_y_continuous(breaks=  and as ylim(min, max) to have the - also limits-c()
# last tick at the extremity of the axis.
# other parameter is n. of ticks
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

#w3_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/type_cutoff0.05-non_responder_3Q.vs.responder_1Q.deseq2.tsv"
w3 <- read.table(w3_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(w3) <- paste0(colnames(w3), "_w3")
w3$geni <- rownames(w3)

#m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/type_cutoff0.05-resistant.vs.sensitive.deseq2.tsv"
m29 <- read.table(m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(m29) <- paste0(colnames(m29), "_m29")
m29$geni <- rownames(m29)

merged <- merge(w3, m29, by="geni")

sink(log_f, append=TRUE)
"Correlation logFC w3 and magnifici29"
cor.test(merged$log2FoldChange_w3, merged$log2FoldChange_m29)
sink()
# Pearson's product-moment correlation
# 
# data:  merged$log2FoldChange_w3 and merged$log2FoldChange_w6
# t = 92.66, df = 17904, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5593329 0.5791338
# sample estimates:
#       cor 
# 0.5693159 

y_breaks <- guess_ticks(merged$log2FoldChange_w3, fixed_min=-2, fixed_max=3)
x_breaks <- guess_ticks(merged$log2FoldChange_w6, fixed_min=-2, fixed_max=3)

p <- ggplot(merged, aes(x=log2FoldChange_w3, y=log2FoldChange_m29)) + 
  geom_point(size = 1)+
  geom_smooth(method=lm, se=FALSE)+ # ratio is ... 1 becomes 0.939
  unmute_theme+theme(legend.position="none") + xlab('log2FoldChange_w3')+ylab('log2FoldChange_m29')+
  scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
  scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))

ggsave(p, file=plot_cor_png, width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")

# p <- ggplot(merged, aes(x=log2FoldChange_w3, y=log2FoldChange_m29)) + 
#   geom_point(size = 1)+
#   geom_smooth(method=lm, se=FALSE)+ # ratio is ... 1 becomes 0.939
#   unmute_theme+theme(legend.position="none")+theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ 
#   theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
#   scale_y_continuous(breaks=y_breaks, limits=c(min(y_breaks),max(y_breaks)), expand=c(0,0))+
#   scale_x_continuous(breaks=x_breaks, limits=c(min(x_breaks),max(x_breaks)), expand=c(0,0))

# ggsave(p, file=plot_cor_eps, width=89*(death_conversion_dpi96), height=89*(death_conversion_dpi96), units="mm")

#h_w3_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_H_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv"
# h_w3 <- read.table(h_w3_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# #h_m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_H_type_cutoff0.05-resistant.vs.sensitive.tsv"
# h_m29 <- read.table(h_m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# # h_w3 <- h_w3 %>% filter(p.adjust < 0.05)
# # 
# # h_m29 <- h_m29 %>% filter(p.adjust < 0.05)
# ws <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE",
#         "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
#         "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_MYC_TARGETS_V1","HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN")

# h_w3 <- h_w3 %>% filter(ID %in% ws)
# h_m29 <- h_m29 %>% filter(ID %in% ws)
# names(h_w3)[names(h_w3) == 'enrichmentScore'] <- "enrichmentScore_w3"
# names(h_m29)[names(h_m29) == 'enrichmentScore'] <- "enrichmentScore_m29"

# merged_h <- merge(h_w3, h_m29, by="ID")
# merged_h <- merged_h[c("ID", "enrichmentScore_w3", "enrichmentScore_m29")]

# sink(log_f, append=TRUE)
# "Correlation enrichment Score gsea H"
# cor.test(merged_h$enrichmentScore_w3, merged_h$enrichmentScore_m29)
# sink()

# # Pearson's product-moment correlation
# # 
# # data:  merged_h$enrichmentScore_w3 and merged_h$enrichmentScore_m29
# # t = 2.4594, df = 6, p-value = 0.04916
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.00772091 0.94258852
# # sample estimates:
# #       cor 
# # 0.7085389 

# write.table(merged_h, file=cor_gsea_h, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# #c2_w3_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_C2_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv"
# c2_w3 <- read.table(c2_w3_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)  

# #c2_m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C2_type_cutoff0.05-resistant.vs.sensitive.tsv"
# c2_m29 <- read.table(c2_m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# # c2_w3 <- c2_w3 %>% filter(p.adjust < 0.05)
# # 
# # c2_m29 <- c2_m29 %>% filter(p.adjust < 0.05)

# ws <- c("WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA",
#         "REACTOME_FATTY_ACID_METABOLISM","KEGG_FATTY_ACID_METABOLISM",
#         "WP_TCA_CYCLE_AKA_KREBS_OR_CITRIC_ACID_CYCLE","REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
#         "REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION","KIM_MYC_AMPLIFICATION_TARGETS_DN",
#         "LIN_APC_TARGETS","SANSOM_APC_TARGETS_DN")
 
# c2_w3 <- c2_w3 %>% filter(ID %in% ws)
# c2_m29 <- c2_m29 %>% filter(ID %in% ws)
# names(c2_w3)[names(c2_w3) == 'enrichmentScore'] <- "enrichmentScore_w3"
# names(c2_m29)[names(c2_m29) == 'enrichmentScore'] <- "enrichmentScore_m29"
 
# merged_c2 <- merge(c2_w3, c2_m29, by="ID")
# merged_c2 <- merged_c2[c("ID", "enrichmentScore_w3", "enrichmentScore_m29")]
 
# write.table(merged_c2, file=cor_gsea_c2, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
# sink(log_f, append=TRUE)
# "Correlation enrichment Score gsea c2"
# cor.test(merged_c2$enrichmentScore_w3, merged_c2$enrichmentScore_m29)
# sink()
# # Pearson's product-moment correlation
# # 
# # data:  merged_c2$enrichmentScore_w3 and merged_c2$enrichmentScore_m29
# # t = 2.4642, df = 7, p-value = 0.0432
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.03184686 0.92636879
# # sample estimates:
# #       cor 
# # 0.6815536


# #c6_w3_f <- "/scratch/trcanmed/DE_RNASeq/dataset/chemio_jul23/GSEA_results_C6_type_cutoff0.05-non_responder_3Q.vs.responder_1Q.tsv"
# c6_w3 <- read.table(c6_w3_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)  

# #c6_m29_f <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/GSEA_results_C6_type_cutoff0.05-resistant.vs.sensitive.tsv"
# c6_m29 <- read.table(c6_m29_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# #c6_w3 <- c6_w3 %>% filter(p.adjust < 0.05)
 
# #c6_m29 <- c6_m29 %>% filter(p.adjust < 0.05)
# ws <- c("LEF1_UP.V1_UP","LEF1_UP.V1_DN","BCAT.100_UP.V1_UP")
 
# c6_w3 <- c6_w3 %>% filter(ID %in% ws)
# c6_m29 <- c6_m29 %>% filter(ID %in% ws)
# names(c6_w3)[names(c6_w3) == 'enrichmentScore'] <- "enrichmentScore_w3"
# names(c6_m29)[names(c6_m29) == 'enrichmentScore'] <- "enrichmentScore_m29"
 
# merged_c6 <- merge(c6_w3, c6_m29, by="ID")
# merged_c6 <- merged_c6[c("ID", "enrichmentScore_w3", "enrichmentScore_m29")]
 
# write.table(merged_c6, file=cor_gsea_c6, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
# sink(log_f, append=TRUE)
# "Correlation enrichment Score gsea C6"
# cor.test(merged_c6$enrichmentScore_w3, merged_c6$enrichmentScore_m29)
# sink()
# # Pearson's product-moment correlation
# # 
# # data:  merged_c6$enrichmentScore_w3 and merged_c6$enrichmentScore_m29
# # t = 11.103, df = 1, p-value = 0.05719
# # alternative hypothesis: true correlation is not equal to 0
# # sample estimates:
# #       cor 
# # 0.9959683

# # sign_w3 <- rbind(h_w3, c2_w3, c6_w3)
# # sign_w3_v <- sign_w3$Description
# # 
# # sign_m29 <- rbind(h_m29, c2_m29, c6_m29)
# # sign_m29_v <- sign_m29$Description
# # 
# # list_sign_gsea <- list(w3 = sign_w3_v, m29 = sign_m29_v)
# # 
# # length(intersect(sign_w3_v, sign_m29_v))
# # 
# # ggvenn(list_sign_gsea, c("w3", "m29")) 
