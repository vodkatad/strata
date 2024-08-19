## deg chemio aug24 removal groups

library(tidyverse)


metadata_o_f <- snakemake@input[["metadata"]] 
casi_f <- snakemake@input[["folfiri"]]
remove_f <- snakemake@input[["caseremove"]]
meta <- snakemake@output[["sample"]]


#casi_f <- "/scratch/trcanmed/DE_RNASeq/local/share/data/chemio_def_jul23/CHEMIO_WATERFALL_PLOT_Eugy_Luglio2023.tsv"
casi <- read.table(casi_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
casi$ff30 <- NULL
casi <- casi[order(casi$X3WKS, decreasing = TRUE),]
casi$ntile <- ntile(casi$X3WKS, 4)
casi$ntile3 <- ntile(casi$X3WKS, 3)

ggplot(casi, aes(x = reorder(CASE, -X3WKS), y = X3WKS, fill = as.factor(ntile))) + geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, size = 5))

ggplot(casi, aes(x = reorder(CASE, -X3WKS), y = X3WKS, fill = ntile3)) + geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, size = 5))

quarti <- quantile(casi$X3WKS, c(1/3, 2/3))

casi$quartile <- NA

casi$quartile <- ifelse(casi$X3WKS < quarti[1], 1, ifelse(casi$X3WKS > quarti[2], 3, 2))                

casi <- casi %>% filter(quartile == 1 | quartile == 3)

#remove_f <- "/scratch/trcanmed/DE_RNASeq/dataset/new_chemio_groups/removefromDEG.tsv"
remove <- read.table(remove_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
remove <- rownames(remove)

casi <- casi %>% filter(!CASE %in% remove)

save.image('pippo.Rdata')

#metadata_o_f <- "/scratch/trcanmed/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meda_f <- read.table(metadata_o_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meda_f$RNA_marker <- NULL
meda_f$RNA_PC <- NULL
meda_f$METHYL_L <- NULL
meda_f$FRA_L <- NULL
meda_f$w3_cetuxi <- NULL
meda_f$w3_irino <- NULL

meda_f <- filter(meda_f, grepl("LMX_BASALE", type))
meda_f$CASE <- substr(meda_f$sample_id_R, 1,7)
meda_f <- meda_f %>% mutate(type = gsub(".1", "", type))
meda_f <- meda_f %>% mutate(sample_id_R = gsub("-2", ".2", sample_id_R))

merged <- merge(casi, meda_f, by = "CASE")
res <- as.data.frame(merged[, c(1, 8, 9, 10)])
rownames(res) <- res$sample_id_R
res$sample_id_R <- NULL
names(res)[names(res) == "CASE"] <- "sample"

for (i in seq(rownames(res))){
  if (res[i, "quartile"] == 1) {
    res[i, "quartile"] <- "responder_1Q"
  } else {
    res[i, "quartile"] <- "non_responder_3Q"
  }
}

names(res)[names(res) == "quartile"] <- "type"

### rimozioni 
### CRC0578 perchè tette
res <- res %>% filter(!sample == "CRC0578")

write.table(res, file = meta, quote = FALSE, sep = "\t", col.names = TRUE)
