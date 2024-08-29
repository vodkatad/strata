#!/usr/bin/env Rscript

library(getopt)

set.seed(42)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'mvalues_in', 'm', 1, 'character',
  'info_in', 'i', 1, 'character',
  'formula_in', 'f', 1, 'character',
  'contrast_in', 'c', 1, 'character',
  'batch_in', 'b', 1, 'character',
  'output', 'o', 1, 'character',
  'output_summ', 's', 1, 'character',
  'output_volcano', 'v', 1, 'character',
  'alpha', 'a', 1, 'numeric',
  'lfc', 'l', 1, 'numeric',
  'what', 'w', 1, 'character',
  'image_out', 'r', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)

library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)

contrast <- opt$contrast_in
print(contrast)

if ( opt$batch_in=='off' ) {
  formula <- paste0(opt$formula_in,' + 0')
} else {
  formula <- paste0(opt$formula_in,' + ',opt$batch_in,' + 0')
}
print(formula)

mvalues <- read.table(gzfile(opt$mvalues_in), sep='\t', quote="", row.names=1, header=TRUE)
print("expression file read")
# mvalues = log2(mvalues+0.1) 
# names(mvalues) <- gsub(".", "-", names(mvalues), fixed=TRUE)
info <- read.table(opt$info_in, sep='\t', quote="", row.names=2, header=TRUE)
mvalues <- mvalues[,rownames(info)]
# info$cluster <- as.character(info$cluster)


if (! all(rownames(info)==colnames(mvalues))) {
    stop("Cannot sort rows of info according to m values!")
}

f1 <- as.formula(formula)
terms <- attr(terms(f1),"term.labels")
if (! all(terms %in% names(info))) { stop(paste0("Cannot find ",terms," in info colnames: please use a known colname")) }
f2 <- as.formula(paste0("~ info$",terms," + 0"))

terms_sum <- paste0("info$",terms) # Add info$ ai termini passati
final_form <- paste0("~ ",paste0(terms_sum, collapse = "+")," + 0") # Concatena col + i termine e aggiunge davanti la ~
f2 <- as.formula(final_form) # Crea formula finale

design <- model.matrix(f2, mvalues)

contrast_obj <- unlist(strsplit(contrast, "-")) # Splitta il contrasto passato per avere i singoli caratteri
found <- FALSE
levels_design <- c()

for (i in seq(1, length(terms))) {
  levels <- levels(info[,terms[i]])
  if (all(contrast_obj %in% levels)) { 
    found <- TRUE
  }
  if (i == 1){
    levels_design <- c(levels_design, levels)
  } else {
    levels_design <- c(levels_design, levels[-1])
  }
}
if (! found) { stop("The contrast you required cannot be used: please use the correct contrast levels") }

colnames(design) <- levels_design

### https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
### prova provata (te lo dico perchè alla fine è la cosa più rapida) puoi prendere i dati di una probe molto signifcativa con logfc grande in un senso e fare il boxplot dei suoi valori confrontando DNMT3B high vs low, che sto tipo di controlli son meglio di qualsiasi mia memoria e/o url
fit <- lmFit(mvalues, design)
cont.matrix <- makeContrasts(contrast, levels=design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

final_table <- topTable(fit2, adjust="BH", sort.by="p", number=nrow(fit2))
colnames(final_table) <- c("log2FC", "AveExpr", "t", "PVal", "adjPValFDR", "B")

alpha <- opt$alpha
lfc <- opt$lfc

summ <- as.data.frame.matrix(summary(decideTests(fit2, method="global", adjust.method="BH", p.value=alpha)))

plot_volcano <- function(df, alpha, lfc, outfile, title) {
  df$Significance <- ifelse(abs(df$log2FC) > lfc & df$adjPValFDR < alpha, "both", ifelse(abs(df$log2FC) > lfc, "LogFC",
                                                                                                  ifelse(df$adjPValFDR < alpha, "Adjusted.Pval", "NS")))
  df$Significance <- factor(df$Significance, levels=c("LogFC", "Adjusted.Pval", "both", "NS"))
  df[df$adjPValFDR ==0,"Adjusted.Pval"] <- .Machine$double.xmin
  p <- ggplot(df, aes(log2FC, -log10(adjPValFDR))) +
    geom_point(aes(col = Significance),size=0.5) + theme_bw() +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#999999"), drop=FALSE) + # red orange green black -> orange blue green gray 
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
        # legend.key.height = unit(1, 'cm'), #change legend key height
        # legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=13), #change legend title font size
        legend.text = element_text(size=11)) + #change legend text font size 
    ggtitle(title)
  
  nsign <- nrow(df[df$Significance=="both",])
  if (nsign > 20) {
    dfUP <- df[df$log2FC>lfc,]
    dfDOWN <- df[df$log2FC<(-lfc),]
    fin <- rbind(dfUP[1:10,], dfDOWN[1:10,])
    p + geom_text_repel(data=fin, aes(label=rownames(fin)))
  } else {
    p + geom_text_repel(data=df[df$sig=="both",], aes(label=rownames(df[df$sig=="both",])))
  }
  ggsave(outfile)
}

title <- paste0("Differential expression results ", opt$what)
plot_volcano(final_table, alpha, lfc, opt$output_volcano, title)

print("writing...")

write.table(final_table, opt$output, sep='\t', quote=FALSE, row.names=TRUE)
write.table(summ, opt$output_summ, sep='\t', quote=FALSE, row.names=TRUE)

save.image(opt$image_out)
