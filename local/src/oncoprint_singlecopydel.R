library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

op_f <- snakemake@output[['op']]
opwf_f <- snakemake@output[['opwf']]
op_data_f <- snakemake@output[['op_data']]
cnmat_f <- snakemake@input[['cnmat']]
wf_f <- snakemake@input[['wf']]
all_genes_f <- snakemake@input[['all_genes']]
log_f <- snakemake@log[['log']]
wil_f <- snakemake@output[['wilcox_genes']]
gistic_event <- as.numeric(snakemake@params[['gistic']])
gistic_event_name <- snakemake@params[['gistic_name']]

# Do we work on the whole cohort or have a selection? (for the suppl on zona grigia)
keep <- snakemake@params[['keep']]
keep_f <- snakemake@params[['keep']]
if (keep != 'all' && !file.exists(keep_f)) {
  stop('The keep params must be an existent file with the wanted models if you use it != all')
}

# wf and complete list of hr genes  
wf <- read.table(wf_f, sep="\t", stringsAsFactors = FALSE)
colnames(wf) <- c('smodel', 'perc')
genes <- read.table(all_genes_f, sep="\t", header=TRUE)

# filter wf if we have a keep
if (keep_f != 'all') {
  keep_samples <- read.table(keep_f, sep="\t", header=F, stringsAsFactors = FALSE)
  wf <- wf[wf$smodel %in% keep_samples$V1,]
}

save.image('stress.Rdata')
cn <- t(read.table(cnmat_f, sep="\t", header=TRUE, row.names=1))
rownames(cn) <- substr(rownames(cn), 0, 7)
cnkeep <- intersect(rownames(cn), wf$smodel)
cn <- cn[rownames(cn) %in% cnkeep,]
wf <- wf[wf$smodel %in% cnkeep,]
cn <- cn[wf$smodel,]

cnhigh <- ifelse(cn== gistic_event, gistic_event_name, 'background')
#cn <- t(apply(mut, 2, as.character))
cn <- t(cnhigh)

mat <- cn

# add genes without any event
genes <- read.table(all_genes_f, sep="\t", header=TRUE)
to_add <- setdiff(genes$gene_symbol, rownames(mat))

tot_genes <- nrow(genes)
tot_addgenes <- length(to_add)
tot_models <- ncol(mat)
to_add_mat <- matrix(rep('background', tot_addgenes*tot_models), nrow=tot_addgenes, ncol=tot_models)
rownames(to_add_mat) <- to_add
colnames(to_add_mat) <- colnames(mat)
mat <- rbind(mat, to_add_mat)

sink(log_f)
print('Missing in CN:')
print(missing)
sink()

info <- which(mat!="background", arr.ind=T)
models <- unique(rownames(info))
altgenes <- unique(colnames(mat)[info[,2]])
sink(log_f)
#'single-copy losses (GISTIC thresholded value = -1) were identified for XX genes in XX models.''
print(paste0('for ', length(altgenes), ' genes'))
print(paste0('in ', length(models), ' models'))
print(altgenes)
sink()


sink(log_f, append=T)
print(dim(mat))
sink()

pal=brewer.pal(9,'YlOrRd')[seq(3,9)]
#pal <- rev(hcl.colors(7, palette = "Sunset"))
col = c('Singlecopy_loss'='blue')

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Singlecopy_loss = function(x, y, w, h) { ## eval here TODO
    grid.rect(x, y, w*0.2, h, 
              gp = gpar(fill = col['Singlecopy_loss'], col = NA))
  }
)

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50,  rgb(red=30, green=85, blue=130, maxColorValue=255), ifelse(x > 35,  rgb(red=165, green=0, blue=25, maxColorValue=255), 8))
  #res[is.na(x)] <- 'black'
  return(res)
}

s <- t(mat!='background')
su <- colSums(s)
names(su) <- rownames(mat)
su <- su[order(-su)]

save.image('stress.Rdata')
op <- oncoPrint(mat, alter_fun = alter_fun, column_order=wf$smodel, col=col, row_order = names(su),
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                 top_annotation = HeatmapAnnotation(Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc), col=NA), height= unit(5, 'cm'), border=FALSE),
                                                    cbar = anno_oncoprint_barplot()), show_pct=FALSE)


pdf(file=opwf_f, width=9.4, height=7, family="sans")
print(op)
graphics.off()

op <- oncoPrint(mat, alter_fun = alter_fun, column_order=wf$smodel, col=col, row_order = names(su),
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=8)), show_pct=FALSE)


pdf(file=op_f, width=9.4, height=7, family="sans")
print(op)
graphics.off()

### wilcoxon
wilcoxon_perc_3 <- function(gene, mydata) {
  percwt <- mydata[mydata[,gene] < thr, 'perc']
  percmut <- mydata[mydata[,gene]  > thr, 'perc']
  if (length(percwt) != 0 && length(percmut) != 0) {
    wt <- wilcox.test(percwt, percmut)
    return(c(wt$p.value, length(percwt), length(percmut)))
  } else {
    return(c(NA, NA, NA))
  }
}

#wil <- as.data.frame(t(sapply(colnames(num), wilcoxon_perc_3, data_merged)))
#colnames(wil) <- c('pvalue', 'nwt', 'nmut')
#wil$padj <- p.adjust(wil$pvalue, method="BH")
write.table(data.frame(TODO='wow'), file=wil_f, quote=FALSE, sep="\t")

save.image(op_data_f)
