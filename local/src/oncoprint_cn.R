library(ComplexHeatmap)
library(ggplot2)

op_f <- snakemake@output[['op']]
opwf_f <- snakemake@output[['opwf']]
opmute_f <- snakemake@output[['opmute']]
op_data_f <- snakemake@output[['op_data']]
binary_f <- snakemake@input[['mutmat']]
wf_f <- snakemake@input[['wf']]

save.image('pippo.Rdata')
mut <- t(read.table(binary_f, sep="\t", header=TRUE, row.names=1))
rownames(mut) <- substr(rownames(mut), 0, 7)
wf <- read.table(wf_f, sep="\t", stringsAsFactors = FALSE)
colnames(wf) <- c('smodel', 'perc')
keep <- intersect(rownames(mut), wf$smodel)
mut <- mut[rownames(mut) %in% keep,]
wf <- wf[wf$smodel %in% keep,]
mut <- mut[wf$smodel,]

cn = t(apply(mut, 2, as.character))
colnames(cn) <- rownames(mut)
col = c("0"="#cccccc", "1"= "#e67272",  "2"="#ff0000", "-1"="#6e72e9", "-2"= "#0000ff")

alter_fun = list(
  '0' = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col['0'], col = NA))
  },
  '1' = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col['1'], col = NA))
  },
  '2' = function(x, y, w, h) {
    grid.rect(x, y, w*1, h, 
              gp = gpar(fill = col['2'], col = NA))
  },
  '-1' = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col['-1'], col = NA))
  },
  '-2' = function(x, y, w, h) {
    grid.rect(x, y, w*1, h, 
              gp = gpar(fill = col['-2'], col = NA))
  }
)

# to order genes depending on total muts
s <- cn
su <- colSums(mut != 0)
su <- su[order(-su)]

op <- oncoPrint(cn, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=8)), show_pct=FALSE)

pdf(op_f, width=10, height=6, family="sans")

print(op)
graphics.off()

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50, 4, ifelse(x > 35, 2, 8))
  #res[is.na(x)] <- 'black'
  return(as.numeric(res))
}
op2 <- oncoPrint(cn, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                  Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc))),
                                                  annotation_name_gp=gpar(fontsize=8), height = unit(4, "cm")), show_pct=FALSE)

pdf(opwf_f, width=10, height=6, family="sans")
print(op2)
graphics.off()

op3 <- oncoPrint(cn, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, show_row_names=FALSE, show_column_names=FALSE,
                 top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                    Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc))),
                                                    annotation_name_gp=gpar(fontsize=0), height = unit(4, "cm")), show_pct=FALSE,
                 heatmap_legend_param=list())

pdf(opmute_f, width=10, height=6, family="sans")
print(op3)
graphics.off()

save.image(op_data_f)
