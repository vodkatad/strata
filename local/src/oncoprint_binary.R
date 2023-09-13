library(ComplexHeatmap)
library(ggplot2)

op_f <- snakemake@output[['op']]
opwf_f <- snakemake@output[['opwf']]
opmute_f <- snakemake@output[['opmute']]
op_data_f <- snakemake@output[['op_data']]
binary_f <- snakemake@input[['mutmat']]
wf_f <- snakemake@input[['wf']]

save.image('pippo.Rdata')
mut <- read.table(binary_f, sep="\t")
wf <- read.table(wf_f, sep="\t", stringsAsFactors = FALSE)
colnames(wf) <- c('smodel', 'perc')

keep <- intersect(rownames(mut), wf$smodel)
mut <- mut[rownames(mut) %in% keep,]
wf <- wf[wf$smodel %in% keep,]
mut <- mut[wf$smodel,]

mat_list <- list(HR=as.matrix(t(mut)))

col = c("HR"= "firebrick")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HR = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["HR"], col = NA))
  }
)

# to order genes depending on total muts
s <- t(mat_list[[1]])
su <- colSums(s)
su <- su[order(-su)]

#oncoPrint(mat_list2, alter_fun = alter_fun, col = col)
op <- oncoPrint(mat_list, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
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
op2 <- oncoPrint(mat_list, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                  Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc))),
                                                  annotation_name_gp=gpar(fontsize=8), height = unit(4, "cm")), show_pct=FALSE)

pdf(opwf_f, width=10, height=6, family="sans")
print(op2)
graphics.off()

op3 <- oncoPrint(mat_list, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, show_row_names=FALSE, show_column_names=FALSE,
                 top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                    Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc))),
                                                    annotation_name_gp=gpar(fontsize=0), height = unit(4, "cm")), show_pct=FALSE,
                 heatmap_legend_param=list())

pdf(opmute_f, width=10, height=6, family="sans")
print(op3)
graphics.off()

save.image(op_data_f)
