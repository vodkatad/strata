library(ComplexHeatmap)
library(ggplot2)

op_f <- snakemake@output[['op']]
opwf_f <- snakemake@output[['opwf']]
opmute_f <- snakemake@output[['opmute']]
op_data_f <- snakemake@output[['op_data']]
binary_f <- snakemake@input[['mutmat']]
wf_f <- snakemake@input[['wf']]


mut <- read.table(binary_f, sep="\t")
wf <- read.table(wf_f, sep="\t", stringsAsFactors = FALSE)
colnames(wf) <- c('smodel', 'perc')
#mut <- as.data.frame(t(mut))

chemionobio <- setdiff(wf$smodel, rownames(mut))
#list <- as.list(chemionobio)
mut2 <- data.frame(matrix(nrow = length(chemionobio), ncol = ncol(mut)))
rownames(mut2) <- chemionobio
colnames(mut2) <- colnames(mut)
mut <- rbind(mut, mut2)
# bionochemio <- setdiff(rownames(mut), wf$smodel)
# mut$case <- rownames(mut)
# mut <- mut %>% filter(!case %in% bionochemio)
# mut$case <- NULL

mat_list <- list(HR=as.matrix(t(mut)))

col = c("wt"="#969595", "Mut"= "firebrick", "None"="#CCCCCC")
alter_fun = list(
  wt = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["wt"], col = NA))
  },
  Mut = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Mut"], col = NA))
  },
  None = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["None"], col = NA))
  }
)



# to order genes depending on total muts
s <- t(mat_list[[1]])
su <- colSums(s, na.rm = TRUE)
su <- su[order(-su)]

mat_list <-mat_list[[1]]
mat_list[is.na(mat_list)] <- 3
mat_list <- ifelse(mat_list==1, "Mut", ifelse(mat_list==0, "wt", "None"))



ordine_casi <- wf$smodel
ordered_mat_list <- mat_list[, ordine_casi]
print(length(ordine_casi))

#oncoPrint(mat_list2, alter_fun = alter_fun, col = col)
op <- oncoPrint(ordered_mat_list, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=8)), show_pct=FALSE)

#pdf(op_f, width=10, height=6, family="sans")

#print(op)
#graphics.off()

setEPS()
postscript(op_f, width=9.4, height=7, family="sans")
print(op)
graphics.off()



get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50, 4, ifelse(x > 35, 2, 8))
  #res[is.na(x)] <- 'black'
  return(as.numeric(res))
}
op2 <- oncoPrint(ordered_mat_list, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                 top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                    Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc))),
                                                    annotation_name_gp=gpar(fontsize=8), height = unit(4, "cm")), show_pct=FALSE)

# pdf(opwf_f, width=10, height=6, family="sans")
# print(op2)
# graphics.off()

setEPS()
postscript(opwf_f, width=9.4, height=7, family="sans")
print(op2)
graphics.off()

op3 <- oncoPrint(ordered_mat_list, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, show_row_names=FALSE, show_column_names=FALSE, show_pct = FALSE)#,
                 #top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                 #                                  Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc))),
                 #                                  annotation_name_gp=gpar(fontsize=0), height = unit(4, "cm")), show_pct=FALSE,
                 #heatmap_legend_param=list())

# pdf(opmute_f, width=10, height=6, family="sans")
# print(op3)
# graphics.off()

setEPS()
postscript(opmute_f, width=9.4, height=7, family="sans")
print(op3)
graphics.off()
save.image(op_data_f)
