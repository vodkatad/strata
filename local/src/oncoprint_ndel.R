library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

op_f <- snakemake@output[['op']]
opwf_f <- snakemake@output[['opwf']]
opmute_f <- snakemake@output[['opmute']]
op_data_f <- snakemake@output[['op_data']]
binary_f <- snakemake@input[['mutmat']]
wf_f <- snakemake@input[['wf']]
all_genes_f <- snakemake@input[['all_genes']]
log_f <- snakemake@log[['log']]

mut <- read.table(binary_f, sep="\t")
wf <- read.table(wf_f, sep="\t", stringsAsFactors = FALSE)
colnames(wf) <- c('smodel', 'perc')

keep <- intersect(rownames(mut), wf$smodel)
mut <- mut[rownames(mut) %in% keep,]
wf <- wf[wf$smodel %in% keep,]
mut <- mut[wf$smodel,]

mat <- t(mut)
# to order genes depending on total muts
s <- t(mat[[1]])
su <- colSums(s)
su <- su[order(-su)]

mat <- ifelse(mat==1, 'Del1', ifelse(mat==2 ,'Del2', ifelse(mat==3, 'Del3', ifelse(mat==4, 'Del4', ifelse(mat==0, 'background', ifelse(mat==0.5, 'Del0', 'Del5'))))))


# add genes without any non syn mut
genes <- read.table(all_genes_f, sep="\t", header=TRUE)
to_add <- setdiff(genes$gene_symbol, rownames(mat))

tot_genes <- nrow(genes)
tot_addgenes <- length(to_add)
tot_models <- ncol(mat)
to_add_mat <- matrix(rep('background', tot_addgenes*tot_models), nrow=tot_addgenes, ncol=tot_models)
rownames(to_add_mat) <- to_add
colnames(to_add_mat) <- colnames(mat)
mat <- rbind(mat, to_add_mat)

## TODO numbers Livio
sink(log_f)
#In total, XX models (XX%) harbored HR gene mutations, with most cases involving mutations in more than one gene. 
#These mutations were monoallelic in XX cases (XX%) and biallelic in the remaining XX cases (XX%), 
#with XX of them (XX%) consistently predicted as deleterious variants by independent algorithms
print('Mut are TRUE:')
nn <- apply(mut, 1, function(x){any(x!=0)})
print(table(nn))
print('Tot models:')
print(tot_models)
print('frac:')
print(sum(nn)/tot_models)
sink()

pal=brewer.pal(9,'YlOrRd')[seq(3,9)]
col = c("Del0"= pal[1], "Del1"= pal[2], "Del2"= pal[3], "Del3"= pal[4], "Del4"= pal[5], "Del5"=pal[6])

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Del0 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Del0"], col = NA))
  },
  Del1 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Del1"], col = NA))
  },
  Del2 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Del2"], col = NA))
  },
  Del3 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Del3"], col = NA))
  },
  Del4 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Del4"], col = NA))
  },  
  Del5 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Del5"], col = NA))
  }
)


save.image('pippo.Rdata')
#oncoPrint(mat_list2, alter_fun = alter_fun, col = col)
op <- oncoPrint(mat, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=8)), show_pct=FALSE)
#                height=unit(129.5, 'mm'), width=unit(182.46, unit='mm'))

#pdf(op_f, width=10, height=7, family="sans")
#pdf(op_f, family="sans")

#print(op)
#graphics.off()

setEPS()
postscript(op_f, width=10.05, height=7, family="sans")
print(op)
graphics.off()

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50,  rgb(red=30, green=85, blue=130, maxColorValue=255), ifelse(x > 35,  rgb(red=165, green=0, blue=25, maxColorValue=255), 8))
  #res[is.na(x)] <- 'black'
  return(res)
}
op2 <- oncoPrint(mat, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc), col=NA), height= unit(5, 'cm'), border=FALSE),
                                                                      cbar = anno_oncoprint_barplot()), show_pct=FALSE)#, height=unit(129.5, 'mm'), width=unit(182.46, unit='mm'))

#pdf(opwf_f, width=10, height=7, family="sans")
#pdf(opwf_f, family="sans")
#print(op2)
#graphics.off()
setEPS()
postscript(opwf_f, width=10.05, height=7, family="sans")
print(op2)
graphics.off()


op3 <- oncoPrint(mat, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, show_row_names=FALSE, show_column_names=FALSE,
                 top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                    annotation_name_gp=gpar(fontsize=0), height = unit(4, "cm")), show_pct=FALSE,
                 heatmap_legend_param=list())#,  height=unit(129.5, 'mm'), width=unit(182.46, unit='mm'))

#pdf(opmute_f, width=10, height=6, family="sans")
#pdf(opmute_f, family="sans")
#print(op3)
#graphics.off()
setEPS()
postscript(opmute_f, width=10.05, height=7, family="sans")
print(op3)
graphics.off()


save.image(op_data_f)
