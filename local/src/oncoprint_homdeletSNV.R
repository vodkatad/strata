library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

op_f <- snakemake@output[['op']]
opwf_f <- snakemake@output[['opwf']]
op_data_f <- snakemake@output[['op_data']]
vaf_f <- snakemake@input[['vaf']]
mutmat_f <- snakemake@input[['mutmat']]
wf_f <- snakemake@input[['wf']]
all_genes_f <- snakemake@input[['all_genes']]
log_f <- snakemake@log[['log']]
wil_f <- snakemake@output[['wilcox_genes']]
thr <- as.numeric(snakemake@params[['thr']])

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

### loading and setup of VAF data, threshold and adding missing genes (without muts)
data <- read.table(vaf_f, sep="\t", header=TRUE, stringsAsFactors = F)

# get the max vaf when there is more than one mutation on the same gene
num <- data.frame(matrix(0, nrow=nrow(data), ncol(data)), stringsAsFactors = F)
rownames(num) <- rownames(data)
colnames(num) <- colnames(data)
for (i in seq(1, nrow(data))) {
  for (j in seq(1, ncol(data))) {
    if (!is.na(data[i,j]) && grepl(',', data[i,j])) {
      vals <- strsplit(data[i,j], split=',')
      print(max(as.numeric(unlist(vals))))
      num[i,j] <- max(as.numeric(unlist(vals)))
    } else if (!is.na(data[i,j]) && data[i,j] != '') {
      num[i,j] <- as.numeric(data[i,j])
    } else {
      num[i,j] <- 0
    }
  }
}

data_merged <- merge(num, wf, by.x="row.names", by.y='smodel')

if (keep != "all") {
  num <- num[rownames(num) %in% data_merged$Row.names,]
}
### wilcoxon
wilcoxon_perc_3 <- function(gene, mydata) {
  percwt <- mydata[mydata[,gene] < thr, 'perc']
  percmut <- mydata[mydata[,gene]  > thr, 'perc']
  if (length(percwt) != 0 && length(percmut) != 0) {
    print(gene)
    wt <- wilcox.test(percwt, percmut)
    return(c(wt$p.value, length(percwt), length(percmut)))
  } else {
    return(c(NA, NA, NA))
  }
}

wil <- as.data.frame(t(sapply(colnames(num), wilcoxon_perc_3, data_merged)))
colnames(wil) <- c('pvalue', 'nwt', 'nmut')
wil$padj <- p.adjust(wil$pvalue, method="BH")
write.table(wil, file=wil_f, quote=FALSE, sep="\t")


# add genes without any non syn mut
mat <- t(num)
to_add <- setdiff(genes$gene_symbol, rownames(mat))

tot_genes <- nrow(genes)
tot_addgenes <- length(to_add)
tot_models <- ncol(mat)
to_add_mat <- matrix(rep(0, tot_addgenes*tot_models), nrow=tot_addgenes, ncol=tot_models)
rownames(to_add_mat) <- to_add
colnames(to_add_mat) <- colnames(mat)
mat <- rbind(mat, to_add_mat)


mat3 <- data.frame(matrix(0, nrow=nrow(mat), ncol(mat)), stringsAsFactors = F)
rownames(mat3) <- rownames(mat)
colnames(mat3) <- colnames(mat)
for (i in seq(1, nrow(mat))) {
  for (j in seq(1, ncol(mat))) {
    if (mat[i,j] > thr) {
      mat3[i,j] <- 'HomozigMut'
    } else {
      mat3[i,j] <- 'Background'
    }
  }
}


mat_vaf <- t(mat3)

#### loading of nDel data
mut <- read.table(mutmat_f, sep="\t")
save.image('stress.Rdata')

keep <- intersect(rownames(mut), wf$smodel)
mut <- mut[rownames(mut) %in% keep,]
wf <- wf[wf$smodel %in% keep,]
mut <- mut[wf$smodel,]
mat <- t(mut)

mat_ndel <- ifelse(mat==1, 'Del1', ifelse(mat==2 ,'Del2', ifelse(mat==3, 'Del3', ifelse(mat==4, 'Del4', ifelse(mat==0, 'background', ifelse(mat==0.5, 'Del0', 'Del5'))))))

# add missing genes
to_add <- setdiff(genes$gene_symbol, rownames(mat_ndel))
tot_genes <- nrow(genes)
tot_addgenes <- length(to_add)
tot_models <- ncol(mat)
to_add_mat <- matrix(rep('background', tot_addgenes*tot_models), nrow=tot_addgenes, ncol=tot_models)
rownames(to_add_mat) <- to_add
colnames(to_add_mat) <- colnames(mat)
mat_ndel <- rbind(mat_ndel, to_add_mat)

# We need to get data from mat_ndel only where we have a HomozigMut in mat_vaf

# align the two matrixes
mat_vaf <- t(mat_vaf)
stopifnot(all(dim(mat_vaf)==dim(mat_ndel)))

mat_vaf <- mat_vaf[rownames(mat_ndel), colnames(mat_ndel)]

mat_ndel_hetSNV <- matrix(rep('background', nrow(mat_ndel)*ncol(mat_ndel)), nrow=nrow(mat_ndel), ncol=ncol(mat_ndel))
mat_ndel_hetSNV[which(mat_vaf == "HomozigMut")] <- mat_ndel[which(mat_vaf == "HomozigMut")]
rownames(mat_ndel_hetSNV) <- rownames(mat_ndel)
colnames(mat_ndel_hetSNV) <- colnames(mat_ndel)
# to order genes depending on total shown alterations
s <- t(mat_ndel_hetSNV!='background')
su <- colSums(s)
names(su) <- rownames(mat_ndel_hetSNV)
su <- su[order(-su)]


pal=brewer.pal(9,'YlOrRd')[seq(3,9)]
#pal <- rev(hcl.colors(7, palette = "Sunset"))
col = c("Del0"= pal[1], "Del1"= pal[2], "Del2"= pal[3], "Del3"= pal[4], "Del4"= pal[5], "Del5"=pal[6])

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Del0 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, 
              gp = gpar(fill = col["Del0"], col = NA))
  },
  Del1 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, 
              gp = gpar(fill = col["Del1"], col = NA))
  },
  Del2 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, 
              gp = gpar(fill = col["Del2"], col = NA))
  },
  Del3 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, 
              gp = gpar(fill = col["Del3"], col = NA))
  },
  Del4 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, 
              gp = gpar(fill = col["Del4"], col = NA))
  },  
  Del5 = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, 
              gp = gpar(fill = col["Del5"], col = NA))
  }
)

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50,  rgb(red=30, green=85, blue=130, maxColorValue=255), ifelse(x > 35,  rgb(red=165, green=0, blue=25, maxColorValue=255), 8))
  #res[is.na(x)] <- 'black'
  return(res)
}

sink(log_f)
print(dim(mat_ndel_hetSNV))
sink()


op <- oncoPrint(mat_ndel_hetSNV, alter_fun = alter_fun, column_order=wf$smodel, col=col, row_order = names(su),
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                 top_annotation = HeatmapAnnotation(Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc), col=NA), height= unit(5, 'cm'), border=FALSE),
                                                    cbar = anno_oncoprint_barplot()), show_pct=FALSE)


pdf(file=opwf_f, width=9.4, height=7, family="sans")
print(op)
graphics.off()

op <- oncoPrint(mat_ndel_hetSNV, alter_fun = alter_fun, column_order=wf$smodel, col=col, row_order = names(su),
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=8)), show_pct=FALSE)


pdf(file=op_f, width=9.4, height=7, family="sans")
print(op)
graphics.off()

save.image(op_data_f)
