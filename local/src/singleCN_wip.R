library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
load('/mnt/trcanmed/snaketree/prj/strata/dataset/figures/noia.Rdata')

cn <- cn[wf$smodel,]
cnhigh <- ifelse(cn== 2, 'Gain', ifelse(cn==-2, 'HomDel', ifelse(cn==-1, 'Del', 'WT')))
#cn <- t(apply(mut, 2, as.character))
cn <- t(cnhigh)

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

## add genes without cn data and order in the same way
missing <- setdiff(rownames(mat), rownames(cn))
tot_addgenes <- length(missing)
tot_models <- ncol(cn)
to_add_mat <- matrix(rep('WT', tot_addgenes*tot_models), nrow=tot_addgenes, ncol=tot_models)
rownames(to_add_mat) <- missing
colnames(to_add_mat) <- colnames(cn)
cn <- rbind(cn, to_add_mat)
cn <- cn[rownames(mat),]
stopifnot(all(rownames(mat)==rownames(cn)))
stopifnot(all(colnames(mat)==colnames(cn)))

mat[cn=="Gain"] <- "Amplification" # https://docs.cbioportal.org/user-guide/faq/#what-do-amplification-gain-deep-deletion-shallow-deletion-and--2--1-0-1-and-2-mean-in-the-copy-number-data
mat[cn=="HomDel"] <- "HomDel"
# Single dels are instead together with some mutations
# unique(mat[cn=="Del"])

mat2 <- mat
w <- F
for (i in seq(1, nrow(mat))) {
  for (j in seq(1, ncol(mat))) {
    if (cn[i,j] == 'Del') {
      if (mat[i,j] != 'background') {
        w <- T
      }
      mat2[i,j] <- paste0(mat2[i,j], ',Del')
    }
  }
  if (w) { 
    print(rownames(mat)[i])
  }
  w <- F
}

pal=brewer.pal(9,'YlOrRd')[seq(3,9)]
#pal <- rev(hcl.colors(7, palette = "Sunset"))
col = c("Del0"= pal[1], "Del1"= pal[2], "Del2"= pal[3], "Del3"= pal[4], "Del4"= pal[5], "Del5"=pal[6], 'Amplification'='mediumorchid1', 'HomDel'='aquamarine', 'Del'='blue')

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
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Amplification"], col = NA))
  },  
  HomDel = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["HomDel"], col = NA))
  },
  Del = function(x, y, w, h) {
    grid.rect(x, y, w*0.2, h, 
              gp = gpar(fill = col["Del"], col = NA))
  }
)

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50,  rgb(red=30, green=85, blue=130, maxColorValue=255), ifelse(x > 35,  rgb(red=165, green=0, blue=25, maxColorValue=255), 8))
  #res[is.na(x)] <- 'black'
  return(res)
}
op2 <- oncoPrint(mat2, alter_fun = alter_fun, col = col, row_order = names(su), column_order=wf$smodel,
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                 top_annotation = HeatmapAnnotation(Irinotecan = anno_barplot(wf$perc,  gp = gpar(fill = get_recist(wf$perc), col=NA), height= unit(5, 'cm'), border=FALSE),
                                                    cbar = anno_oncoprint_barplot()), show_pct=FALSE)#, height=unit(129.5, 'mm'), width=unit(182.46, unit='mm'))

## 3/5 deleterious are somehow relevant?

data_merged <- merge(mut, wf, by.x="row.names", by.y='smodel')
stopifnot(nrow(mut)==nrow(data_merged))

wilcoxon_perc <- function(gene, mydata) {
  percwt <- mydata[mydata[,gene] == 0, 'perc']
  percmut <- mydata[mydata[,gene] >= 3, 'perc']
  if (length(percwt) != 0 && length(percmut) != 0) {
    wt <- wilcox.test(percwt, percmut)
    return(c(wt$p.value, length(percwt), length(percmut)))
  } else {
    return(c(NA, NA, NA))
  }
}

wil <- as.data.frame(t(sapply(colnames(mut), wilcoxon_perc, data_merged)))
colnames(wil) <- c('pvalue', 'nwt', 'nmut')
wil$padj <- p.adjust(wil$pvalue, method="BH")
wil <- wil[order(wil$pvalue),]

