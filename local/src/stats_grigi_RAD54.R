setwd('/mnt/cold1/snaketree/prj/strata/dataset/figures/')
load('w3_hr_oncoprintndel_cn.Rdata')

m29 <- read.table('m29')$V1
grigi <- read.table('m29_noRAD')$V1
nongrigi <- setdiff(m29, grigi)

wholedata <- mat

# RAD54B amplificato nei grigi non grigi
mat['RAD54B',] == "Amplification"
colnames(mat)[mat['RAD54B',] == "Amplification"]
rad54b_ampl <- colnames(mat)[mat['RAD54B',] == "Amplification"]
tabrad <- wf[wf$smodel %in% rad54b_ampl,]

tabrad$m29 <- ifelse(tabrad$smodel %in% m29, 'm29', 'other')
tabrad$grigi <- ifelse(tabrad$smodel %in% m29 & tabrad$smodel %in% grigi, 'grigi', 'non_grigi')
tabrad[tabrad$m29=="other" & tabrad$m29=="non_grigi", 'grigi'] <- ''

# fisher R-S nei grigi singlecopyloss
load('w3_hr_11E2_singlecopydel.Rdata')
labels_grigi <- read.table('/mnt/trcanmed/snaketree/prj/strata/local/share/data/grigi_annot', sep="\t")
colnames(labels_grigi) <- c('smodel', 'cet')

singlecn <- t(mat)
m <- merge(singlecn, labels_grigi, by.x="row.names", by.y="smodel")
row.names(m) <- m$Row.names
m$Row.names <- NULL

fisher_two <- function(gene, mydata) {
  fid <- table(mydata[, gene], mydata[, 'cet'])
  if (all(dim(fid) == c(2, 2))) {
    fi <- fisher.test(fid)
    return(fi$p.value)
  } else {
    return(NA)
  }
}

wil <- as.data.frame(sapply(colnames(m), fisher_two, m))
colnames(wil) <- c('pvalue')
wil$padj <- p.adjust(wil$pvalue, method="BH")
wil <- wil[order(wil$pvalue),, drop=F]

# fisher grigi-m29 per muts

universe <- wholedata[, colnames(wholedata) %in% m29]
universe[grepl('Del', universe)] <- 'SNV'
universe[!grepl('SNV', universe)] <- 'background'

uni <- as.data.frame(t(universe))
uni$grigi <- ifelse(rownames(uni) %in% grigi, 'grigi', 'non_grigi')

fisher_tremendi <- function(gene, mydata) {
  fid <- table(mydata[, gene], mydata[, 'grigi'])
  if (all(dim(fid) == c(2, 2))) {
    fi <- fisher.test(fid)
    return(fi$p.value)
  } else {
    return(NA)
  }
}

wil <- as.data.frame(sapply(colnames(uni), fisher_tremendi, uni))
colnames(wil) <- c('pvalue')
wil$padj <- p.adjust(wil$pvalue, method="BH")
wil <- wil[order(wil$pvalue),, drop=F]
