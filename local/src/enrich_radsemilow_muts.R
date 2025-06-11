setwd('/mnt/trcanmed/snaketree/prj/strata/dataset/figures')
load('w6_hr_oncoprintndel_cn_grigi2.Rdata')

fra <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/IHC_m29_merged.tsv"
fra <- read.table(fra, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
fra$model <- rownames(fra)

sm29 <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/samples_data"
sm29 <- read.table(sm29, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sm29$batch <- NULL
sm29 <- sm29 %>% filter(!duplicated(model))

m29_fra <- merge(fra, sm29, by="model")

m29_fra_nohighrad <- m29_fra %>% filter(!model %in% c("CRC0029", "CRC0151", "CRC0479", "CRC0204"))

m29_fra_nohighrad <- m29_fra_nohighrad[,c('model', 'type')]

data <- t(mat)
data[grepl('Del', data)] <- 'Del'
m <- merge(data, m29_fra_nohighrad, by.x="row.names", by.y="model")
rownames(m) <- m$Row.names
m$Row.names <- NULL

fisher <- function(gene, mydata, test) {
  ct <- table(mydata$type, mydata[,gene] == test)
  if (ncol(ct) == 2) {
    ft <- fisher.test(ct)
    return(ft$p.value)
  } else {
    return(NA)
  }
}

fisha <- as.data.frame(sapply(colnames(data), fisher, m, 'Amplification'))
fishd <- as.data.frame(sapply(colnames(data), fisher, m, 'Del'))
fisha <- fisha[!is.na(fisha[,1]),, drop=F]
fisha <- fisha[order(fisha[,1]),, drop=F]

fishd <- fishd[!is.na(fishd[,1]),, drop=F]
fishd <- fishd[order(fishd[,1]),, drop=F]
###############
load('w6_hr_11E2_singlecopydel2.Rdata')
fra <- "/mnt/trcanmed/snaketree/prj/strata/dataset/figures/IHC_m29_merged.tsv"
fra <- read.table(fra, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
fra$model <- rownames(fra)

sm29 <- "/scratch/trcanmed/DE_RNASeq/dataset/magnifici29/samples_data"
sm29 <- read.table(sm29, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sm29$batch <- NULL
sm29 <- sm29 %>% filter(!duplicated(model))

m29_fra <- merge(fra, sm29, by="model")

m29_fra_nohighrad <- m29_fra %>% filter(!model %in% c("CRC0029", "CRC0151", "CRC0479", "CRC0204"))

m29_fra_nohighrad <- m29_fra_nohighrad[,c('model', 'type')]
data <- t(mat)
m <- merge(data, m29_fra_nohighrad, by.x="row.names", by.y="model")
rownames(m) <- m$Row.names
m$Row.names <- NULL

fishl <- as.data.frame(sapply(colnames(data), fisher, m, 'Singlecopy_loss'))
fishl <- fishl[!is.na(fishl[,1]),, drop=F]
fishl <- fishl[order(fishl[,1]),, drop=F]

