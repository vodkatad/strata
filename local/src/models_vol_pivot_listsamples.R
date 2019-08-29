library(ggplot2, quietly=TRUE)
library(WriteXLS, quietly=TRUE)

n_mutated <- snakemake@params[["cutoff"]]
volData <- snakemake@input[["volumes"]]
pivotData <- snakemake@input[["pivot"]]
msi <- snakemake@input[["msi"]]
allData <- snakemake@input[["all_samples"]]
debug <- snakemake@params[["debug"]]
linearOut <- snakemake@output[["linear"]]
logisticOut <- snakemake@output[["logistic"]]
linearOutXls <- snakemake@output[["linearxls"]]
plotdir <- snakemake@params[["plotdir"]]

dir.create(plotdir, showWarnings = FALSE)

mutdata <- read.table(pivotData, sep="\t", header=TRUE)
averages <- read.table(volData, sep="\t", header=TRUE)
all <- read.table(allData, sep="\t", header=FALSE)
rownames(averages) <- averages$case
rownames(mutdata) <- mutdata[,1]
mutdata[,1] <- NULL

missing <- setdiff(all$V1, rownames(mutdata))
wanted <- length(missing)*ncol(mutdata)
nomut <- matrix(rep(0, wanted), ncol=ncol(mutdata))
nomut_df <- as.data.frame(nomut)
rownames(nomut_df) <- missing
colnames(nomut_df) <- colnames(mutdata)
# we add as completely not mutated the list of studied case with no reported muts
mutdata <- rbind(mutdata, nomut_df)

feasible <- rownames(mutdata)[rownames(mutdata) %in% averages$case]
mutdatafeas <- mutdata[rownames(mutdata) %in% feasible,, drop=FALSE]

nmf <- colSums(mutdatafeas)
#summary(mm)
wanted <- names(nmf[nmf>=n_mutated])
mutdatafeas <- mutdatafeas[,colnames(mutdatafeas) %in% wanted, drop=FALSE]

# Fix  TODO output real n. of MSI with mutations!
msi <- read.table(msi, sep="\t", header=FALSE)
averages[,"msi"] <- rep(0, nrow(averages))
averages[rownames(averages) %in% msi[,1],"msi"] <- 1

av <- averages[rownames(averages) %in% feasible,, drop=FALSE]
cat("MSI:\t")
cat(nrow(av[av$msi==1,]))
cat("\n")

if (debug == "yes") {
  save.image(file=paste0(linearOut,'.debug','.RData'))
}

onemodel <- function(mut, avg, samples) {
  muts <- data.frame(status=mut, row.names=samples)
  m <- merge(avg, muts, by="row.names")
  model <- lm("perc~status+msi",data=m)
  sm <- summary(model)
  pvF <- NA
  pvT <- NA
  if (length(sm$fstatistic)==3) {
    f <- sm$fstatistic[1]
    df1 <- sm$fstatistic[2]
    df2 <- sm$fstatistic[3]
    pvF <- pf(f,df1,df2,low=F)
    coef <- sm$coefficients[,4]
    if (any(names(coef) == "status")) {
        pvT <- coef[names(coef) == "status"]
    }
  }
  medianmut <- median(m[m$status==1, "perc"])
  medianwt <- median(m[m$status==0, "perc"])
  
  return(c(sm$r.squared, pvF, pvT, medianmut, medianwt))
}

res <- as.data.frame(t(as.data.frame(apply(mutdatafeas, 2, onemodel, averages, rownames(mutdatafeas)))))
colnames(res) <- c("R2","pval_global", "pval_gene", "medianmut", "medianwt")
res$adj_pval_global <- p.adjust(res$pval_global, method="BH")
res$adj_pval_gene <- p.adjust(res$pval_gene, method="BH")
res <- res[order(res$pval_gene),]
write.table(res, file=linearOut, sep="\t", quote=FALSE)

#mres <- merge(res, annot, by.x="row.names", by.y="gs")
mres <- res
#if (snakemake@wildcards[[1]] != "genepathways") {
if (!grepl("pathways" ,snakemake@wildcards[[1]])) {
    library ("AnnotationDbi")
    library("org.Hs.eg.db") # very generic indeed! Bravo Elena!
    if (grepl("cnv",snakemake@wildcards[[1]])) {
        gs <- sapply(strsplit(rownames(mres),"_"), function(x) {x[1]})
        mres$desc <- mapIds(org.Hs.eg.db,keys=gs, column="GENENAME", keytype="SYMBOL", multiVals="first") # first chooses the first value, for description 
    } else {
        mres$desc <- mapIds(org.Hs.eg.db,keys=row.names(mres), column="GENENAME", keytype="SYMBOL", multiVals="first") # first chooses the first value, for description 
    }
    # should be good enough.
    # we need to put the last column of mres as second one, the number of columns is fixed but I'm lazy to check by hand
    last <- ncol(mres)
    others <- seq(1, last-1)
    mres <- mres[,c(last, others)]
}
WriteXLS(mres, ExcelFileName = linearOutXls, row.names = TRUE, col.names = TRUE) 

plots <- function(wanted, mut, avg) {
    idx <- which(colnames(mut)==wanted)
    data <- data.frame(mut=mut[,idx], row.names = rownames(mut))
    m <- merge(avg, data, by="row.names")
    m$mut <- as.factor(m$mut)
    ggplot(m, aes(x = mut, y=perc, fill=mut)) + geom_boxplot() +geom_jitter()+scale_fill_manual(values=c("grey","red"))+theme_bw()+theme(text = element_text(size=15))+ylab("Delta %volume")+xlab(paste0("Mutated ", wanted))
    ggsave(paste0(plotdir, "/", wanted, "_boxplot.png"))
    #table(m$class, m$mut)
    ggplot(m, aes(y=perc,x=reorder(case, -perc),fill=mut))+geom_col()+ylab("Delta %volume")+xlab("Case")+theme_bw()+theme(axis.text.x = element_text(size=15, angle = 90, hjust = 1, vjust=0.5))+scale_fill_manual(values=c("grey","red"))
    ggsave(paste0(plotdir, "/", wanted, "_waterfall.png"))
}

wanted <- head(rownames(res), n=20)
garbage <- lapply(wanted, plots, mutdatafeas, averages)

averages$class <- ifelse(averages$perc < -50, 'OR', ifelse(averages$perc > 35, 'PD', 'SD'))
averages$class <- as.factor(averages$class)
classmodel <- function(mut, avg, samples) {
  muts <- data.frame(status=mut, row.names=samples)
  m <- merge(avg, muts, by="row.names")
  model <- glm("class~status", data=m, family=binomial(link="logit"))
  sm <- summary(model)
  pv <- 1 - pchisq(sm$null.deviance-sm$deviance, df=sm$df.null-sm$df.residual)
  coef <- sm$coefficients[,4]
  pvT <- NA
  if (any(names(coef) == "status")) {
    pvT <- coef[names(coef) == "status"]
  }
  return(c(sm$aic, pv, pvT))
}
res2 <- as.data.frame(t(as.data.frame(apply(mutdatafeas, 2, classmodel, averages, rownames(mutdatafeas)))))
#head(res2[order(res2[,2]),])
colnames(res2) <- c("AIC","pval_global","pval_gene")
res2$adj_pval_global <- p.adjust(res2$pval_global, method="BH")
res2$adj_pval_gene <- p.adjust(res2$pval_gene, method="BH")
write.table(res2, file=logisticOut, sep="\t", quote=FALSE)

if (debug == "yes") {
  save.image(file=paste0(linearOut,'.debug','.RData'))
}

