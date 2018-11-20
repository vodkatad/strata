#library(ggplot2, quietly=TRUE)

n_mutated <- snakemake@params[["cutoff"]]
volData <- snakemake@input[["volumes"]]
pivotData <- snakemake@input[["pivot"]]
debug <- snakemake@params[["debug"]]
linearOut <- snakemake@output[["linear"]]
logisticOut <- snakemake@output[["logistic"]]

mutdata <- read.table(pivotData, sep="\t", header=TRUE)
averages <- read.table(volData, sep="\t", header=TRUE)
# rownames(averages) <- averages$case
rownames(mutdata) <- mutdata[,1]
mutdata[,1] <- NULL

feasible <- rownames(mutdata)[rownames(mutdata) %in% averages$case]
mutdatafeas <- mutdata[rownames(mutdata) %in% feasible,]

nmf <- colSums(mutdatafeas)
#summary(mm)
wanted <- names(nmf[nmf>=n_mutated])
mutdatafeas <- mutdatafeas[,colnames(mutdatafeas) %in% wanted]


if (debug == "yes") {
  save.image(file=paste0(linearOut,'.debug','.RData'))
}

onemodel <- function(mut, avg, samples) {
  muts <- data.frame(status=mut, row.names=samples)
  m <- merge(avg, muts, by="row.names")
  model <- lm("perc~status",data=m)
  sm <- summary(model)
  pv <- NA
  if (length(sm$fstatistic)==3) {
    f <- sm$fstatistic[1]
    df1 <- sm$fstatistic[2]
    df2 <- sm$fstatistic[3]
    pv <- pf(f,df1,df2,low=F)
  }
  return(c(sm$r.squared, pv))
}

res <- as.data.frame(t(as.data.frame(apply(mutdatafeas, 2, onemodel, averages, rownames(mutdatafeas)))))
colnames(res) <- c("R2","pval")
res$adj_pval <- p.adjust(res$pval, method="bonferroni")
#head(res[order(res[,2]),])
write.table(res, file=linearOut, sep="\t", quote=FALSE)

averages$class <- ifelse(averages$perc < -50, 'OR', ifelse(averages$perc > 35, 'PD', 'SD'))
averages$class <- as.factor(averages$class)
classmodel <- function(mut, avg, samples) {
  muts <- data.frame(status=mut, row.names=samples)
  m <- merge(avg, muts, by="row.names")
  model <- glm("class~status", data=m, family=binomial(link="logit"))
  sm <- summary(model)
  pv <- 1 - pchisq(sm$null.deviance-sm$deviance, df=sm$df.null-sm$df.residual)
  return(c(sm$aic, pv))
}
res2 <- as.data.frame(t(as.data.frame(apply(mutdatafeas, 2, classmodel, averages, rownames(mutdatafeas)))))
#head(res2[order(res2[,2]),])
colnames(res2) <- c("AIC","pval")
res2$adj_pval <- p.adjust(res2$pval, method="bonferroni")
write.table(res2, file=logisticOut, sep="\t", quote=FALSE)


if (debug == "yes") {
  save.image(file=paste0(linearOut,'.debug','.RData'))
}

