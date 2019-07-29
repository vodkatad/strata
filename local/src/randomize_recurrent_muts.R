library(scales)
library(ggplot2)
library(RColorBrewer)
library(reshape)

datafile <- snakemake@input[["data"]]
possiblemutfile <- snakemake@input[["possiblemut"]]
plotrand <- snakemake@output[["plotrand"]]
plotaf <- snakemake@output[["plotaf"]]
table <- snakemake@output[["table"]]
nrand <- snakemake@params[["nrand"]]
debug <- snakemake@params[["debug"]]

if (debug == "yes") {
  save.image(file=paste0(table,'.debug','.RData'))
}

d <- read.table(datafile, header=TRUE, row.names=1)
countdf <- as.data.frame(table(apply(d,1, function(x) {sum(x>0)})))
write.table(countdf, file=table, sep="\t", quote=FALSE)
countdf$count <- as.numeric(countdf$Var1)
md <- melt(d)
mdnz <- md[md$value != 0,]
ggplot(mdnz, aes(value, fill = variable)) + geom_histogram() + facet_wrap(~variable)
ggsave(plotaf)

nave <- apply(d,2, function(x) {sum(x>0)})
#possiblemut <- 47805824
pm <- read.table(possiblemutfile, header=FALSE)
possiblemut <- pm[1,1]
#47805824
nsamples <- ncol(d)
# 
randAndrea <- function(nave, nsamples) {
  onesample <- function(ave_mut, possiblemut) {
    size <- rnorm(1, mean(ave_mut), sd=(ave_mut))
    while (size < 1) {
      size <- rnorm(1, mean(ave_mut), sd=(ave_mut))
    }
    #print(size)
    muts <- sample.int(possiblemut, size=size)
    res <- data.frame(row.names=muts, value=rep(1, size))
    names(res)[names(res) == "value"] <- replic
    replic <<- replic+1
    res
  }
  replic <- 0
  samples <- replicate(nsamples, onesample(nave, possiblemut), simplify=F)

  mymrcs <- function (dfs, ...) 
  {
      if (length(dfs) == 2) {
        m <- merge(dfs[[1]], dfs[[2]], all = TRUE, sort = FALSE, ...)
        rownames(m) <- m$Row.names
        m$Row.names <- NULL
      }
      else {
          m  <- merge(dfs[[1]], Recall(dfs[-1], ...), all = TRUE, sort = FALSE, 
              ...)
          rownames(m) <- m$Row.names
          m$Row.names <- NULL
      }
      m
  }
  
  
  mrec <- mymrcs(samples, by="row.names")
  mrec[is.na(mrec)] <- 0
  counts_randA <- apply(mrec,1, function(x) {sum(x>0)})
  res <- as.data.frame(table(counts_randA))
  #https://stackoverflow.com/questions/23121654/r-why-is-merge-recurse-failing
  #Reduce(function(x, y) x[y], list(dt1, dt2, dt3)) ?
  n <- nrow(res)
  res[,1] <- as.numeric(res[,1])
  while (nrow(res) < nsamples) {
    #print(n)
    n <- n + 1
    res <- rbind(res, c(n, 0))
  }
  res
}

manyrandsA <- replicate(nrand, randAndrea(nave, nsamples), simplify=FALSE)
togpA <- do.call(rbind, manyrandsA)
colnames(togpA)[1] <- "count"                                                        

onerandpupido <- function(alldata, samples, muts) {
  samplenoz <- function(alld, samples) {
    pick <- sample(alld, samples)
    while (all(pick == 0)) {
      pick <- sample(alld, samples)
    }
    pick
  }
  r <- t(replicate(muts, samplenoz(alldata, samples)))
  count <- apply(r,1, function(x) {sum(x>0)})
  res <- as.data.frame(table(count))
  res$count <- as.numeric(res$count)
  if (nrow(res) != samples) {
    res <- rbind(res, c(samples, 0))
  }
  res
}

alld <- unlist(d)
manyrands <- replicate(nrand, onerandpupido(alld, ncol(d), nrow(d)), simplify=FALSE)
togp <- do.call(rbind, manyrands)
togp$rand <- "real"
togpA$rand <- "rand"
wholedata <- rbind(togp, togpA)
cbPalette <- c("royalblue1", "seagreen4")
ggplot(wholedata, aes(as.factor(count), Freq))  + geom_boxplot(aes(colour=rand)) + geom_point(data=countdf, aes(x=as.factor(count), y=Freq),colour = "darkred", size=3)+theme_bw()+xlab("N. of mutated samples")+ylab("N. of muts")+ theme(legend.position = "bottom", text = element_text(size=15), axis.text=element_text(size=15))+scale_colour_manual(values=cbPalette, labels = c("Completely random mutations", "Shuffle real data") )+labs(colour="")
ggsave(plotrand)

if (debug == "yes") {
  save.image(file=paste0(table,'.debug','.RData'))
}
