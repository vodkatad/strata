library(ggplot2, quietly=TRUE)
# irinotecan.txt derives from 'WATERFALL_PLOT_def.xlsx' raw_data sheet, with desired colums. need at least two (baseline-ref and what
# we study, eg. w0 - w5).
# data@clo:~/Dropbox/work/strata/data$ awk -F"\t" '$1!="" && $1!="Genealogy ID"{print}' irinotecan.txt  > filter1

pdx_id_col <- snakemake@params[["pdx"]]
ref_week <- snakemake@params[["ref"]]
target_week <- snakemake@params[["target"]]
inputData <- snakemake@input[["volumes"]]
outputPlot <- snakemake@output[["waterfall"]]
outputMeans <- snakemake@output[["means"]]
debug <- snakemake@params[["debug"]]

#TODO makes sense to have a "possible" filter on n of mice? or a priori?
#n_5week_data <- sapply(levels(data$case), function(x) { d <- data[data$case==x & !is.na(data$vol_5w),]; nrow(d) })
#d <- as.data.frame(n_5week_data)
#ggplot(d, aes(x=-n_5week_data))+geom_histogram(aes(y = cumsum(..count..)), binwidth = 1, boundary = 0,
#                                               color = "black", fill = "white")
#wanted <- names(n_5week_data[n_5week_data >=4])


data <- read.table(inputData, sep="\t", header=TRUE)
#colnames(data) <- c("pdx_long", "vol0", "vol_5w", "vol_6w")
data$case <- as.factor(substr(data[,pdx_id_col], 1,7))


avg <- sapply(levels(data$case), function(x) { 
  d <- data[data$case==x,]; 
  c(mean(d[,ref_week], na.rm=TRUE), mean(d[,target_week], na.rm=TRUE))
  })
averages <- as.data.frame(t(avg))
colnames(averages) <- c(ref_week,target_week)
averages <- averages[!is.na(averages[,target_week]),]
#=(F7/C7-1)*100
averages$perc <- ((averages[,target_week]/averages[,ref_week])-1)*100
averages$case <- rownames(averages)
ggplot(averages, aes(y=perc,x=reorder(case, -perc)))+geom_col()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(outputPlot)
write.table(averages, file=outputMeans, sep="\t", row.names =  FALSE, col.names=TRUE, quote=FALSE)

if (debug == "yes") {
  save.image(file=paste0(outputMeans,'.debug','.RData'))
}
