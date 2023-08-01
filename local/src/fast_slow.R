library(ggplot2)
data <- read.table('/mnt/trcanmed/snaketree/prj/strata/local/share/data/shrink_data_chemio.tsv', sep="\t", header=TRUE)
data$sickness <- ifelse(data$Sickness != "", 'Y', 'N')
data$size <- ifelse(data$Size != "", 'Y', 'N')

data$Shrink <- factor(data$Shrink, levels=c('non-shrink', 'slow', 'fast'))
ggplot(data=data, aes(x=Shrink, fill=Response.Class.3wks))+geom_bar()+theme_bw()
ggplot(data=data, aes(x=Shrink, fill=Response.Class.6wks))+geom_bar()+theme_bw()

ggplot(data=data, aes(x=Shrink, y=n.animal.3WKS))+geom_boxplot()+geom_jitter(aes(color=sickness))+theme_bw()

ggplot(data=data, aes(x=Shrink, y=n.animal.6WKS))+geom_boxplot()+geom_jitter(aes(color=sickness))+theme_bw()

ggplot(data=data, aes(x=Shrink, y=n.animal.6WKS))+geom_boxplot()+geom_jitter(aes(color=size))+theme_bw()


ggplot(data=data, aes(x=Shrink, y=n.animal.3WKS))+geom_boxplot()+geom_jitter(aes(color=sickness))+theme_bw()


f30 <- read.table('/mnt/trcanmed/snaketree/prj/strata/local/share/data/ff30.txt', sep="\t", header=TRUE, stringsAsFactors = F)

m <- merge(data, f30, by.x="Model", by.y="CASE")
nrow(m)
nrow(f30)
nrow(data)

data <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/treats/february2023/irino_feb23.txt',  sep="\t", header=TRUE, stringsAsFactors = F)
m <- merge(data, f30, by="CASE")
nrow(m)
nrow(f30)
nrow(data)

#> data[data$CASE=="",]
#CASE MUT X3WKS X4WKS X5WKS X6WKS Response.Class.3wks KRAS NRAS PI3KCA X4ple.WT BRAF KRAS.1 NRAS.1 PI3KCA.1 X4ple.WT.1 BRAF.1 KRAS.2 NRAS.2 PI3KCA.2 X4ple.WT.2 BRAF.2 KRAS.3 NRAS.3 PI3KCA.3
#126             NA    NA    NA    NA                       NA   NA     NA       NA   NA     NA     NA       NA         NA     NA     NA     NA       NA         NA     NA     NA     NA       NA
#127             NA    NA    NA    NA                       NA   NA     NA       NA   NA     NA     NA       NA         NA     NA     NA     NA       NA         NA     NA     NA     NA       NA
#X4ple.WT.3 BRAF.3  X X.1 X.2 X.3 X.4 X.5 X.6 X.7 X.8 X.9 X.10 X.11 X.12 X.13 X.14 X.15 X.16
#126         NA     NA NA  NA  NA  NA  NA  NA  NA  NA  NA  NA   NA   NA   NA   NA   NA   NA   NA
#127         NA     NA NA  NA  NA  NA  NA  NA  NA  NA  NA  NA   NA   NA   NA   NA   NA   NA   NA
# Unknown garbage?

m$recist3w <-  as.character(ifelse(m$X3WKS < -50, 'PR', ifelse(m$X3WKS > 35, 'PD','SD')))
m$recist6w <-  as.character(ifelse(m$X6WKS < -50, 'PR', ifelse(m$X6WKS > 35, 'PD','SD')))
m[is.na(m$recist6w),'recist6w'] <- 'PD'
table(m$Response.Class.3wks==m$recist3w)

m[m$recist3w == 'PD' & m$recist6w == 'PR',]

mm <- m[m$recist3w != m$recist6w,]

table(mm$recist3w, mm$recist6w)
mmf30 <- mm[mm$ff30=="ff30",c('CASE', 'MUT','X3WKS','X4WKS','X5WKS','X6WKS', 'recist3w', 'recist6w')]
