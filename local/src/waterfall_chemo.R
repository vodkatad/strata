# what if we add HR from wes to ours?
  
#d <- read.table('www', sep="\t", stringsAsFactors=F)
#colnames(d) <- c('omodel', 'class')
#d$model <- gsub('CRC', '', d$omodel)
#pad <- 4 - nchar(d$model)
#for (i in seq(1, length(pad))) {
#  d$model[i] <- paste0('CRC',  paste0(rep('0', pad[i]), collapse=''), d$model[i])
#}
#head(d)

# not in suppl of Velcu but in DTB, do not know why.
toremove <- c( 'CRC0257',
               'CRC0307',
               'CRC0316',
               'CRC0539',
               'CRC0598' )

# shallow biobanca
biobanca <- read.table('/scratch/trcanmed/biobanca/local/share/data/shallowseq/map_lmx', sep="\t", header=F, stringsAsFactors = F)
biobanca$model <- substr(biobanca$V2, 0, 7)
# already requested shallow
already <- read.table('/mnt/trcanmed/snaketree/prj/strata/local/share/data/xenos_withoutshallow_withHRtargetedAndResponse.tsv', sep="\t", header=T, stringsAsFactors = F)

length(intersect(already$short_genealogy, biobanca$model))

w3 <- read.table('/mnt/trcanmed/snaketree/prj/strata/dataset/figures/w3_waterfall.tsv', sep="\t", header=F, stringsAsFactors = F)
colnames(w3) <- c('model', 'perc')
hr <- read.table('/mnt/trcanmed/snaketree/prj/strata/dataset/figures/all_cases_hr.tsv', sep="\t", header=F, stringsAsFactors = F)
colnames(hr) <- c('model')
hr$seq <- 'hr'

m <- merge(w3, hr, by="model", all.x=T)
m[is.na(m$seq), 'seq'] <-'miss'

wes <- read.table('/scratch/trcanmed/pdxopedia/local/share/data/old_seq/wes_cases', sep="\t", header=F, stringsAsFactors = F)
colnames(wes) <- c('model')
wes$seq2 <- 'wes'

m2 <- merge(m, wes, by="model", all.x=T)
m2[is.na(m2$seq2), 'seq2'] <-'miss'

m2$seqq <- ifelse(m2$seq == "hr", 'hr', ifelse(m2$seq2=="wes", 'wes', 'miss'))


#### list for suppl 4
total_cases_with_mut <- m2[m2$seqq != "miss",]
total_cases_with_mut2 <- total_cases_with_mut[!total_cases_with_mut$model %in% toremove,]
write.table(total_cases_with_mut2, file="/mnt/trcanmed/snaketree/prj/strata/local/share/data/appello_mut.tsv", sep="\t", quote=FALSE, row.names=FALSE)
##########
library(ggplot2)
ggplot(data=m2, aes(x=reorder(model, -perc), fill=seqq, y=perc))+geom_col()+theme_bw()+scale_fill_manual(values=c('forestgreen','grey','darkgoldenrod'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# add shallow to MR
biobanca$biob <-"biob"
biobanca$V1 <- NULL
biobanca$V2 <- NULL
m3 <- merge(m2, biobanca, by="model", all.x=TRUE)
m3[is.na(m3$biob), 'biob'] <-'miss'

# already
colnames(already) <- 'model'
already$already <- 'already'
m4 <- merge(m3, already, by="model", all.x=TRUE)
m4[is.na(m4$already), 'already'] <-'miss'
# removetoremove

toadd <- m4[m4$seqq=="wes"& m4$biob== "miss",]

toadd2 <- toadd[!toadd$model %in% toremove,]
setdiff(toadd$model, toadd2$model)


m2$addvWES <- ifelse(m2$model %in% toadd2$model, 'yes', 'no')
ggplot(data=m2, aes(x=reorder(model, -perc), fill=addvWES, y=perc))+geom_col()+theme_bw()+scale_fill_manual(values=c('grey','darkgoldenrod'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


m4$addvWES <- ifelse(m4$model %in% toadd2$model, 'yes', 'no')

mutdata <- m4[m4$seqq!="miss",'model']

new <- m4[m4$already!="miss" | m4$addvWES=="yes", 'model']
biob <- m4[m4$biob=="biob", 'model']
length(intersect(new, biob))
length(intersect(mutdata, c(new, biob)))

length(setdiff(mutdata, c(new, biob)))

mm <- setdiff(c(new, biob), mutdata)
m4[m4$model %in% mm,]
