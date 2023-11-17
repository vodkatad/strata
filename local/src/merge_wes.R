#egrassi@godot:/mnt/trcanmed/snaketree/prj/strata/local/share/data/recuperi_wes$ cat altered_snp.bed altered_indel.bed  > all_hg18
#egrassi@godot:/mnt/trcanmed/snaketree/prj/strata/local/share/data/recuperi_wes$ cat snps_hg38.bed indels_hg38.bed  > all_hg38
setwd('/mnt/trcanmed/snaketree/prj/strata/local/share/data/recuperi_wes')
hg18 <- read.table('all_hg18', sep="\t", header=FALSE, stringsAsFactors = FALSE)
hg38 <- read.table('all_hg38', sep="\t", header=FALSE, stringsAsFactors = FALSE)

hg18$id <- paste0(hg18$V1,":", hg18$V3, ":", hg18$V4)
hg18$id <- gsub('-',":", hg18$id)
hg38$id <- paste0(hg38$V1,":", hg38$V3, ":", hg38$V4)
hg38$id <- gsub('-',":", hg38$id)
map <- data.frame(hg18=hg18$id, hg38=hg38$id)


muts <- read.table('altered_def.tsv', sep="\t", header=FALSE, stringsAsFactors = FALSE)
cases <- read.table('cases_def.tsv', sep="\t", header=FALSE, stringsAsFactors = FALSE)

data <- data.frame(smodel=muts$V2, mut=muts$V7, VAF=muts$V16, stringsAsFactors = FALSE)

s1 <- strsplit(data$mut, "_", fixed=TRUE)
chr <- sapply(s1, function(x) {x[[1]][1]})
coord <- sapply(s1, function(x) {x[[2]][1]})
ref <- sapply(s1, function(x) {x[[3]][1]})
alt <- sapply(s1, function(x) {if (length(x) == 4) { x[[4]][1] }else {''}} )

s2 <- strsplit(coord, "-", fixed=TRUE)
coorde <- sapply(s2, function(x) {x[[2]][1]})

data$id <- paste0(chr,":", coorde, ":", ref, ":", alt)
data$mut <- NULL

m <- merge(data, map, by.x="id", by.y="hg18")

stopifnot(nrow(m)==nrow(data))

library(reshape)
wide <- cast(data=m, formula = "hg38 ~ smodel", value='VAF')
wide[is.na(wide)] <- 0
colnames(wide)[1] <- 'ID'

toadd <- setdiff(cases$V1, unique(data$smodel))
nom <- data.frame(matrix(0, nrow=nrow(wide), ncol=length(toadd)))
colnames(nom) <- toadd
wes <- cbind(wide, nom)

hr <- read.table('../new_gatk/merged.table_nomultiallele_CRC00578', sep="\t", header=TRUE, stringsAsFactors = FALSE)

intersect(hr$id, wes$id)
intersect(colnames(hr), colnames(wes))

nom <- data.frame(matrix(0, nrow=nrow(hr), ncol=ncol(wes)-1))
nom$id <- hr$ID
nom <- nom[, c(ncol(nom), seq(1, ncol(nom)-1))]
colnames(nom) <- colnames(wes)
wes_added <- rbind(wes, nom)


nom <- data.frame(matrix(0, nrow=nrow(wes), ncol=ncol(hr)-1))
nom$id <- wes$ID
nom <- nom[, c(ncol(nom), seq(1, ncol(nom)-1))]
colnames(nom) <- colnames(hr)
hr_added <- rbind(hr, nom)

mtot <- merge(wes_added, hr_added, by="ID")

nrow(wes) + nrow(hr) == nrow(mtot)
ncol(wes) + ncol(hr) - 1== ncol(mtot)

write.table(mtot, file='/mnt/trcanmed/snaketree/prj/strata/local/share/data/merge_hr_wes/merged.tsv', sep="\t", quote=FALSE, row.names=FALSE)
