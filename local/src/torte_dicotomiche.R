data_markers <- read.table('/mnt/trcanmed/snaketree/prj/strata/dataset/figures/POLD1_RAD51_quartili_PFS_fwup.tsv', sep="\t", header=T, stringsAsFactors = F)
resp <- read.table('/mnt/trcanmed/snaketree/prj/strata/local/share/data/IRIS_responses.tsv', sep="\t", header=T, stringsAsFactors = F)

m <- merge(data_markers, resp, by="ID")

stopifnot(nrow(m)==82)

m$recist <- m$ORR
m$ORR <- ifelse(m$recist %in% c('PR', 'CR'), 'yes', 'no')

m$RAD51_dico <- ifelse(m$RAD51_class %in% c('0', '1+'), 'low', 'high')
m$POLD1_dico <- ifelse(m$POLD1_class %in% c('0', '1+'), 'low', 'high')

rad51 <- table(m$ORR, m$RAD51_dico)

pold1 <- table(m$ORR, m$POLD1_dico)

fisher.test(rad51)
chisq.test(rad51, simulate.p.value = T)

fisher.test(pold1)
chisq.test(pold1, simulate.p.value = T)

library(reshape)
library(RColorBrewer)
#mp <- melt(as.data.frame(rad51))
mp <- melt(as.data.frame(pold1))
title <- 'ORR'

df <- mp[mp$Var2 == "high",]
ggplot(df, aes(x = "", y = value, fill = Var1)) +
  geom_col(color = "black") +
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c('darkred','darkblue'))+theme_bw(base_size=15)+ggtitle('High')+
  labs(fill=title)

df <- mp[mp$Var2 == "low",]
ggplot(df, aes(x = "", y = value, fill = Var1)) +
  geom_col(color = "black") +
  geom_text(aes(label = value),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c('darkred','darkblue'))+theme_bw(base_size=15)+ggtitle('Low')+
  labs(fill=title)
