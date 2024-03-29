library(OmnipathR)
library(org.Hs.eg.db)
library(readr)

full_pkn <- as.data.frame(import_omnipath_interactions())
full_pkn <- full_pkn[!is.na(full_pkn$references),]

clean_PKN <- full_pkn[full_pkn$consensus_stimulation == 1 | full_pkn$consensus_inhibition == 1,]

clean_PKN$sign <- clean_PKN$consensus_stimulation - clean_PKN$consensus_inhibition

clean_PKN <- clean_PKN[,c(3,4,16)]

clean_PKN_supp <- clean_PKN[clean_PKN$sign == 0,]
clean_PKN_supp$sign <- -1
clean_PKN[clean_PKN$sign == 0,"sign"] <- 1

clean_PKN <- as.data.frame(rbind(clean_PKN, clean_PKN_supp))

names(clean_PKN) <- c("source","target","sign")
write_csv(x = clean_PKN, file = "results/clean_omnipath_PKN.csv")

      