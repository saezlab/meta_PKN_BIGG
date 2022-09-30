library(readr)

clean_omnipath_PKN <- as.data.frame(read_csv("results/clean_omnipath_PKN.csv"))

recon3D_connected <- as.data.frame(read_csv("results/recon3D_connected.csv"))

STITCH_filtered <- as.data.frame(read_csv("results/STITCH_filtered.csv"))

meta_PKN <- as.data.frame(rbind(clean_omnipath_PKN,STITCH_filtered))
recon3D_connected$sign <- 1

meta_PKN <- as.data.frame(rbind(meta_PKN, recon3D_connected))

meta_PKN[is.na(meta_PKN$source) | is.na(meta_PKN$target),]

meta_PKN <- unique(meta_PKN)

save(meta_PKN, file = "results/meta_PKN.RData")

meta_network <- meta_PKN
meta_network <- meta_network[,c(1,3,2)]
names(meta_network) <- c("source", "interaction", "target")

#probably erroneous interaction
meta_network <- meta_network[-which(meta_network$source == "PRKCA" & meta_network$target == "SRC"),]

#probably erroneous interaction
meta_network <- meta_network[-which(meta_network$source == "LTC4S"),] #I don't know where this interaction comes from, the sources are wrong (https://www.nature.com/articles/onc2008228)

meta_network <- meta_network[!(grepl("CAD_reverse",meta_network$source) | grepl("CAD_reverse",meta_network$target)) ,] #redHuman confirms that the reaction is actually not reversible

# save(meta_network, file = "../cosmosR/data/meta_network.RData")
