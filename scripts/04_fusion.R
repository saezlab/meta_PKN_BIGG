library(readr)

clean_omnipath_PKN <- as.data.frame(read_csv("results/clean_omnipath_PKN.csv"))

recon3D_connected <- as.data.frame(read_csv("results/recon3D_connected.csv"))

STITCH_filtered <- as.data.frame(read_csv("results/STITCH_filtered.csv"))

meta_PKN <- as.data.frame(rbind(clean_omnipath_PKN,STITCH_filtered))
recon3D_connected$sign <- 1

meta_PKN <- as.data.frame(rbind(meta_PKN, recon3D_connected))

meta_PKN[is.na(meta_PKN$source) | is.na(meta_PKN$target),]

save(meta_PKN, file = "results/meta_PKN.RData")

