library(readr)
library(stringr)

load("~/Dropbox/meta_PKN_BIGG/results/meta_PKN.RData")

metab_ttop_tumour_vs_healthy <- as.data.frame(
  read_csv("data/cosmos_inputs_test/cosmos_paper/metab_ttop_tumour_vs_healthy.csv"))

metab_name_to_HMDB <- as.data.frame(
  read_csv("support/cosmos_test/metab_name_to_HMDB.csv"))

top_metabs <- metab_ttop_tumour_vs_healthy[metab_ttop_tumour_vs_healthy$P.Value <= 0.1,]

names(top_metabs)[1] <- "name"
top_metabs <- merge(top_metabs,metab_name_to_HMDB, by = "name")

top_metabs <- top_metabs[,c(8,2:7)]
top_metabs$HMDB <- paste("Metab__",top_metabs$HMDB, sep = "")

comps <- unique(c(str_extract(meta_PKN$source, "_[a-z]$"), str_extract(meta_PKN$target, "_[a-z]$")))
comps <- comps[!is.na(comps)]

top_df_list <- list()
for(comp in comps)
{
  df <- top_metabs
  df[,1] <- paste(df[,1],comp,sep = "")
  top_df_list[[comp]] <- df
}

top_df <- do.call(rbind, top_df_list)

write_csv(top_df, file = "data/cosmos_inputs_test/cosmos_paper/metab_ready.csv")
