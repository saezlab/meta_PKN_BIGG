library(readr)
library(stringr)

recon3D_bigg_no_cofactor <- as.data.frame(read_csv("data/recon3D_bigg_no_cofactor.csv"))

HMDB_full_mapping_long <- as.data.frame(
  read_csv("support/HMDB_full_mapping_long.csv"))

HMDB_full_mapping_long_vec <- paste("Metab__",HMDB_full_mapping_long$good_HMDB_id,sep = "")
names(HMDB_full_mapping_long_vec) <- paste("Metab__",HMDB_full_mapping_long$legacy_HMDB_id,sep = "")

HMDB_full_mapping_long_vec <- c(HMDB_full_mapping_long_vec,"Metab__HC00342" = "Metab__HMDB0000072")
for(i in 1:2)
{
  recon3D_bigg_no_cofactor[,i] <- sapply(recon3D_bigg_no_cofactor[,i], function(x,HMDB_full_mapping_long_vec){
    suffixe <- str_extract(x, "_[a-z]$")
    metab <- gsub("_[a-z]$","",x)
    
    if(metab %in% names(HMDB_full_mapping_long_vec))
    {
      metab <- paste(HMDB_full_mapping_long_vec[metab],suffixe,sep = "")
      return(metab)
    }
    return(x)
}, HMDB_full_mapping_long_vec = HMDB_full_mapping_long_vec)
}

###Build the connectors to omnipath
elements <- unique(as.character(unlist(recon3D_bigg_no_cofactor)))
elements <- elements[!grepl("Metab__",elements)]
elements <- elements[!grepl("__[0-9]",elements)]

connectors_list <- sapply(elements, function(x){
  # prefixe <- gsub("__.*","",x)
  # suffixe <- ifelse(grepl("_reverse$",x),"_reverse","")
  genes <- gsub(".*__","",x)
  genes <- gsub("_reverse$","",genes)
  if(grepl("_",genes))
  {
    elements <- str_split(string = genes, pattern = "_")[[1]]
    if(length(elements) < 10)
    {
      genes_connector_list <- sapply(elements,function(gene)
        {
        return(c(gene,x))
      })
      return(t(genes_connector_list))
    } 
  } else
  {
    return(c(genes,x))
  }
})

connectors_df <- as.data.frame(do.call(rbind,connectors_list))
names(connectors_df) <- c("source","target")

clean_omnipath_PKN <- as.data.frame(read_csv("results/clean_omnipath_PKN.csv"))

connectors_df <- connectors_df[which(connectors_df$source %in% clean_omnipath_PKN$source | connectors_df$source %in% clean_omnipath_PKN$target),]

recon3D_connected <- as.data.frame(rbind(recon3D_bigg_no_cofactor,connectors_df))

write_csv(recon3D_connected, file = "results/recon3D_connected.csv")
