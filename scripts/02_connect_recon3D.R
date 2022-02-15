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


###

enzyme_reacs <- unique(c(recon3D_bigg_no_cofactor$source, recon3D_bigg_no_cofactor$target))
enzyme_reacs <- enzyme_reacs[grepl("^Gene",enzyme_reacs)]

enzyme_reacs_reverse <- enzyme_reacs[grepl("_reverse",enzyme_reacs)]
enzyme_reacs <- enzyme_reacs[!grepl("_reverse",enzyme_reacs)]

new_df_list <- sapply(enzyme_reacs, function(enzyme_reac, recon3D_bigg_no_cofactor){
  df <- recon3D_bigg_no_cofactor[which(recon3D_bigg_no_cofactor$source == enzyme_reac | recon3D_bigg_no_cofactor$target == enzyme_reac),]
  if(dim(df)[1] < 2)
  {
    return(NA)
  } else
  {
    if(dim(df)[1] < 3)
    {
      return(df)
    } else
    {
      for(i in 1:dim(df)[1])
      {
        if(grepl("Metab__",df[i,1]))
        {
          counterpart <- which(gsub("_[a-z]$","",df[,2]) == gsub("_[a-z]$","",df[i,1]))
          if(length(counterpart) > 0)
          {

            df[i,2] <- paste(df[i,2],paste("_TRANSPORTER",i,sep = ""),sep = "")
            df[counterpart,1] <- paste(df[counterpart,1],paste("_TRANSPORTER",i,sep = ""),sep = "")
          }
        }
      }
      return(df)
    }
  }
}, recon3D_bigg_no_cofactor = recon3D_bigg_no_cofactor)
new_df <- as.data.frame(do.call(rbind,new_df_list))


new_df_list <- sapply(enzyme_reacs_reverse, function(enzyme_reac_reverse, recon3D_bigg_no_cofactor){
  df <- recon3D_bigg_no_cofactor[which(recon3D_bigg_no_cofactor$source == enzyme_reac_reverse | recon3D_bigg_no_cofactor$target == enzyme_reac_reverse),]
  if(dim(df)[1] < 2)
  {
    return(NA)
  } else
  {
    if(dim(df)[1] < 3)
    {
      return(df)
    } else
    {
      for(i in 1:dim(df)[1])
      {
        if(grepl("Metab__",df[i,1]))
        {
          counterpart <- which(gsub("_[a-z]$","",df[,2]) == gsub("_[a-z]$","",df[i,1]))
          if(length(counterpart) > 0)
          {
            transporter <- gsub("_reverse","",df[i,2])
            transporter <- paste(transporter,paste(paste("_TRANSPORTER",i,sep = ""), "_reverse",sep = ""),sep = "")
            df[i,2] <- transporter
            df[counterpart,1] <- transporter
          }
        }
      }
      return(df)
    }
  }
}, recon3D_bigg_no_cofactor = recon3D_bigg_no_cofactor)
new_df_reverse <- as.data.frame(do.call(rbind,new_df_list))

recon3D_bigg_no_cofactor <- as.data.frame(rbind(new_df, new_df_reverse))

recon3D_bigg_no_cofactor <- recon3D_bigg_no_cofactor[complete.cases(recon3D_bigg_no_cofactor),]

###Build the connectors to omnipath
elements <- unique(as.character(unlist(recon3D_bigg_no_cofactor)))
elements <- elements[!grepl("Metab__",elements)]
elements <- elements[!grepl("__[0-9]",elements)]

connectors_list <- sapply(elements, function(x){
  # prefixe <- gsub("__.*","",x)
  # suffixe <- ifelse(grepl("_reverse$",x),"_reverse","")
  genes <- gsub(".*__","",x)
  genes <- gsub("_TRANSPORTER[0-9]+","",genes)
  genes <- gsub("_reverse$","",genes)
  print(genes)
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
