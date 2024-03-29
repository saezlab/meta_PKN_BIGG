library(readr)
library(metaboliteIDmapping)

### Import the high confidence interactions of STITCH
STITCH <- as.data.frame(read_delim("data/9606.actions.v5.0.tsv", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE))

STITCH <- STITCH[STITCH$mode == "activation" | STITCH$mode == "inhibition",]

####THIS SECTION TO REMOVE TEXTMINING, download this file from STITCHdb website
X9606_protein_chemical_links_detailed_v5_0 <- as.data.frame(read_delim("data/9606.protein_chemical.links.detailed.v5.0.tsv", 
                                                                       "\t", escape_double = FALSE, trim_ws = TRUE))
threshold <- 700

not_text_mining <- X9606_protein_chemical_links_detailed_v5_0[X9606_protein_chemical_links_detailed_v5_0$combined_score >= threshold,]
rm(X9606_protein_chemical_links_detailed_v5_0)

not_text_mining <- not_text_mining[not_text_mining$experimental >= threshold | not_text_mining$database >= threshold,]

not_text_mining$ID <- paste(not_text_mining$chemical, not_text_mining$protein , sep = "_")
not_text_mining$ID_reverse <- paste(not_text_mining$protein, not_text_mining$chemical, sep = "_")

### We only care about allosteric interactions
STITCH <- STITCH[STITCH$a_is_acting,]
STITCH$ID <- paste(STITCH$item_id_a, STITCH$item_id_b, sep = "_")
STITCH <- STITCH[STITCH$ID %in% not_text_mining$ID | STITCH$ID %in% not_text_mining$ID_reverse,]
STITCH <- STITCH[,-7]
####END OF REMOVE TEXTMINING

### We need to map the Ensembl Ids to NCBI gene ids

## We make a vector of uniue ensembl IDs
prots <- unique(c(STITCH$item_id_a, STITCH$item_id_b))
prots <- prots[grepl("9606[.]ENSP", prots)]

## We remove the taxon ID so the IDs can be mapped
prots <- as.data.frame(cbind(prots, gsub("9606[.]","",prots)))

##We use biomart package to map the ensembl ids to NCBI ones
library(biomaRt)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

G_list <- getBM(filters = "ensembl_peptide_id", 
                attributes = c("ensembl_peptide_id",'hgnc_symbol','entrezgene_id', "description"),
                values = prots$V2, mart = ensembl)
names(G_list)[1] <- "V2"

## We put all the ID versions into one dataframe
prots <- merge(prots, G_list, by = "V2")
prots <- prots[prots$hgnc_symbol != "",]

## We create a named vector to make the id conversion more efficient
prots_vec <- prots$hgnc_symbol
names(prots_vec) <- prots$prots

## Ids are converted in the STITCH interaciton dataframe
for(i in 1:2)
{
  for(j in 1:length(STITCH[,1]))
  {
    if(STITCH[j,i] %in% names(prots_vec))
    {
      STITCH[j,i] <- prots_vec[STITCH[j,i]]
    }
    else
    {
      STITCH[j,i] <- gsub("CID[a-z]0*","Metab__",STITCH[j,i])
    }
  }
}

### We converted the interactions of STITCH into a SIF format
STITCH$sign <- ifelse(STITCH$action == "inhibition", -1, 1)
STITCH <- STITCH[grepl("Metab__",STITCH$item_id_a),]

STITCH <- STITCH[,c(1,7,2)]
names(STITCH) <- c("source","sign","target")

CIDs <- unique(as.character(unlist(STITCH[,c(1,3)])))
CIDs <- CIDs[grepl("Metab__",CIDs)]
CIDs <- gsub("Metab__","",CIDs)


##Convert CID to HMDB Id when available
metabolitesMapping <- metabolitesMapping[which(metabolitesMapping$CID %in% CIDs),]
metabolitesMapping <- metabolitesMapping[!is.na(metabolitesMapping$HMDB),]

metabolitesMapping_vec <- paste("Metab__",metabolitesMapping$HMDB, sep = "")
names(metabolitesMapping_vec) <- paste("Metab__",metabolitesMapping$CID, sep ="")

for(i in c(1,3))
{
  STITCH[,i] <- sapply(STITCH[,i], function(x, metabolitesMapping_vec)
    {
    if(x %in% names(metabolitesMapping_vec))
    {
      return(metabolitesMapping_vec[x])
    } else
    {
      return(x)
    }
  }, simplify = T, metabolitesMapping_vec = metabolitesMapping_vec)
}

clean_omnipath_PKN <- as.data.frame(read_csv("results/clean_omnipath_PKN.csv"))
omni_proteins <- unique(as.character(unlist(clean_omnipath_PKN[,c(1,2)])))

STITCH <- STITCH[which(STITCH$target %in% omni_proteins),]
STITCH <- unique(STITCH)

STITCH$source <- paste(STITCH$source, "_c", sep = "")
STITCH <- STITCH[,c(1,3,2)]

### Save the SIF network as csv
write_csv(STITCH,"results/STITCH_filtered.csv")
