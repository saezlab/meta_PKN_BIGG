library(readr)
library(org.Hs.eg.db)
library(dplyr)

RNA_ttop_tumorvshealthy <- as.data.frame(
  read_csv("data/cosmos_inputs_test/cosmos_paper/RNA_ttop_tumorvshealthy.csv"))

mapping_vec <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, as.character(RNA_ttop_tumorvshealthy$ID), 'SYMBOL', 'ENTREZID')
RNA_ttop_tumorvshealthy$ID <- mapping_vec


RNA_ttop_tumorvshealthy <- RNA_ttop_tumorvshealthy %>% group_by(ID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
RNA_ttop_tumorvshealthy <- as.data.frame(RNA_ttop_tumorvshealthy)

write_csv(RNA_ttop_tumorvshealthy, file = "data/cosmos_inputs_test/cosmos_paper/RNA_ttop_ready.csv")
