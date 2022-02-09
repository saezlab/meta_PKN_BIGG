library(readr)
library(org.Hs.eg.db)
library(stringr)
library(cosmosR)

load("~/Dropbox/meta_PKN_BIGG/results/meta_PKN.RData")

metab_input_COSMOS <- as.data.frame(read_csv("data/cosmos_inputs_test/cosmos_paper/metab_ready.csv"))
signaling_input_COSMOS <- as.data.frame(read_csv("data/cosmos_inputs_test/cosmos_paper/signaling_ready.csv"))

ttop_rna <- as.data.frame(read_csv("data/cosmos_inputs_test/cosmos_paper/RNA_ttop_ready.csv"))
ttop_rna <- ttop_rna[complete.cases(ttop_rna),]
RNA_input <- ttop_rna[,"t"]
names(RNA_input) <- ttop_rna$ID
#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options()

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable
my_options$solver <- "cplex" #or cbc
my_options$timelimit <- 1800
my_options$mipGAP <- 0.05
my_options$threads <- 6

metab_input_COSMOS_vec <- metab_input_COSMOS$t
names(metab_input_COSMOS_vec) <- metab_input_COSMOS$HMDB

signaling_input_COSMOS_vec <- signaling_input_COSMOS$NES
names(signaling_input_COSMOS_vec) <- signaling_input_COSMOS$ID

metab_input_COSMOS_vec <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input_COSMOS_vec, meta_PKN)
signaling_input_COSMOS_vec <- cosmosR:::filter_input_nodes_not_in_pkn(signaling_input_COSMOS_vec, meta_PKN)

meta_PKN <- meta_PKN[,c(1,3,2)]
names(meta_PKN) <- c("source", "interaction", "target")

comp <- "_c"

metab_input_COSMOS_vec <- metab_input_COSMOS_vec[grepl(comp,names(metab_input_COSMOS_vec))]

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_PKN,
                                                      signaling_data = signaling_input_COSMOS_vec,
                                                      metabolic_data = metab_input_COSMOS_vec,
                                                      diff_expression_data = RNA_input,
                                                      maximum_network_depth = 8,
                                                      remove_unexpressed_nodes = T,
                                                      CARNIVAL_options = my_options
                                                      
)

my_options$timelimit <- 7200

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)
