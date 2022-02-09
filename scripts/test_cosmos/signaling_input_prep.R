library(readr)

load("~/Dropbox/meta_PKN_BIGG/results/meta_PKN.RData")

kinase_activities <- as.data.frame(
  read_csv("data/cosmos_inputs_test/cosmos_paper/phospho_kinase_activities.csv"))
names(kinase_activities) <- c("ID","NES")

TF_activities <- as.data.frame(
  read_csv("data/cosmos_inputs_test/cosmos_paper/RNA_TF_scores.csv"))
names(TF_activities) <- c("ID","NES")

kinase_activities_top <- kinase_activities[order(abs(kinase_activities$NES), decreasing = T),]
kinase_activities_top <- kinase_activities_top[1:20,]

TF_activities_top <- TF_activities[order(abs(TF_activities$NES), decreasing = T),]
TF_activities_top <- TF_activities_top[1:80,]

signaling_ready <- as.data.frame(rbind(kinase_activities_top,TF_activities_top))

write_csv(signaling_ready, file = "data/cosmos_inputs_test/cosmos_paper/signaling_ready.csv")
