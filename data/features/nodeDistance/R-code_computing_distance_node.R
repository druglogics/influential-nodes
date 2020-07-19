#This script is used to calculate minimal distance of any node in the network to a drug target and output node

#import R-package for network analysis

if (!require("dodgr")) {
  install.packages("dodgr", dependencies = TRUE)
  library(dodgr)
}

data = read.csv("../models/insilicomodel_tab.sif", sep = "\t", header = FALSE)
colnames(data) = c("from", "d", "to")
data = data[,c("from","to", "d")]

#set the weight of all interactions to one
data$d = 1

#compute the distance of any node to all other nodes in network
distance = as.data.frame(dodgr_dists (data))

result = distance[,c("MAP3K7", "AKT_f", "MAPK14", "GSK3_f", "MEK_f", "PIK3CA", "CTNNB1", "JNK_f", "RSK_f",
                     "IKBKB", "TGFBR1", "ACVR1", "JAK_f", "CK1_f", "MYC", "PTEN", "STAT3", "PDPK1", "ROCK1", "SYK",
  "Prosurvival", "Antisurvival")]
result = result[c(1:144),]

#calculate shortest distance to a target or output
result$Shortest_distance_target = apply(result[, 1:20], 1, min)
result$Shortest_distance_output = apply(result[, 21:22], 1, min)

write.table(result, "distance_nodes_target_outputs.txt",sep = "\t", quote = FALSE)
