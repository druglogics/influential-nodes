#calculate nr of negative/postive in/out-going edges

data = read.csv("../../models/insilicomodel_tab.sif", sep = "\t", header = F)
data$V1 = as.character(data$V1)
data$V2 = as.character(data$V2)

node = as.data.frame(unique(data[,1]))

node = unique(rbind(as.data.frame(unique(data[,3]))))

colnames(node) = c("Node")

#Otdegree inhibit/activate
for(i in 1:nrow(node)){
  #identify node
  x = node[i,1]
  #calculate number of occurences
  node[i,2] = nrow(subset(data, data$V1 == x & data$V2 == "inhibit"))
}

for(i in 1:nrow(node)){
  #identify node
  x = node[i,1]
  #calculate number of occurences
  node[i,3] = nrow(subset(data, data$V1 == x & data$V2 == "activate"))
}

#Indegree inhibit/activate
for(i in 1:nrow(node)){
  #identify node
  x = node[i,1]
  #calculate number of occurences
  node[i,4] = nrow(subset(data, data$V3 == x & data$V2 == "inhibit"))
}

for(i in 1:nrow(node)){
  #identify node
  x = node[i,1]
  #calculate number of occurences
  node[i,5] = nrow(subset(data, data$V3 == x & data$V2 == "activate"))
}

colnames(node) = c("Node", "Outdegree.inhibit", "Outdegree.activate", "Indegree.inhibit", "Indegree.activate")


#read in other network feature files and combined them
#Cytoscape
cytoscape = read.csv("../Cytoscape/20190201_Cytoscape_features.csv", sep = ",", header = T)

##keep columns of interest
cytoscape = cytoscape[,c("canonicalName", "ClosenessCentrality", "BetweennessCentrality", 
                         "Drug.target",	"Indegree",	"Outdegree",	"NeighborhoodConnectivity",
                         "ClusteringCoefficient",	"AverageShortestPathLength")]
colnames(cytoscape)[1] = "Node"
cytoscape$Drug.target = as.logical(cytoscape$Drug.target)
cytoscape$Drug.target[is.na(cytoscape$Drug.target)] = FALSE

#add HGNC symbols
HGNC = read.csv("../../features/20190130_nodes_to_HGNC.txt", sep = "\t")
cytoscape = merge(cytoscape, HGNC, by = "Node", all = T)

#Shortest_distance_target and Shortest_distance_output
distance = read.csv("../nodeDistance/distance_nodes_target_outputs.txt", sep = "\t")
distance$Node = rownames(distance)

cytoscape = merge(cytoscape, distance[,c("Node", "Shortest_distance_target", "Shortest_distance_output")], by = "Node", all = T)

#PCI
PCI = read.csv("../nodeDistance/PCI.txt", sep = "\t", skip = 1)
colnames(PCI)[1] = "Node"

cytoscape = merge(cytoscape, PCI, by = "Node", all = T)

#add outdegree/indegree activate/inhibit
cytoscape = merge(cytoscape, node, by = "Node", all = T)


#output node feature file
write.table(cytoscape,"../node_features.txt", sep = "\t", quote = F, row.names = F)
