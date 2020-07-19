#Rscript to generate subgraphs that can be used to explore putative synergy mechanisms
library("tidyr")
library("dplyr")
library("igraph")

#read in AGS data for influential nodes per inversion and fixation
AGS_fix = read.csv("taskNodeAssessment/results/predictions/AGS_fixed.txt", sep = "\t")
AGS_inverted = read.csv("taskNodeAssessment/results/predictions/AGS_inverted.txt", sep = "\t")
AGS_fix$class = "fixed"
AGS_inverted$class = "inverted"

#WT model predictions
AGS_wt = read.csv("taskNodeAssessment/results/predictions/AGS_WT.txt", sep = "\t")

#combine fixed and inverted models and compare to WT
data = rbind(AGS_inverted, AGS_fix)
data = merge(data, AGS_wt, by = "drugs")
data = data[data$classification.x != data$classification.y,]

#read in network
network = read.csv("data/models/insilicomodel_tab.sif", sep = "\t", header = FALSE)
colnames(network) = c("from", "d", "to")
network = network[,c("from","to", "d")]

# Generate colors based on activity type:
g1 = graph.data.frame(network,directed=T)

plot(g1, layout=layout_with_kk(g1), edge.arrow.size=.2, vertex.size=15, 
     vertex.frame.color="black", vertex.label.color="black", 
     vertex.label.cex=0.8, edge.curved=0.2)

# Get the list of influential nodes and induce subgraphs
##define combination
comb = "PI-D1"

subv = subset(data[data$drugs == comb & data$classification.x != "n/a",], select = c("mutated_node"))

#combined with influential nodes
subv = separate(subv, "mutated_node", c("Nodes", "Activity"), sep = "%")
subv$Activity = as.numeric(subv$Activity)

subv = subv %>%
  group_by(Nodes) %>%
  summarize(Activity = mean(Activity, na.rm=FALSE))

#generate subgraph
g2 =  induced.subgraph(graph=g1,vids=subv$Nodes)

subv = subv[match(vertex_attr(g2)$name, subv$Nodes),]

vertex_attr(g2)$type = subv$Activity

vertex_attr(g2)$color[vertex_attr(g2)$type == "1"] ="green"
vertex_attr(g2)$color[vertex_attr(g2)$type == "0"] ="red"
vertex_attr(g2)$color[vertex_attr(g2)$type == "0.5"] ="orange"
vertex_attr(g2)$color[is.na(vertex_attr(g2)$type)] ="grey"
vertex_attr(g2)$color[vertex_attr(g2)$type == "changed"] = "lightblue"

plot(g2, edge.arrow.size=.2, 
     vertex.size=15, 
     vertex.frame.color="black", vertex.label.color="black", 
     vertex.label.cex=0.8, edge.curved=0.2,
     edge.color=ifelse(edge_attr(g2)$d =="activate", "green", "red"),
     layout=layout.fruchterman.reingold) #layout_as_tree, layout_nicely, layout.fruchterman.reingold

#select nodes with no outgoing edges
subv_2 = as.data.frame(unlist(V(g2)[ degree(g2, mode = "out") == 0]$name))

#specifcy whether they are influential upon KO (%0) or KI (%1)
##OBS needs to be manually specified
subv_2 = as.data.frame(c("AKT_f%1", "AKT_f%0", "RSK_f%0"))

#check which nodes are different in their activity when these nodes are knocked in or out
#compared to WT
WT_stable_state = read.csv("taskNodeAssessment/results/simulationsStates/AGS_WT_growth.txt", sep = "\t")
WT_stable_state = as.data.frame(strsplit(as.character(WT_stable_state[WT_stable_state$drugs == comb,5]), split = ""))
colnames(WT_stable_state) = c("WT")

WT_stable_state[,2] = read.csv("data/full_node_list.txt")

fixed_stable_state = read.csv("taskNodeAssessment/results/simulationsStates/AGS_fixed_growth.txt", sep = "\t", stringsAsFactors = F)
inverted_stable_state = read.csv("taskNodeAssessment/results/simulationsStates/AGS_inverted_growth.txt", sep = "\t", stringsAsFactors = F)
mutated_stable_state = rbind(fixed_stable_state, inverted_stable_state)

extract_nodes = function(mutated_stable_state, WT_stable_state, combination = comb){
  stable_state = data.frame(ncol = 1)
  stable_state_select = data.frame(matrix(ncol = 2))
  
  for(nodes in as.character(subv_2[,1])){
    stable_state = as.data.frame(strsplit(as.character(
      mutated_stable_state[mutated_stable_state$drugs == combination & mutated_stable_state$mutated_node == nodes,6]),
      split = ""))
    colnames(stable_state) = c("Mutated")
    
    data = cbind(stable_state, WT_stable_state)
    
    data$Test = as.numeric(data$Mutated) == as.numeric(data$WT)
    
    colnames(stable_state_select) = colnames(data[,c(1,3)])
    
    stable_state_select = rbind(stable_state_select, data[,c(1,3)][data$Test == FALSE,])
  }
  
  return(unique(stable_state_select[complete.cases(stable_state_select),]))
}

#select which nodes change in combinations and single drugs
stable_state_select = extract_nodes(mutated_stable_state, WT_stable_state, combination = comb)

stable_state_select$Activity = c("changed")

#remove influential nodes from nodes merely changing activity
for(i in stable_state_select$Nodes){
  stable_state_select_V2 = stable_state_select[!grepl(i, subv$Nodes),]
}

subv = rbind(subv, stable_state_select_V2[,c("Nodes", "Activity")])

#generate new subgraph with selected nodes
#generate subgraph
g3 =  induced.subgraph(graph=g1,vids=subv$Nodes)

#remove nodes with no outgoing edges
remove_sink_nodes = function(g3){
  subv_3 = subv[subv$Nodes != "Antisurvival" & subv$Nodes != "Prosurvival",]
  g4 = induced.subgraph(graph=g1,vids=subv_3$Nodes)
  for(i in 1:length(V(g3))){
      if(any(degree(g4, mode = "out") == 0)){
    g3 = delete.vertices(g3, 
                V(g3)[ degree(g3, mode = "out") == 0 &
                                !grepl("Antisurvival", V(g3)$name) &
                                !grepl("Prosurvival", V(g3)$name) ])
  }
  
    else{
    
      }
  }

  
  return(g3)
  
}

g3 = remove_sink_nodes(g3)

subv = subv[match(vertex_attr(g3)$name, subv$Nodes),]

vertex_attr(g3)$type = subv$Activity
#color nodes according to beeing influential
vertex_attr(g3)$color[vertex_attr(g3)$type == "1"] ="green"
vertex_attr(g3)$color[vertex_attr(g3)$type == "0"] ="red"
vertex_attr(g3)$color[vertex_attr(g3)$type == "0.5"] ="orange"
vertex_attr(g3)$color[is.na(vertex_attr(g3)$type)] ="grey"
vertex_attr(g3)$color[vertex_attr(g3)$type == "changed"] = "lightblue"


#generate subgraph structure
plot(g3, edge.arrow.size=.2, 
     vertex.size=15, 
     vertex.frame.color="black", vertex.label.color="black", 
     vertex.label.cex=0.8, edge.curved=0.2,
     edge.color=ifelse(edge_attr(g3)$d =="activate", "green", "red"),
     layout=layout.davidson.harel)


