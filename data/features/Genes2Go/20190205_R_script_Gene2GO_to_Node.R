#Genes to Nodes for GO terms
library(tidyr)
setwd("M:/PhD/projects/Influential Node Project/influentialnodes/data/cascade/Features/Genes2Go")

data = read.csv("result_analysis_selected_GO_terms.txt", sep = "\t", skip = 2, header = FALSE, stringsAsFactors = F)

nodes = read.csv("20190130_nodes_to_HGNC.txt", sep = "\t", stringsAsFactors = F)

#merge GO ID with term
GO_terms = as.data.frame(list(data[1,]))
GO_terms[2,] = as.data.frame(list(data[2,]))
GO_terms = as.data.frame(t(GO_terms))

GO_terms[,3] = paste(GO_terms[,1], GO_terms[,2], sep="-")

colnames(data) = GO_terms[,3]
data = data[-2,]
data = data[-1,]

colnames(data)[colnames(data)=="GO_Ids-Gene_Ids/GO_Terms"] = "HGNC.symbols"

rownames(nodes) = nodes[,1]
res = data.frame(Node= character(0), HGNC.symbols = character(0), stringsAsFactors = F)
for (node in nodes[,1]) {
  hgnc.symbols = unlist(strsplit(as.character(nodes[node, ][2]), split = ", "))
  # print(paste0(node, " and symbols ", hgnc.symbols))
  
  for (hgnc.symbol in hgnc.symbols) {
    res = rbind(res, c(node, hgnc.symbol), stringsAsFactors = F)
  }
}

colnames(res) = c("Nodes", "HGNC.symbols")

res = merge(res, data[1:197,], by = "HGNC.symbols", all.y = TRUE)

#remove HGNC symbols and keep only genes
test = aggregate(res[,2:21], by = list(res$Nodes), max)
test = test[,-1]

write.table(test, "20190205_Gene2GO_node_list.txt", quote = FALSE, sep = "\t", row.names = F)




