#R script to make a venn diagram of binarized influential node data per cell line
# install.packages('VennDiagram')

library(VennDiagram)

data = read.csv("../Classification_influential_nodes.txt", sep ="\t", stringsAsFactors = F)

#select data for diagram
data.venn = list()
data.venn$AGS = data$X[data$AGS == 1]
data.venn$`DU-145` = data$X[data$DU145 == 1]
data.venn$`COLO 205` = data$X[data$COLO205 == 1]
data.venn$`SW-620` = data$X[data$SW620 == 1]

#save venn to file
venn.diagram(data.venn, "Venn_diagramm_influential_nodes", height = 3000, width = 3000,
             resolution = 600, imagetype = "png", units = "px",
             lwd = 2, fill=rainbow(4), apha = 0.5, filename = NULL,
             cex = 0.8)


venn.plot = venn.diagram(data.venn, lwd = 2,
                         fill=rainbow(4), apha = 0.5, filename = NULL,
                         cex = 0.8) 
pdf(file="Venn_diagram_influential_nodes.pdf")
grid.draw(venn.plot)
dev.off()
