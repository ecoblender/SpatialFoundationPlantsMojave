
library(ggplot2)
graph <- read.csv("joincountoutput.csv")
ggplot(graph, aes(fill = Result, x = species, y = Value)) + geom_bar(position = "dodge",stat = "identity") + facet_grid(. ~Transect) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, size = 14, hjust = 1, vjust = 0.25), axis.text.y = element_text(size =14)) + xlab("Species") + scale_x_discrete(labels=c("A. salsola", "Ephedra sp.", "L. tridentata - \n A. salsola", "S.mexicana")) + ylab("Join Count")

windows()
