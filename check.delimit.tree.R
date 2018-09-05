#Rscript to check a delimitation result.
library(ape)

tr <- read.tree("./sistrurus.delimit.tre")

tr$node.label <- substr(tr$node.label, 1,6)
plot(tr, show.node.label=T)

