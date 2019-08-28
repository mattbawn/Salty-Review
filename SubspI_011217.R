
library(phangorn)
library(ape)
library(ggtree)

options(stringsAsFactors = FALSE)

my.tree <- read.tree("RAxML_bestTree.SubspI_II")

my.tree <- drop.tip(my.tree, "Reference")

my.tips <- my.tree$tip.label

my.tips <- gsub("mysnps_SubspI_", "", my.tips) 
my.tips <- gsub("mysnps_patric_", "", my.tips) 

my.nos <- c(1:15,20,21,22,31,55) 

my.strain <- c("Agona", "Cholerasesuis", "Dublin", "Enteritidis",
               "Gallinarium", "Heidelberg", "Java", "Typhimurium",
               "Newport", "Paratyphi A", "Paratyphi B", "Paratyphi C",
               "Schwarzengrund", "Typhi", "Weltevreden", "D23580",
               "DT2", "14028", "SL1344", "LT2")


my.data <- as.data.frame(cbind(my.nos, my.strain))

my.tree <- midpoint(my.tree)

my.tree$tip.label <- my.data$my.strain[match(my.tips, my.data$my.nos)]

my.tree <- drop.tip(my.tree, "Typhimurium")

par(bg="transparent")
 p <- ggtree(my.tree,layout = "circular", size=2)

# p <- ggtree(my.tree, size=2)

# Tip labels have to be explicitly added
p <- p + geom_tiplab(size=8)
ggsave("Salmonella-tree-3.pdf", width = 2000, height = 2000, units = "mm",
       dpi = 600, bg = "transparent", limitsize = FALSE)
