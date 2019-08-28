options(stringsAsFactors = FALSE)
options(digits=10)
library(RColorBrewer)
library(colorspace)

r134.name.data <- read.table("samples-overview-representative_host.csv", 
                             stringsAsFactors = FALSE, sep = ",", header = TRUE)


patric.name.data <- read.table("PATRIC_genome.csv", 
                               colClasses ="character",  
                            stringsAsFactors = FALSE, sep = ",", header = TRUE,
                            comment.char = "")

patric.nums <- readLines("patric.txt")


r134.mlst <- read.table("R134_mlst_summary.txt", 
                        stringsAsFactors = FALSE, sep = " ", header = F)


patric.mlst <- read.table("patric_mlst_summary.txt", 
                          stringsAsFactors = FALSE, sep = " ", header = F)

DT160.mlst <- read.table("DT160_mlst_summary.txt", 
                         stringsAsFactors = FALSE, sep = " ", header = F)

BII.mlst <- read.table("BII_mlst_summary.txt", 
                       stringsAsFactors = FALSE, sep = " ", header = F)

r134.mlst$V1 <- gsub("report_mlst", "mysnps_Rep", r134.mlst$V1)
patric.mlst$V1 <- gsub("mlst_report", "mysnps_patric", patric.mlst$V1)
DT160.mlst$V1 <- gsub("mlst_report", "mysnps_DT160", DT160.mlst$V1)
BII.mlst$V1 <- gsub("mlst_report", "mysnps_BII", BII.mlst$V1)

r134.mlst$V3 <- r134.name.data$strain.y[as.numeric(lapply(strsplit(r134.mlst$V1, "_"), '[',3))]
patric.mlst$V3 <- patric.name.data$Strain[as.numeric(lapply(strsplit(patric.mlst$V1, "_"), '[',3))]
DT160.mlst$V3 <- "NZ DT160"
BII.mlst$V3 <- "North America BII"


mlst.all <- rbind(r134.mlst, patric.mlst, DT160.mlst, BII.mlst)

colnames(mlst.all) <- c("tip", "ST", "Strain")

patric.nums <- gsub(".fna", "", patric.nums)

patric.nums <- t(as.data.frame(strsplit(patric.nums, "_")))

patric.nums <- as.data.frame(patric.nums)

rownames(patric.nums) <- NULL

colnames(patric.nums) <- c("strain", "number")

patric.nums$strain <- as.character(patric.nums$strain)

patric.name.data <- patric.name.data[patric.name.data$Genome.ID %in% patric.nums$strain,]

patric.merged <- merge(patric.name.data, patric.nums, by.x = "Genome.ID", by.y = "strain")

library(ggtree)
library(phangorn)

my.tree <- read.tree("RAxML_bestTree.raxml-pat_R134_DT160_BII")

my.origional.tree <- my.tree

# Remove Reference from tree

my.tree <- drop.tip(my.tree, "Reference")

my.tips <- my.tree$tip.label

old.tips <- my.tips

r134 <- grep("mysnps_Rep", my.tips)

patric <- grep("mysnps_patric", my.tips)

my.tips[r134] <- sapply(my.tips[r134], function(x) as.numeric(gsub("mysnps_Rep_","", x)))

my.tips[patric] <- sapply(my.tips[patric], function(x) as.numeric(gsub("mysnps_patric_","", x)))

my.tips[r134] <- r134.name.data$V2[as.numeric(my.tips[r134])]

my.tips[patric] <- patric.merged$number[as.numeric(my.tips[patric])]

my.tips[r134] <- r134.name.data$strain.y[as.numeric(my.tips[r134])]

my.tips[patric] <- patric.merged$Strain[as.numeric(my.tips[patric])]

my.tree$tip.label <- old.tips

my.tree <- midpoint(my.tree)


complete.patric <- patric.merged[patric.merged$Sequencing.Status == "complete",]

test <- paste0("mysnps_patric_",complete.patric$number)

mlst.all$comp <- NA

# mlst.all$tip[! mlst.all$tip %in% test] <- ""



par(bg="transparent")
p <- ggtree(my.tree, size=2, layout="circular")

r134.name.data$source[r134.name.data$source == "cattle "] <- "cattle"
r134.name.data$source[r134.name.data$source == "pig "] <- "pig"
r134.name.data$source[r134.name.data$source == "duck "] <- "duck"

# cols2= c(brewer.pal(15, "Paired"))
cols1=cols2= c(brewer.pal(8, "Accent"))
cols2=rainbow_hcl(22)
cols3=colorRampPalette(brewer.pal(8, "Accent"))(22)
cols4=colorRampPalette(brewer.pal(8, "Blues"))(18)
# tree.tips <- my.tree$tip.label




source.data <- as.data.frame(cbind(r134.name.data$V2, r134.name.data$source, r134.name.data$Year))

source.data$V1 <- paste0("mysnps_Rep_", source.data$V1)

my.null <- rep(NA, length(mlst.all$tip))
my.year <- rep(NA, length(mlst.all$tip))

metadata <- as.data.frame(cbind(mlst.all$tip, my.null, my.year))

metadata$my.null[match(source.data$V1, metadata$V1)] <- source.data$V2
metadata$my.year[match(source.data$V1, metadata$V1)] <- source.data$V3
rownames(source.data) <- source.data$V1

source.data <- source.data[,-1]

mlst.all$Status <- NA

mlst.all$Status[mlst.all$tip %in% paste0("mysnps_patric_", patric.merged$number)] <- patric.merged$Sequencing.Status

p1 <- gheatmap(p, metadata,  offset = 0.002, width=0.2, colnames = F)

p2 <- p %<+% mlst.all + geom_tiplab(aes(label=Strain, color=ST, angle=angle), size=4, align =TRUE , show.legend = TRUE) + geom_tippoint(aes(shape=Status, color=Status),size=5)

p3 <- p2 +theme(legend.position="right")
# ggsave("Typhimurium_mlst-circle.pdf", width = 1000, height = 1000, units = "mm",
#         dpi = 600, bg = "transparent", limitsize = FALSE)
rownames(metadata) <- metadata[,1]
metadata <- metadata[, -1]

metadata$my.null[metadata$my.null == "turkey" | metadata$my.null == "chicken" ] <- "poultry"
metadata$my.null[metadata$my.null == "feed" | metadata$my.null == "environment" ] <- "environment"
metadata$my.null[metadata$my.null == "dog" | metadata$my.null == "cat" |  metadata$my.null == "sheep" |  metadata$my.null == "horse" ] <- "other"
metadata$my.null[metadata$my.null == "pigeon" | metadata$my.null == "bird" ] <- "wild bird"

metadata[,2] <- NA

metadata$my.null[grep("mysnps_BII",rownames(metadata))] <- "pig"
metadata$my.null[grep("mysnps_DT160",rownames(metadata))] <- "cattle"

patric.pig <- patric.merged$number[grep("sus",patric.merged$Host.Name)]

patric.pig <- paste0("mysnps_patric_", patric.pig)
metadata$my.null[which(rownames(metadata) %in% patric.pig)] <- "pig"

patric.cattle <- patric.merged$number[grep("taurus",patric.merged$Host.Name)]

patric.cattle <- paste0("mysnps_patric_", patric.cattle)
metadata$my.null[which(rownames(metadata) %in% patric.cattle)] <- "cattle"

metadata.2 <- as.data.frame(cbind(metadata, mlst.all$ST))

metadata.2$`mlst.all$ST`[metadata.2$`mlst.all$ST` == "ND"] <- NA
metadata.2$`mlst.all$ST`[metadata.2$`mlst.all$ST` == "Novel*"] <- NA


p4 <- p2 %>% gheatmap(metadata.2,  offset = 0.01, width=0.2, colnames = F) +
  scale_fill_manual(values=c(cols1, cols3))
ggsave("Typhimurium_mlst-circle-smaller.pdf", width = 1000, height = 1000, units = "mm",
        dpi = 600, bg = "transparent", limitsize = FALSE)