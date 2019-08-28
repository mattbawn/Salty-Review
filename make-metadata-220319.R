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
patric.nums <- readLines("patric.txt")


patric.sgi4 <- read.table("Ty_SGI4.csv", 
                          stringsAsFactors = FALSE, sep = ",", header = TRUE)
patric.sgi4$name <- gsub("ALL_REPORTS/report_", "mysnps_patric_", patric.sgi4$name)
patric.sgi4$name <- (gsub(".tsv", "", patric.sgi4$name))

R134.sgi4 <- read.table("R134_SGI4.csv", 
                          stringsAsFactors = FALSE, sep = ",", header = TRUE)
R134.sgi4$name <- gsub("ALL_REPORTS/report_", "mysnps_Rep_", R134.sgi4$name)
R134.sgi4$name <- (gsub(".tsv", "", R134.sgi4$name))

BII.sgi4 <- read.table("mather_SGI4.csv", 
                        stringsAsFactors = FALSE, sep = ",", header = TRUE)
BII.sgi4$name <- gsub("ALL_REPORTS/report_", "mysnps_BII_", BII.sgi4$name)
BII.sgi4$name <- (gsub(".tsv", "", BII.sgi4$name))




patric.nums <- gsub(".fna", "", patric.nums)
patric.nums <- t(as.data.frame(strsplit(patric.nums, "_")))
patric.nums <- as.data.frame(patric.nums)
rownames(patric.nums) <- NULL
colnames(patric.nums) <- c("strain", "number")
patric.nums$strain <- as.character(patric.nums$strain)

patric.name.data <- patric.name.data[patric.name.data$Genome.ID %in% patric.nums$strain,]
patric.merged <- merge(patric.name.data, patric.nums, by.x = "Genome.ID", by.y = "strain")

BII.reads <- readLines("BII_ERR.txt")
DT160.reads <- readLines("DT160_ERR.txt")

DT160.number <- lapply(strsplit(DT160.reads, "_"), '[', 3)
BII.number <- lapply(strsplit(BII.reads, "_"), '[', 3)

DT160.acession <- lapply(strsplit(DT160.reads, "_"), '[', 1)
BII.acession <- lapply(strsplit(BII.reads, "_"), '[', 1)

BII.data <- data.frame(cbind(BII.number, BII.acession))
DT160.data <- data.frame(cbind(DT160.number, DT160.acession))

BII.data$BII.number <- gsub(".fastq.gz", "", BII.data$BII.number)
DT160.data$DT160.number <- gsub(".fastq.gz", "", DT160.data$DT160.number)


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

par(bg="transparent")
p <- ggtree(my.tree, size=2, layout="circular")
# p <- ggtree(my.tree, size=2)

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

p1 <- gheatmap(p, metadata,  offset = 0.004, width=0.5, colnames = F, color = NULL)

p2 <- p %<+% mlst.all + geom_tiplab(aes(label=Strain, color=ST, angle=angle), size=2, align =TRUE , show.legend = TRUE) + geom_tippoint(aes(shape=Status, color=Status),size=5)

p3 <- p2 +theme(legend.position="right")



meta.patric <- grep("patric", metadata$V1)
meta.R134 <- grep("Rep", metadata$V1)
meta.BII <- grep("BII", metadata$V1)
meta.DT160 <- grep("DT160", metadata$V1)



new.metadata <- metadata
new.metadata$strain <- NA
new.metadata$acession <- NA
new.metadata$ID <- NA
new.metadata$ST <- NA
new.metadata$V1 <- gsub("mysnps_", "", new.metadata$V1)

# new.metadata$strain[meta.R134] <- r134.name.data$strain.x[gsub("Rep_", "", new.metadata$V1[meta.R134]) %in%  r134.name.data$V2]

new.metadata$acession[meta.R134] <- r134.name.data$acession[match(as.character(gsub("Rep_", "", new.metadata$V1[meta.R134])), as.character(r134.name.data$V2))]
new.metadata$strain[meta.R134] <- r134.name.data$strain.x[match(as.character(gsub("Rep_", "", new.metadata$V1[meta.R134])), as.character(r134.name.data$V2))]
new.metadata$acession[meta.patric] <- patric.merged$Assembly.Accession[match(as.character(gsub("patric_", "", new.metadata$V1[meta.patric])), as.character(patric.merged$number))]
new.metadata$strain[meta.patric] <- patric.merged$Strain[match(as.character(gsub("patric_", "", new.metadata$V1[meta.patric])), as.character(patric.merged$number))]
new.metadata$strain[meta.patric] <- patric.merged$Strain[match(as.character(gsub("patric_", "", new.metadata$V1[meta.patric])), as.character(patric.merged$number))]
new.metadata$ID[meta.patric] <- patric.merged$Genome.ID[match(as.character(gsub("patric_", "", new.metadata$V1[meta.patric])), as.character(patric.merged$number))]
new.metadata$ST[meta.BII] <- BII.mlst$V2[match(as.character(gsub("BII_", "", new.metadata$V1[meta.BII])), as.character(gsub("mysnps_BII_", "", BII.mlst$V1)))]
new.metadata$acession[meta.BII] <- BII.data$BII.acession[match(as.character(gsub("BII_", "", new.metadata$V1[meta.BII])), as.character(BII.data$BII.number))]

new.metadata$acession[meta.DT160] <- DT160.data$DT160.acession[match(as.character(gsub("DT160_", "", new.metadata$V1[meta.DT160])), as.character(DT160.data$DT160.number))]
new.metadata$ST[meta.DT160] <- DT160.mlst$V2[match(as.character(gsub("DT160_", "", new.metadata$V1[meta.DT160])), as.character(gsub("mysnps_DT160_", "", DT160.mlst$V1)))]
new.metadata$ST[meta.patric] <- patric.mlst$V2[match(as.character(gsub("patric_", "", new.metadata$V1[meta.patric])), as.character(gsub("mysnps_patric_", "", patric.mlst$V1)))]
new.metadata$ST[meta.R134] <- r134.mlst$V2[match(as.character(gsub("Rep_", "", new.metadata$V1[meta.R134])), as.character(gsub("mysnps_Rep_", "", r134.mlst$V1)))]

colnames(new.metadata) <- c("snippy-label", "source", "year", "strain", "acession", "patric-ID", "ST" )
new.metadata$acession <- unlist(new.metadata$acession)
write.table(new.metadata, "Figure_3_metadata.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

# new.metadata <- merge(new.metadata, patric.nums, by.x = "number", by.y = "number", all = TRUE )
# new.metadata$number[meta.R134] <- gsub("mysnps_Rep_", "", new.metadata$V1[meta.R134])
# new.metadata[meta.R134,] <- merge(new.metadata[meta.R134,], r134.name.data, by.x = "number", by.y = "V2", all = F )

# ggsave("Typhimurium_mlst-circle.pdf", width = 1000, height = 1000, units = "mm",
#         dpi = 600, bg = "transparent", limitsize = FALSE)
# rownames(metadata) <- metadata[,1]
# metadata <- metadata[, -1]
# 
# metadata$my.null[metadata$my.null == "turkey" | metadata$my.null == "chicken" ] <- "poultry"
# metadata$my.null[metadata$my.null == "feed" | metadata$my.null == "environment" ] <- "environment"
# metadata$my.null[metadata$my.null == "dog" | metadata$my.null == "cat" |  metadata$my.null == "sheep" |  metadata$my.null == "horse" ] <- "other"
# metadata$my.null[metadata$my.null == "pigeon" | metadata$my.null == "bird" ] <- "wild bird"
# 
# metadata[,2] <- NA
# 
# metadata$my.null[grep("mysnps_BII",rownames(metadata))] <- "pig"
# metadata$my.null[grep("mysnps_DT160",rownames(metadata))] <- "cattle"
# 
# patric.pig <- patric.merged$number[grep("sus",patric.merged$Host.Name)]
# 
# patric.pig <- paste0("mysnps_patric_", patric.pig)
# metadata$my.null[which(rownames(metadata) %in% patric.pig)] <- "pig"
# 
# patric.cattle <- patric.merged$number[grep("taurus",patric.merged$Host.Name)]
# 
# patric.cattle <- paste0("mysnps_patric_", patric.cattle)
# metadata$my.null[which(rownames(metadata) %in% patric.cattle)] <- "cattle"
# 
# metadata.2 <- as.data.frame(cbind(metadata, mlst.all$ST))
# 
# metadata.2$`mlst.all$ST`[metadata.2$`mlst.all$ST` == "ND"] <- NA
# metadata.2$`mlst.all$ST`[metadata.2$`mlst.all$ST` == "Novel*"] <- NA


metadata.2 <- as.data.frame(rbind(patric.sgi4, R134.sgi4, BII.sgi4))

meta.rownames <- metadata.2$name

metadata.2 <- metadata.2[,-1, drop=FALSE]

rownames(metadata.2) <- meta.rownames

sgi4.genes <- read.table("SGI4.clusters.tsv", sep = "\t", header = F)

colnames(metadata.2) <- gsub(".match", "", colnames(metadata.2))

colnames(metadata.2) <- sgi4.genes$V2[sgi4.genes$V1 %in% colnames(metadata.2)]

metadata.2 <- metadata.2[, order(colnames(metadata.2))]

metadata.2[metadata.2 == "no"] <- NA
# metadata.2[metadata.2 == "yes"] <- 1


metadata.2$SGI4_00066[metadata.2$SGI4_00066 == "yes"] <- "silE"
metadata.2$SGI4_00067[metadata.2$SGI4_00067 == "yes"] <- "cusS"
metadata.2$SGI4_00068[metadata.2$SGI4_00068 == "yes"] <- "cusR"
metadata.2$SGI4_00069[metadata.2$SGI4_00069 == "yes"] <- "cusC"
metadata.2$SGI4_00070[metadata.2$SGI4_00070 == "yes"] <- "cusF"
metadata.2$SGI4_00071[metadata.2$SGI4_00071 == "yes"] <- "cusB"
metadata.2$SGI4_00072[metadata.2$SGI4_00072 == "yes"] <- "cusA"
metadata.2$SGI4_00074[metadata.2$SGI4_00074 == "yes"] <- "silP"
metadata.2$SGI4_00077[metadata.2$SGI4_00077 == "yes"] <- "silE"
metadata.2$SGI4_00078[metadata.2$SGI4_00078 == "yes"] <- "copA"
metadata.2$SGI4_00079[metadata.2$SGI4_00079 == "yes"] <- "copB"
metadata.2$SGI4_00080[metadata.2$SGI4_00080 == "yes"] <- "copC"
metadata.2$SGI4_00081[metadata.2$SGI4_00081 == "yes"] <- "copD"
metadata.2$SGI4_00082[metadata.2$SGI4_00082 == "yes"] <- "copR"
metadata.2$SGI4_00083[metadata.2$SGI4_00083 == "yes"] <- "cusS"

cols.1=unlist(metadata.2)

cols.1[cols.1 == "yes"] <- "lightgrey"
cols.1[cols.1 == "silE"] <- "orange"
cols.1[cols.1 == "silP"] <- "orange"
cols.1[cols.1 == "cusS"] <- "forestgreen"
cols.1[cols.1 == "cusR"] <- "forestgreen"
cols.1[cols.1 == "cusC"] <- "forestgreen"
cols.1[cols.1 == "cusF"] <- "forestgreen"
cols.1[cols.1 == "cusB"] <- "forestgreen"
cols.1[cols.1 == "cusA"] <- "forestgreen"
cols.1[cols.1 == "copA"] <- "blue"
cols.1[cols.1 == "copB"] <- "blue"
cols.1[cols.1 == "copC"] <- "blue"
cols.1[cols.1 == "copD"] <- "blue"
cols.1[cols.1 == "copR"] <- "blue"
cols.1[cols.1 == "copS"] <- "blue"

names(cols.1) <- unlist(metadata.2)

p4 <- p2 %>% gheatmap(metadata.2,  offset = 0.02, width=0.2, colnames = F) %>%
  scale_x_ggtree()  + scale_fill_manual(values=cols.1)
# + scale_fill_manual(values=c(cols1, cols3))
# ggsave("Typhimurium_mlst-circle-smaller.pdf", width = 1000, height = 1000, units = "mm",
#         dpi = 600, bg = "transparent", limitsize = FALSE)