library(readr)
library(plyr)
library(ggplot2)

#Focal Genes for graph
sarg <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab")
mitg <- c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
cecg <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23")

#Programmer data
genexp <- read_table2("/projectnb/bf528/users/dachshund/project5_gott/programmer/P0_1_cufflinks/genes.fpkm_tracking")
aggregateTable <- genexp[c("tracking_id",  "FPKM")] #to start the aggregate table
n_id <- genexp[c("gene_short_name", "tracking_id")] #to add the gene short name later
names(aggregateTable) <- c("tracking_id", "P0_1") #rename column names
#test <- read_table2("/project/bf528/project_2/data/fpkm_matrix.csv")

#get other samples genes
setwd("/project/bf528/project_2/data/samples")
files <- list.dirs(recursive=F)

for(name in files){
  fname <- paste(substring(name, 3), "genes.fpkm_tracking", sep="/")
  if (file.exists(fname)){
    temp <- read_table2(fname)
    temp <- temp[c("tracking_id", "FPKM")]
    names(temp) <- c("tracking_id", substring(name, 3))
    aggregateTable <- merge(aggregateTable, temp)
  }
}

aggregateTable <- merge(aggregateTable, n_id)[-1]


transform <- function(data){
  filtered <- aggregateTable[aggregateTable$gene_short_name %in% data,]
  filtered <- reshape(filtered, idvar= "gene_short_name", varying =colnames(filtered)[-(length(filtered))], 
                  v.name="FPKM", times=colnames(filtered)[-(length(filtered))], direction="long")
  names(filtered) <- c("Gene", "Group", "FPKM")
  return(filtered)
}

#change directory where files will be saved
setwd("/projectnb/bf528/users/dachshund/project5_gott/biologist")

ggplot(data = transform(sarg), aes(x=Gene, y=FPKM, group=Group, colour=Group)) +
  geom_line() +
  geom_point() +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1)) +
  ggtitle("Sarcomere")
ggplot(data = transform(mitg), aes(x=Gene, y=FPKM, group=Group, colour=Group)) +
  geom_line() +
  geom_point() +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1)) +
  ggtitle("Mitochondria")
ggplot(data = transform(cecg), aes(x=Gene, y=FPKM, group=Group, colour=Group)) +
  geom_line() +
  geom_point() +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Cycle")


