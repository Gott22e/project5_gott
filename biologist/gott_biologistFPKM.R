library(readr)
library(plyr)
library(ggplot2)

#Focal Genes for graph
sarg <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab")
mitg <- c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh") #Mpc1 not found?
cecg <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "Cdc27", "E2f1", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23")#Bora not found

#Programmer data
genexp <- read_table2("/projectnb/bf528/users/dachshund/project5_gott/programmer/P0_1_cufflinks/genes.fpkm_tracking")
cuffDiff <- read_table2("/projectnb/bf528/users/dachshund/project5_gott/programmer/cuffdiff_out/gene_exp.diff")

#Isolate necessary data from tables
#for comparing FPKM values
buildTable2 <- genexp[c("tracking_id","locus", "FPKM")] #to start the aggregate table
n_id <- unique(genexp[c("gene_short_name", "tracking_id", "locus")]) #to add the gene short name later
names(buildTable2) <- c("tracking_id","locus", "P0_1") #rename column names

ut <- data.frame(table(temp$locus))
ut <- ut[ut$Freq >1,]
colnames(ut) <- c("locus", "Freq")
gt <- genexp[genexp$locus %in% ut$locus,]

#Get top 1000 DE genes
cuffDiff <- cuffDiff[order(cuffDiff$q_value),]
top1k <- cuffDiff[1:1000,c("locus")]

#test <- read_table2("/project/bf528/project_2/data/fpkm_matrix.csv")

#get other samples genes
setwd("/project/bf528/project_2/data/samples")
files <- list.dirs(recursive=F)

#open files and append desired data to buildTable
for(name in files){
  name <- file[2]
  fname <- paste(substring(name, 3), "genes.fpkm_tracking", sep="/")
  if (file.exists(fname)){
    temp <- read_table2(fname)
    temp <- temp[c("tracking_id","locus", "FPKM")]
    names(temp) <- c("tracking_id","locus", substring(name, 3))
    buildTable2 <- merge(buildTable, temp, all.x = T, all.y = T, on= c("tracking_id","locus"))
  }
}


aggregateTable2 <- merge(buildTable2, n_id, by =c("tracking_id","locus") )[-1]
t1k <- aggregateTable[aggregateTable$locus %in% top1k$locus,]

#clean aggregate
aggregateTable <- aggregateTable[c(1:9)]

#clean t1k
t1k <- na.omit(merge(top1k, t1k, all.x=T))
t1k <- t1k[t1k$gene_short_name != "Scn4a",]
rownames(t1k) <- t1k$gene_short_name
numMatrix <- as.matrix(t1k[-c(1, 10)])

transform <- function(data){
  filtered <- aggregateTable[aggregateTable$gene_short_name %in% data,]
  filtered <- reshape(filtered, idvar= "gene", varying =colnames(filtered)[-9], 
                  v.name="FPKM", times=colnames(filtered)[-9], direction="long")
  names(filtered) <- c("Gene", "Group", "FPKM")
  return(filtered)}

#change directory where files will be saved
setwd("/projectnb/bf528/users/dachshund/project5_gott/biologist")

#plots
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

#make heatmap
heatmap(numMatrix)


