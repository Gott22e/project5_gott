library(readr)
#directory to get file(s) from
setwd("/projectnb/bf528/users/dachshund/project5_gott/programmer/cuffdiff_out/")
cuffDiff <- read_table2("gene_exp.diff") #import file

#change directory where files will be saved
setwd("/projectnb/bf528/users/dachshund/project5_gott/analyst")

#Order the table by q-value
cuffDiff <- cuffDiff[order(cuffDiff$q_value),]
#Top 10, only the requested values
top10 <- cuffDiff[1:10, c("gene", "value_1", "value_2", "log2(fold_change)", "p_value", "q_value")]

#histogram of distribution of log2(fold_change) for all genes
hist(cuffDiff$`log2(fold_change)`, main = "Histogram of Log2(fold_change)", 
     xlab = "Log2(fold_change)", col = c('purple', 'violet'))

#subset significant genes
sig <- subset(cuffDiff, cuffDiff$significant == "yes")
#histogram of log2(fold_change) for significant genes only
hist(sig$`log2(fold_change)`, main = "Significant Log2(fold_change) Frequency", 
     xlab = "Log2(fold_change)", col = c('red', 'pink'))

#subset upregulated genes
sigUp <- subset(sig, sig$`log2(fold_change)` > 0)
#subset downregulated genes
sigLow <- subset(sig, sig$`log2(fold_change)` < 0)

#write gene ids to file
write.table(sigUp$gene_id, "upRegulatedGeneID.csv", col.names = F, row.names = F, quote= F)
write.table(sigLow$gene_id, "downRegulatedGeneID.csv", col.names = F, row.names = F, quote= F)

#write top10 table to file
write.table(top10, "top10Table.csv", row.names = F, quote= F)

#subset significant genes
sigp <- subset(cuffDiff, cuffDiff$p_value < 0.01)
dim(sigp) #793 genes
#histogram of log2(fold_change) for significant genes only
hist(sigp$`log2(fold_change)`, main = "Significant Log2(fold_change) Frequency", 
     xlab = "Log2(fold_change)", col = c('blue', 'skyblue'))

#subset upregulated genes
sigpUp <- subset(sigp, sigp$`log2(fold_change)` > 0)
dim(sigpUp)
#subset downregulated genes
sigpLow <- subset(sigp, sigp$`log2(fold_change)` < 0)
dim(sigpLow)

