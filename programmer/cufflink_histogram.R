#Sets working directory; Reads table
setwd("/projectnb/bf528/users/dachshund/project5_gott/programmer/")
cd <- read.table("P0_1_cufflinks/genes.fpkm_tracking")

colnames(cd) <- cd[1,]#Set column names 
cd <- cd [-1,] #delete row of column names

table <- cd[cd$FPKM != 0,]
#Prepare the data

table$FPKM <- as.numeric(table$FPKM) #Change to numeric
#Create numeric lists; to.eX means to 10*X power.

#create pdf with histograms
pdf('cufflink.density.pdf')

hist(log(table[table$FPKM > 1,]$FPKM), main = "log(FPKM) where FPKM > 1", xlab="log(FPKM)", col = 'green')
hist(10**(table[table$FPKM <= 1,]$FPKM), main = "10^(FPKM) where FPKM <= 1", xlab="10^FPKM", col = 'blue')

dev.off()
