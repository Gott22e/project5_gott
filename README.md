# Gott Project 5 #

This project attempted to replicate some of the studies by O'Meara et al by identifying the differentially expressed genes between newborn and adult mouse myocytes. This repository contains the scripts used to analyse the data cited below

O’Meara, C. C., Wamstad, J. A., Gladstone, R. A., Fomovsky, G. M., Butty, V. L., Shrikumar, A., Gannon, J. B., Boyer, L. A., & Lee, R. T. (2015). Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration. Circulation Research, 116(5), 804–815. https://doi.org/10.1161/CIRCRESAHA.116.304269

### Data Curator ###
fastq.qsub - 
Qsub script to split the raw data and generate qc metrics

### Programmer ###
cufflink_histogram.R Requires R 4.0.2
Creates a histogram of the cufflinks FPKM output

p_hits.qsub Requires Python 3.7.9, Samtools 0.1.19, and RSeQC 3.0.0
Runs geneBody_coverage.oy, inner_distance.py, and bam_stat.py on the Tophat accepted_hits.bam output

run_cuffdiff.qsub Requires Cufflinks 2.2.1
Runs cuffdiff from Cufflinks on the Tophat output for P0_2 and P0_2 vs Ad_1 and Ad_2 to determine differentially expressed genes

run_cufflinks.qsub Requires Cufflinks 2.2.1
Runs Cufflinks on the Tophat output accepted_hits.bam. Quantifies FPKM values

run_tophat.qsub- Requires Python 2.7.16, Samtools 0.1.19, Bowtie2 2.4.2 , and Boost 1.69.0
Runs Tophat on the .fastq files produced from fastq.qsub. Aligns the reads


### Analyst ###
gott_p2_analysisFIltering.R Requires R 4.0.2 and package 'readr' 1.4.0
Determines significant genes from the Cuffdiff output, plots histograms of the values, and creates lists of the significant up-regulated and down-regulated genes


### Biologist ###
gott_biologistFPKM.R Requires R 4.0.2 and packages 'readr' 1.4.0, 'dplyr' 0.7.6, and 'ggplot2' 1.4.0
Creates tables from the output of Cufflinks across all 8 samples. Selects the top 1000 genes differentially expressed genes from the Cuffdiff output to make a heatmap and narrows a different table to to plot specific gene FPKM frequencies across the samples
