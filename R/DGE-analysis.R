#Script for loading and preprocessing RNAseq-data + 
#performing differential expression analysis between
#low and high survival patients

install.packages('dplyr')
library(dplyr)
library(tibble)
info_GBM=read.delim("infoSummary_clinical_GBM.txt",stringsAsFactors = F)
counts <- read.csv('RNAdataSummary_count_GBM.txt', stringsAsFactors = F, check.names = F, sep='\t')
N=nrow(counts)
M=nrow(GBM_data)
counts_unique = as.data.frame(counts[1:19571,])
O=1
for (j in 8:M){
  for (i in 1:N){
    if (as.character(GBM_data[j,1]) == as.character(counts[i,1])){
      counts_unique[O,] = counts[i,]
      O = O+1
      print(O)
      break
    }
  }
}
# make the column names of the two files compatible
colnames(counts_unique) = gsub('\\..*', '', colnames(counts_unique))
counts_ds$caseSubmitterID = gsub('*-11R-1850-01', '', counts_ds$caseSubmitterID)
info_GBM = info_GBM %>% mutate(new= gsub('\\..*', '', fileName))
all(info_GBM$new %in% colnames(counts_unique))
setdiff(colnames(counts_unique), info_GBM$new)

#remove all rows with None in fileName
info_GBM = info_GBM[info_GBM$fileName != 'None',]

#make sure ensembl id sticks
rownames(counts_unique) = counts_unique[,1]
counts_unique = counts_unique[,-1]

#merge count-data and clinical information to one data frame
counts2 <- as.data.frame(t(counts_unique))
counts2$new <- rownames(counts2)
merged <- merge.data.frame(info_GBM, counts2, by='new')

### differential expression analysis ###
source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
library('DESeq2')
biocLite("biomaRt")
library(biomaRt)

#fixing a column named Living Days
N=nrow(merged)
for (i in 1:N){
  if (merged[i,12] == 'None'){
    merged[i,12] = merged[i,4]
  }
}
colnames(merged)[12]<-"LivingDays"
merged$LivingDays = as.numeric(merged$LivingDays)

#find out cutoff for low and high
quantile(merged$LivingDays, c(1/3,2/3)) #cutoff should be below 211 and above 454

#extract low and high only
low_and_high = merged[merged$LivingDays<211 | merged$LivingDays>454,]
#replace with 1 (high surival) and 0 (low survival)

M=nrow(low_and_high)
for (i in 1:M){
  if (low_and_high[i,12] > 454){
    low_and_high[i,12] = 1
  }
  else{
    low_and_high[i,12] = 0
  }
}
low_and_high$LivingDays <- as.factor(low_and_high$LivingDays)

#create condition vector
samples_ds = data.frame(low_and_high[,12])
colnames(samples_ds) = "LivingDays"
rownames(samples_ds) = samples_ds[,1]
samples_ds = samples_ds[,-(1)]

#create countData
counts_ds = low_and_high[,c(16,18:60505)]
rownames(counts_ds) = counts_ds[,1]
counts_ds = counts_ds[,-1]
counts_ds = as.data.frame(t(counts_ds))
ds=DESeqDataSetFromMatrix(countData=counts_ds_new_temp, colData=samples_ds_temp, design=~LivingDays)
ds=DESeq(ds)
res <- results(ds)

#removing genes with low counts across all samples
nrow(res)
sum( is.na(res$pvalue) )
res = res[ ! is.na(res$pvalue), ]
nrow(res)
res = res[ order(res$padj), ]

#annotation
ensembl = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",
dataset="hsapiens_gene_ensembl")
genemap = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
filters = "ensembl_gene_id", values = rownames(res), mart = ensembl )
symbols <- tapply(genemap$hgnc_symbol, genemap$ensembl_gene_id, paste, collapse="; ")
res$symbol <- symbols[ rownames(res) ]

#filter based on p-value
sig <- res[ which(res$pvalue < 0.05), ]

#separate into up/downregulated, order, convert, write to table
sig_up = sig[which(sig$log2FoldChange>0),]
sig_down = sig[which(sig$log2FoldChange<0),]
sig_up <- sig_up[ order(sig_up$padj), ]
sig_down <- sig_down[ order(sig_down$padj), ]
sig_up <- as.data.frame(sig_up)
sig_down <- as.data.frame(sig_down)

write.table(sig_up, "results_deseq_20180313_up.txt", sep="\t")
write.table(sig_down, "results_deseq_20180313_down.txt", sep="\t")
