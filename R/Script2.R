#Analysis of gag Treestructure- 04/07/2022

getwd()
setwd("/Users/peggy/Desktop/R scripts for the PANGEA project/Gag Final Treestructure")
#LOAD THE FOLLOWING LIBRARIES
library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)
library(chisq.posthoc.test)
library(dplyr)
library(ape)
library(treedater)
library(treestructure)
library(msa)
library(Biostrings)
library(ggpubr)

geneRegion <- "gag"



# read in alignment
alignment_name <- "gag_finalAli_wCGR.fas"
data <- read.dna(alignment_name, format="fasta", as.matrix=FALSE)


# read in tree, add on dates
tr <- read.tree(paste("RAxML_bestTree.ML_BW_gag", ".nwk", sep=""))
tr2 <- root(tr, outgroup = tr$tip.label[grep("AY371157", tr$tip.label)])
seq <- read.dna(paste("gag_finalAli_wCGR.fas", sep=""), format="fasta", as.character = F, as.matrix=F)
sts <- as.numeric(unlist(lapply(tr2$tip.label, function(x) {strsplit(x, "_")[[1]][2]})))
print(missingDates <- tr2$tip.label[is.na(sts)])
###### check only the outgroup is missing date before running next line
sts[is.na(sts)] <- 2001.5
sts <- as.numeric(sts)
names(sts) <- tr2$tip.label
set.seed(123)
dt2 <- dater(tr2, sts, s = length(data[[1]]), ###grab the length of the first sequence
             numStartConditions = 1, maxit = 1, meanRateLimits = c(0.001,  0.004), ncpu = 8)

#looking for duplication of sequences in the phylogenies 
#cleanTree <-dt2[!duplicated(dt2$tip.label)]-------




# save the time tree
saveRDS(  dt2, file = paste0("editNames_BWonly_dated", geneRegion, "_rooted_dated.rds") )

dt2 <- readRDS(paste0("editNames_BWonly_dated", geneRegion, "_rooted_dated.rds"))

dt3 <- drop.tip(dt2, tr$tip.label[grep("AY371157", tr$tip.label)])

set.seed(123)
ts <- trestruct(dt3, minCladeSize = 1000, minOverlap = -Inf, nsim = 100, # increase nsim later
                level = 0.01, ncpu = 8, verbosity = 1, debugLevel = 1)
ts # need to check number of clusters and their sizes before going any further
saveRDS(  ts, file = paste0("best_tree", geneRegion, "_rooted_dated.rds") )
ts1 <- readRDS(paste0("best_tree", geneRegion, "_rooted_dated.rds"))
ts_df_2 <- as.data.frame(ts)

write.csv(ts_df_2, paste(geneRegion,"_treestructure_splits1.csv", sep=""), row.names=F)

ts_df_2 <- read.csv(paste(geneRegion,"_treestructure_splits1.csv", sep=""))

length(grep("CGR", ts_df_2$taxon))

for (i in 1:max(as.numeric(ts_df_2$partition))){
  print(i)
  mySeq <- ts_df_2$taxon[ts_df_2$partition==i] 
  print(length(mySeq))
  #myNewAli <- seq[match(mySeq, names(seq))]
  #write.FASTA(myNewAli, paste(geneRegion, "_partition", i, "_v1_1.fas", sep=""))
  write.table(mySeq, paste0("seq_partition_",i,".txt"), col.names = F, row.names = F, quote=F)
}


#realignment 

taxa <- attributes(data)$names 


for (i in 1:max(as.numeric(ts_df_2$partition))) { # a loop for each partition
  
  
  # list of names of sequences you want in the final alignment
  list_name <- paste0("seq_partition_", i, ".txt")
  list <- read.table(list_name, header=FALSE, as.is=TRUE)
  
  
  # name of alignment which will just contain the seqs from list
  outname <- paste0 ("seq_partition_", i,".fas")
  
  # choose between one or the other line below
  inds <- match(list[,1],taxa)# exact matches
  #inds <- unlist(lapply(list[,1], function(x) grep(x, taxa)))
  
  missingSeq <- list[which(is.na(inds)),]
  print("These sequences are missing:")
  print(missingSeq)
  
  inds2 <- inds[!is.na(inds)]
  
  # print(data[inds])
  
  write.dna(data[inds2], file=outname, format="fasta", nbcol=-1, colsep="")
}
