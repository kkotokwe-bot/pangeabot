setwd("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Treestructure Final Analysis_01122021")
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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")

geneRegion <- "pol"



mainepi_2 <- read.csv("BW_June2021_forKena_editBWnames_wVillage_edit.csv",TRUE, ",", stringsAsFactors = F)
allnames_2 <- read.csv("pol_treestructure_splits1.csv",TRUE, ",", stringsAsFactors = F)
csvfinal_2 <- merge(mainepi_2, allnames_2, by.x = "BWid",
                    by.y = "taxon", all.x = TRUE, all=FALSE)

#removing the sequences that dont have partitions
new_demographicscombination <- csvfinal_2[!is.na(csvfinal_2$partition),]
demogFinal0 <- new_demographicscombination %>%     group_by(BWid) %>%   filter(n()==1)

##Removing non- BCPP sequences(only mochudi in this cohort)
demogFinal1<-demogFinal0[-grep("mochudi", demogFinal0$study_FINAL),]
#demogFinal2 <-demogFinal1[-grep("bhp012", demogFinal1$study_FINAL),]
#demogFinal2 <-demogFinal1 [-grep("CGR", demogFinal1$BWid),]
#demogFinal2<- demogFinal1[-grep("CGR",attributes(demogFinal1)$BWid)]

#remove the duplicates IDs and reference sequences for analysis 
demogFinal2 <- demogFinal1[!is.na(demogFinal1$study_FINAL),]
demogFinal3 <-demogFinal2[!duplicated(demogFinal2$BWid),]



#Write the CSV file for the combined demographics 
write.csv(demogFinal3 , paste0("polfinal", "_file.csv",sep=""))

library (data.table)
DC <- data.table(demogFinal3 )
DC[, districts := ifelse(community %in% c("ranaka", "digawana","molapowabojang","mmankgodi", "mmathethe"), "southern",
                         ifelse(community %in% c("lerala", "sefophe", "maunatlala", "ramokgonami","mmadinare", "shoshong","nkange","sebina","mathangwane", "mmandunyane","gumare","rakops","sefhare","tsetsebjwe","nata"), "central",
                                ifelse(community %in% c("masunga", "tati_siding"), "north east",
                                       ifelse(community %in% c("otse"), "south east",
                                              ifelse(community %in% c("bokaa","oodi" ), "kgatleng",
                                                     ifelse(community %in% c("letlhakeng","lentsweletau","metsimotlhabe"), "kweneng",
                                                            ifelse(community %in% c("gweta", "shakawe"), "ngamiland", NA)))))))]
#write and read the csv file

write.csv(DC , paste("_combinedpoldemo1", "_file.csv",sep=""))
FINALDEMOFILE <- "DCcombinedpoldemo1file.csv"

#add the districts column 
demogFinal3$Districts <- DC$districts


#Statistical analysis with partitions

df1 <-data.frame(partitionID=0, F=0, M=0, sexNA=0, meanAge=0, central=0, kgatleng=0, kweneng=0, ngamiland=0, northeast=0, southeast=0, southern=0, DistrictsNA=0, chronic=0, recent=0,stagetruepredictNA=0, decimalDate=0, no=0, yes=0,ARTyesnoNA=0, stringsAsFactors = F)
for (i in 1:max(demogFinal3$partition)){
  print(i)
  partition_csv1 <- demogFinal3[demogFinal3$partition==i,]
  
  #sex table
  sex_tab <- round(prop.table(table(partition_csv1$genderFINAL, useNA= "always"))*100,1)
  
  #calculate mean age
  meanAgePartition <- round(mean(partition_csv1$enrollAgeFINAL, na.rm = T),1)
  
  #Districts
  DISTRICTS <- round(prop.table(table(partition_csv1$Districts, useNA = "always"))*100,1)
  
  #stageofinfection
  stagetruepredict<- round(prop.table(table (partition_csv1$stage_true_predict, useNA = "always"))*100,1)
  
  #Year
  decimalDate <- round(mean(partition_csv1$decimalDate, na.rm= T, useNA= "always"),1)
  
  #ARTstatus
  TreatmentStatus <- round(prop.table(table(partition_csv1$ARTyesno, useNA= "always"))*100,1)
  
  df1[i,] <- c(i, as.numeric(sex_tab[1]), as.numeric(sex_tab[2]), as.numeric(sex_tab[3]), meanAge=meanAgePartition, as.numeric(DISTRICTS[1]), as.numeric(DISTRICTS[2]), as.numeric(DISTRICTS[3]), as.numeric(DISTRICTS[4]), as.numeric(DISTRICTS[5]), as.numeric(DISTRICTS[6]), as.numeric(DISTRICTS[7]), as.numeric(DISTRICTS[8]), as.numeric(stagetruepredict[1]), as.numeric(stagetruepredict[2]), as.numeric(stagetruepredict[3]), decimalDate=decimalDate, as.numeric(TreatmentStatus[1]),as.numeric(TreatmentStatus[2]),as.numeric(TreatmentStatus[3]))
}
#add a column of total sequences per partition
df1$Totalsequences <-table (DC$partition)
write.csv(df1 , paste("df11", "_Finaldemographicspartition.csv",sep=""))


###Chi with partitions (Pearson's Chi-squared test)

#gender
table(demogFinal3$partition,demogFinal3$genderFINAL)						
chisq.test(table(demogFinal3$partition, demogFinal3$genderFINAL ))	

#Stage of infection
table(demogFinal3$partition,demogFinal3$stage_true_predict)						
chisq.test(table(demogFinal3$partition, demogFinal3$stage_true_predict ))	

#Treatment status
table(demogFinal3$partition, demogFinal3$ARTyesno)						
chisq.test(table(demogFinal3$partition, demogFinal3$ARTyesno ))	



###Further analysis of Chisq.posthoc.test

#stage of infection
chisq.resultsp <- chisq.test(table(demogFinal3$partition, demogFinal3$stage_true_predict ))
chisq.resultsp$stdres

chisq.posthoc.test(table(demogFinal3$partition, demogFinal3$stage_true_predict ),
                   method = "bonferroni")


#Treatment
chisq.resultsp <- chisq.test(table(demogFinal3$partition, demogFinal3$ARTyesno))
chisq.resultsp$stdres

chisq.posthoc.test(table(demogFinal3$partition, demogFinal3$ARTyesno ),
                   method = "bonferroni")



#How too read CGR in every partition
fas <- read.dna("seq_partition_1.fas", format="fasta")
length(grep("CGR",rownames(fas)))



#FINAL ANOVA and Tukey for Age, decimaldates, districts


#Analysis of age
library(dplyr)
group_by(DC, partition) %>%
  summarise(
    count = n(),
    mean = mean(enrollAgeFINAL, na.rm = TRUE, useNA = "always"),
    sd = sd(enrollAgeFINAL, na.rm = TRUE),
    median = median(enrollAgeFINAL, na.rm = TRUE),
    IQR = IQR(enrollAgeFINAL, na.rm = TRUE)
  )

install.packages("ggpubr")
library("ggpubr")

ggboxplot(DC, x = "partition", y = "enrollAgeFINAL",
          color = "partition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.5,
          order = c("2", "3", "1"),
          ylab = "MedianAge", xlab = "Partitions")
#same thing as above
ggboxplot(DC, x = "partition", y = "enrollAgeFINAL")
#decided to exclude as it is not statistically significant

#trying gglines
ggline(DC, x = "partition", y = "enrollAgeFINAL",
       add = c("mean_se", "jitter"),
       bxp.errorbar = TRUE, bxp.errorbar.width = 0.5,
       order = c("1", "2", "3"),
       ylab =  "MeanAge", xlab = "Partitions")
# Compute the analysis of variance for AGE
res.aov <- aov(enrollAgeFINAL ~ factor(partition), data = DC)
# Summary of the analysis
summary(res.aov)

#tukey analysis
TukeyHSD(res.aov)

#plot Tukey plot for AGE
plot.tukey1 <- TukeyHSD(res.aov)
plot(plot.tukey1)


#proving with other statistical tests
pwc <- DC %>% 
  dunn_test(enrollAgeFINAL ~ partition, p.adjust.method = "bonferroni") 
pwc
#proving with other statistical tests
pwc2 <-  DC %>% 
  wilcox_test(enrollAgeFINAL ~ partition, p.adjust.method = "bonferroni")
pwc2


#Decimaldates anova and tukey
group_by(DC, partition) %>%
  summarise(
    count = n(),
    mean = mean(decimalDate, na.rm = TRUE, useNA = "always"),
    sd = sd(decimalDate, na.rm = TRUE),
    median = median(decimalDate, na.rm = TRUE),
    IQR = IQR(decimalDate, na.rm = TRUE)
  )
#Boxplots for Decimaldate
ggboxplot(DC, x = "partition", y = "decimalDate", 
          color = "partition", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.4,
          order = c("1", "2", "3"),
          ylab = "Year", xlab = "Partition")


## Compute the analysis of variance for DecimalDate
res.aov1 <- aov(decimalDate ~ factor(partition), data = DC)
# Summary of the analysis
summary(res.aov1)

TukeyHSD(res.aov1)

#plot Tukey plot for DecimalDate
plot.tukey <-TukeyHSD(res.aov1)
plot(plot.tukey)


#Annotate the districts for tukey analysis

## Central=1,Southeast=2, Ngamiland=3, Northeast=4, Kgatleng=5, Kweneng=6,  Southern=7 
library (data.table)
DC1 <- data.table(demogFinal3 )
DC1[, districts := ifelse(community %in% c("ranaka", "digawana","molapowabojang","mmankgodi", "mmathethe"), "7",
                          ifelse(community %in% c("lerala", "sefophe", "maunatlala", "ramokgonami","mmadinare", "shoshong","nkange","sebina","mathangwane", "mmandunyane","gumare","rakops","sefhare","tsetsebjwe","nata"), "1",
                                 ifelse(community %in% c("masunga", "tati_siding"), "4",
                                        ifelse(community %in% c("otse"), "2",
                                               ifelse(community %in% c("bokaa","oodi" ), "5",
                                                      ifelse(community %in% c("letlhakeng","lentsweletau","metsimotlhabe"), "6",
                                                             ifelse(community %in% c("gweta", "shakawe"), "3", NA)))))))]
## Compute the analysis of variance for districts
res.aov11 <- aov(partition ~ factor(districts), data = DC1)
# Summary of the analysis
summary(res.aov11)

TukeyHSD(res.aov11)

#plot again using partition as factor
plot.tukey11 <-TukeyHSD(res.aov11)
plot(plot.tukey11)

## Compute the analysis of variance for districts
res.aov110 <- aov(districts ~ factor(partition), data = DC1)
# Summary of the analysis
summary(res.aov110)

TukeyHSD(res.aov110)

#plot Tukey plot for districts
plot.tukey110 <-TukeyHSD(res.aov110)
plot(plot.tukey110)

#combining the pdf graphs
destination = '/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Rplot.pdf'

#open PDF
pdf(file=destination,height = 9, width = 17, onefile = TRUE,)
par( mfrow= c(1,2) )
plot(plot.tukey110)
plot(plot.tukey11)
dev.off()
