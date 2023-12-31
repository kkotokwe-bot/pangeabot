
#FINAL mlesky Script- 3 partitions 
# 23/11/2022 mlesky with bootstraps and CI

setwd("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Analysis for Gag Final/Mlesky")




library(ape)
library(treedater)
library(treestructure)
library(msa)
library(Biostrings)
library(mlesky)
library(vctrs)
library(ggplot2)
library(lubridate)

geneRegion <- "gag"
nboot<- 100
set.seed(5492)

ts_dft <- read.csv(paste(geneRegion,"_treestructure_splits1.csv", sep=""), stringsAsFactors = F)


for (i in 1:max(as.numeric(ts_dft$partition))){
  
  # run IQtree with bootstrap
  prefix <- paste0("seq_partition_",i) 
  mySeq <- read.dna(paste0(prefix, ".fas"), format="fasta", as.matrix= F)
  
  # system(paste0("/Users/peggy/Downloads/iqtree-2.1.3-MacOSX/bin/iqtree2 -m gtr -s ", 
  #              paste0(fasta_noCGR, ".fas")," --prefix ", fasta_noCGR ," -redo --mem 80% -T 4 --fast")) 
  
  system(paste0("/Users/peggy/Downloads/iqtree-2.1.3-MacOSX/bin/iqtree2 -m gtr -s ", 
                paste0(prefix,".fas"), " --prefix ", prefix ," -redo --mem 80% -T 4 --fast -b ",nboot)) # note we are doing bootstraps here
  
  mltree <- read.tree(paste0(prefix, ".treefile"))
  bootTrees= read.tree(paste0(prefix, ".boottrees"))
  
  # run treedater
  
  sts <- as.numeric(unlist(lapply( mltree$tip.label, function(x) {strsplit(x, "_")[[1]][2]})))
  names(sts) <-  mltree$tip.label
  
  mltree_dated <- dater( mltree , sts, s = length(mySeq[[1]]), 
                         numStartConditions = 1, maxit = 1, meanRateLimits = c(0.001,  0.003),
                         ncpu = 1)
  
  
  boot_dated <- bootTrees
  tMRCAs <- c()
  for (j in 1:nboot){
    print(j)
    sts <- as.numeric(unlist(lapply( bootTrees[[j]]$tip.label, function(x) {strsplit(x, "_")[[1]][2]})))
    names(sts) <-  bootTrees[[j]]$tip.label
    
    boot_dated[[j]] <- dater( bootTrees[[j]] , sts, s = length(mySeq[[1]]), 
                              numStartConditions = 1, maxit = 1, meanRateLimits = c(0.001,  0.003),
                              ncpu = 1)
    
    tMRCAs <- c(tMRCAs, boot_dated[[j]]$timeOfMRCA) ## this it the time of MRCAs
  }
  
  write.csv(data.frame("tMRCA"=tMRCAs), paste(geneRegion, prefix, "csv", sep = "." ))
  write.tree(boot_dated, paste(geneRegion, "_partition", i, "_datedTrees.nwk", sep=""))
  
  
  ###dropCGR sequences here from phylogenies  (in ML+ all 100 boot) 
  mltree_dated$intree$edge.length <- mltree_dated$edge.length
  mltree_dated2 <- mltree_dated$intree
  mltree_dated_noCGR <- drop.tip(mltree_dated2, mltree_dated2$tip.label[grep("CGR", mltree_dated2$tip.label)])
  
  boot_dated_noCGR <-  boot_dated
  for (j in 1:nboot){
    boot_dated_noCGR[[j]]$intree$edge.length <- boot_dated_noCGR[[j]]$edge.length
    boot_dated_noCGR[[j]]<- drop.tip(boot_dated_noCGR[[j]]$intree, boot_dated_noCGR[[j]]$intree$tip.label[grep("CGR", boot_dated_noCGR[[j]]$intree$tip.label)])
  }
  
  # run mlesky on mltree + each boot tree
  ml_sts <- as.numeric(unlist(lapply(mltree_dated_noCGR$tip.label, function(x) {strsplit(x, "_")[[1]][2]})))
  names(ml_sts) <- mltree_dated_noCGR$tip.label
  ml_mlesky <- mlskygrid(mltree_dated_noCGR, tau = NULL, tau_lower=.001, tau_upper = 10, 
                         res = 100, ncpu = 1, sampleTimes = ml_sts)
  taxis <- ml_mlesky$time 
  
  
  res = pbmcapply::pbmclapply( 1:length(boot_dated_noCGR), function(irep){
    
    sts <- as.numeric(unlist(lapply(boot_dated_noCGR[[irep]]$tip.label, function(x) {strsplit(x, "_")[[1]][2]})))
    names(sts) <- boot_dated_noCGR[[irep]]$tip.label
    
    
    f1 <- mlskygrid( boot_dated_noCGR[[irep]], sampleTimes = sts, res = ml_mlesky$res,
                     tau = ml_mlesky$tau, tau_tol = ml_mlesky$tau_tol , ncross = ml_mlesky$ncross, 
                     ne0 = median( ml_mlesky$ne ), adapt_time_axis = FALSE, ncpu = 1)
    af <- approxfun( f1$time, f1$ne, rule = 2)
    afgr <- approxfun( f1$time, f1$growthrate , rule =2)
    list(ne = af(taxis),  growthrate = afgr(taxis) )
  }, mc.cores = 1 )
  
  nemat <- do.call( cbind, lapply( res, '[[', 'ne' ) )
  grmat <- do.call( cbind, lapply( res, '[[', 'growthrate' ))
  
  ml_mlesky$ne_ci <- cbind(
    nelb = unname(  apply( nemat, 1, FUN = function(x) quantile(x, c(.025))) )
    , ne = unname(  apply( nemat, 1, FUN = function(x) quantile(x, c(.50))) )
    , neub = unname(  apply( nemat, 1, FUN = function(x) quantile(x, c(.975))) )
  )
  
  ml_mlesky$growthrate_ci <- cbind( 
    grlb = unname(  apply( grmat, 1, FUN = function(x) quantile(x, c(.025))) )
    , gr = unname(  apply( grmat, 1, FUN = function(x) quantile(x, c(.50))) )
    , grub = unname(  apply( grmat, 1, FUN = function(x) quantile(x, c(.975))) )
  )
  
  saveRDS( ml_mlesky, file = paste0(prefix, ".rds" ) )
  
}

pdf(paste0(geneRegion,"_plotsfinal_","_mlesky.pdf",sep=""),height = 9, width = 8, onefile = TRUE,)

f = readRDS(paste0("seq_partition_",1, ".rds" ))
pne = plot(f, logy=F, ggplot=T) + theme_classic() + xlab('')
pne$data$date <- as.Date( date_decimal( pne$data$t ))
(pne1 <- ggplot( data = pne$data , aes(x = date)),  pne2 <- ggplot( data = pne$data , aes(x = date)), pne3 <- ggplot( data = pne$data , aes(x = date)) + 
    
    geom_ribbon(aes(ymin = nelb, ymax = neub), alpha = .50, fill = '#DD6200' , colour='white') + 
    geom_line( aes(y =  nemed), lwd = .8) + 
    theme_classic() + xlab('') + ylab('Effective population size' )
)

f = readRDS(paste0("seq_partition_",2, ".rds" ))
pne = plot(f, logy=F, ggplot=T) + theme_classic() + xlab('')
pne$data$date <- as.Date( date_decimal( pne$data$t ))
(pne2 <- ggplot( data = pne$data , aes(x = date) ) + 
    geom_ribbon(aes(ymin = nelb, ymax = neub), alpha = .50, fill = '#2F4F4F' , colour='white') + 
    geom_line( aes(y =  nemed), lwd = .8) + 
    theme_classic() + xlab('') + ylab('Effective population size' )
)

f = readRDS(paste0("seq_partition_",3, ".rds" ))
pne = plot(f, logy=F, ggplot=T) + theme_classic() + xlab('')
pne$data$date <- as.Date( date_decimal( pne$data$t ))
(pne3 <- ggplot( data = pne$data , aes(x = date) ) + 
    geom_ribbon(aes(ymin = nelb, ymax = neub), alpha = .50, fill = '#DD6200' , colour='white') + 
    geom_line( aes(y =  nemed), lwd = .8) + 
    theme_classic() + xlab('') + ylab('Effective population size' )
)
dev.off() 
ggsave( pne1,pne2,pne3, file = glue(prefix, 'b117.1-pne.png') , width = 3.5, height = 2.8)
ggsave( pne1, file = glue(prefix, 'b117.1-pne.pdf') , width = 3.5, height = 2.8)


#plotting in one pdf file
pdf(paste0(geneRegion,"_partitionplotsgagfinal_22","_mlesky.pdf",sep=""),height = 9, width = 10, onefile = TRUE,)
destination = '/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Analysis for Gag Final/Mlesky/pol_partitionplotspolfinal_221_mlesky.pdf'

#open PDF
pdf(file=destination)

#specify to save plots in 2x2 grid
par(mfrow = c(2,3))

# frst plot with one set from the aplit (1-4)
plt <-readRDS(paste(geneRegion,"_partition", 1, "_mlesky.rds",  sep=""))
plot(plt$time, log(plt$ne), type='l',cex.axis= 2.1, cex.lab = 1.3, lwd=5.0, cex.sub=1, xlim=c(1970,2020), ylab="Effective pop size (gag)", xlab="Time", ylim=c(0,30), col="#2F4F4F", bty="n")
plt <-readRDS(paste(geneRegion,"_partition", 2, "_mlesky.rds",  sep=""))
lines(plt$time, log(plt$ne), type="l",lwd=5.0, col="#FF7F50")
plt <-readRDS(paste(geneRegion,"_partition", 3, "_mlesky.rds",  sep=""))
lines(plt$time, log(plt$ne), type="l",lwd=5.0, col="#8856a7")
legend("topleft", c("Partition 1","Partition 2", "Partition 3"), fill=c("#2F4F4F","#FF7F50","#8856a7"),cex= 2,bty="n")
labels( c("1","2","3"))

setwd("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Treestructure Final Analysis_01122021/Mlesky Final ")

# frst plot with one set from the aplit (1-4)
plt <-readRDS(paste0("seq_partition_",1, ".rds" ))
plot(plt$time, log(plt$ne), type='l',cex.axis= 2.1, cex.lab =1.3, lwd=5.0,cex.sub=1,  xlim=c(1970,2020), ylab="Effective pop size (POL)", xlab="Time", ylim=c(0,40), col="#2F4F4F", bty="n")
plt <-readRDS(paste0("seq_partition_",2, ".rds" ))
lines(plt$time, log(plt$ne), type="l",lwd=5.0, col="#FF7F50")
plt <-readRDS(paste0("seq_partition_",3, ".rds" ))
lines(plt$time, log(plt$ne), type="l",lwd=5.0, col="#8856a7")
legend("topleft", c ("Partition 1","Partition 2", "Partition 3"), fill=c("#2F4F4F","#FF7F50","#8856a7"),cex= 2, bty="n")
labels( c("1","2","3"))
dev.off()
####New Plot ANalysis#####

f1 <- readRDS("seq_partition_1.rds")
f2 <- readRDS("seq_partition_2.rds")
f3 <- readRDS("seq_partition_3.rds")
for (i in 1:3){
  ftemp <- get(paste0("f",i))
  pne = plot(ftemp, logy=F, ggplot=T) + theme_classic() + xlab('')
  pne$data$date <- as.Date( date_decimal( pne$data$t ))
  assign(paste0("pne",i), pne)
}
plot(pne1$data$date, pne1$data$nemed, col="blue", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")
#lines(pne1$data$date, pne1$data$nelb, col="lightblue")
#lines(pne1$data$date, pne1$data$neub, col="lightblue")

polygon(x = c(pne1$data$date, rev(pne1$data$date)),
        y = c(pne1$data$nelb,   rev(pne1$data$neub)),
        col =  yarrr::transparent("blue", trans.val = .9), border=NA)

lines(pne2$data$date, pne2$data$nemed, col="red", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")


polygon(x = c(pne2$data$date, rev(pne2$data$date)),
        y = c(pne2$data$nelb,   rev(pne2$data$neub)),
        col =  yarrr::transparent("red", trans.val = .9), border=NA)

lines(pne3$data$date, pne3$data$nemed, col="darkgreen", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")
polygon(x = c(pne3$data$date, rev(pne3$data$date)),
        y = c(pne3$data$nelb,   rev(pne3$data$neub)),
        col =  yarrr::transparent("darkgreen", trans.val = .9), border=NA)

legend("left", legend = c("Partition 1", "Partition 2", "Partition 3"),
       fill=c("blue", "red", "darkgreen"), bty="n")

#####Combine the graph in one#######
#plotting in one pdf file
pdf(paste("Rplot021 copy.pdf"),height = 2, width = 2, onefile = TRUE,)
destination = ("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Analysis for Gag Final/Mlesky/Rplot021 copy.pdf")

#open PDF
pdf(file=destination)

#specify to save plots in 2x2 grid
par(mfrow = c(3,2))

setwd("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Analysis for Gag Final/Mlesky_1")

# First plot,gag1 
f1 <- readRDS("seq_partition_1.rds")
f2 <- readRDS("seq_partition_2.rds")
f3 <- readRDS("seq_partition_3.rds")
for (i in 1:3){
  ftemp <- get(paste0("f",i))
  pne = plot(ftemp, logy=F, ggplot=T) + theme_classic() + xlab('')
  pne$data$date <- as.Date( date_decimal( pne$data$t ))
  assign(paste0("pne",i), pne)
}
plot(pne1$data$date, pne1$data$nemed, col="blue", type="l", ylim=c(0,1000000),xlab="Time",main="gag region", ylab="Effective population size", sub= "Precision parameter=1")
#lines(pne1$data$date, pne1$data$nelb, col="lightblue")
#lines(pne1$data$date, pne1$data$neub, col="lightblue")

polygon(x = c(pne1$data$date, rev(pne1$data$date)),
        y = c(pne1$data$nelb,   rev(pne1$data$neub)),
        col =  yarrr::transparent("blue", trans.val = .9), border=NA)

lines(pne2$data$date, pne2$data$nemed, col="red", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")


polygon(x = c(pne2$data$date, rev(pne2$data$date)),
        y = c(pne2$data$nelb,   rev(pne2$data$neub)),
        col =  yarrr::transparent("red", trans.val = .9), border=NA)

lines(pne3$data$date, pne3$data$nemed, col="darkgreen", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")
polygon(x = c(pne3$data$date, rev(pne3$data$date)),
        y = c(pne3$data$nelb,   rev(pne3$data$neub)),
        col =  yarrr::transparent("darkgreen", trans.val = .9), border=NA)

legend("left", legend = c("Partition 1", "Partition 2", "Partition 3"),
       fill=c("blue", "red", "darkgreen"), bty="n")

setwd("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Treestructure Final Analysis_01122021/Mlesky Final_1")
# First plot,pol1
f1 <- readRDS("seq_partition_1.rds")
f2 <- readRDS("seq_partition_2.rds")
f3 <- readRDS("seq_partition_3.rds")
for (i in 1:3){
  ftemp <- get(paste0("f",i))
  pne = plot(ftemp, logy=F, ggplot=T) + theme_classic() + xlab('')
  pne$data$date <- as.Date( date_decimal( pne$data$t ))
  assign(paste0("pne",i), pne)
}
plot(pne1$data$date, pne1$data$nemed, col="blue", type="l", ylim=c(0,1000000),xlab="Time",main="pol region", ylab="Effective population size", sub= "Precision parameter=1")
#lines(pne1$data$date, pne1$data$nelb, col="lightblue")
#lines(pne1$data$date, pne1$data$neub, col="lightblue")

polygon(x = c(pne1$data$date, rev(pne1$data$date)),
        y = c(pne1$data$nelb,   rev(pne1$data$neub)),
        col =  yarrr::transparent("blue", trans.val = .9), border=NA)

lines(pne2$data$date, pne2$data$nemed, col="red", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")


polygon(x = c(pne2$data$date, rev(pne2$data$date)),
        y = c(pne2$data$nelb,   rev(pne2$data$neub)),
        col =  yarrr::transparent("red", trans.val = .9), border=NA)

lines(pne3$data$date, pne3$data$nemed, col="darkgreen", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")
polygon(x = c(pne3$data$date, rev(pne3$data$date)),
        y = c(pne3$data$nelb,   rev(pne3$data$neub)),
        col =  yarrr::transparent("darkgreen", trans.val = .9), border=NA)

legend("left", legend = c("Partition 1", "Partition 2", "Partition 3"),
       fill=c("blue", "red", "darkgreen"), bty="n")



setwd("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Analysis for Gag Final/Mlesky")

# second plot, gag2 
f1 <- readRDS("seq_partition_1.rds")
f2 <- readRDS("seq_partition_2.rds")
f3 <- readRDS("seq_partition_3.rds")
for (i in 1:3){
  ftemp <- get(paste0("f",i))
  pne = plot(ftemp, logy=F, ggplot=T) + theme_classic() + xlab('')
  pne$data$date <- as.Date( date_decimal( pne$data$t ))
  assign(paste0("pne",i), pne)
}
plot(pne1$data$date, pne1$data$nemed, col="blue", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size", sub= "Precision parameter=10")
#lines(pne1$data$date, pne1$data$nelb, col="lightblue")
#lines(pne1$data$date, pne1$data$neub, col="lightblue")

polygon(x = c(pne1$data$date, rev(pne1$data$date)),
        y = c(pne1$data$nelb,   rev(pne1$data$neub)),
        col =  yarrr::transparent("blue", trans.val = .9), border=NA)

lines(pne2$data$date, pne2$data$nemed, col="red", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size",)


polygon(x = c(pne2$data$date, rev(pne2$data$date)),
        y = c(pne2$data$nelb,   rev(pne2$data$neub)),
        col =  yarrr::transparent("red", trans.val = .9), border=NA)

lines(pne3$data$date, pne3$data$nemed, col="darkgreen", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")
polygon(x = c(pne3$data$date, rev(pne3$data$date)),
        y = c(pne3$data$nelb,   rev(pne3$data$neub)),
        col =  yarrr::transparent("darkgreen", trans.val = .9), border=NA)

legend("left", legend = c("Partition 1", "Partition 2", "Partition 3"),
       fill=c("blue", "red", "darkgreen"), bty="n")

setwd ("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Treestructure Final Analysis_01122021/Mlesky Final ") 

# second plot pol2
f1 <- readRDS("seq_partition_1.rds")
f2 <- readRDS("seq_partition_2.rds")
f3 <- readRDS("seq_partition_3.rds")
for (i in 1:3){
  ftemp <- get(paste0("f",i))
  pne = plot(ftemp, logy=F, ggplot=T) + theme_classic() + xlab('')
  pne$data$date <- as.Date( date_decimal( pne$data$t ))
  assign(paste0("pne",i), pne)
}
plot(pne1$data$date, pne1$data$nemed, col="blue", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size", sub= "Precision parameter=10")
#lines(pne1$data$date, pne1$data$nelb, col="lightblue")
#lines(pne1$data$date, pne1$data$neub, col="lightblue")

polygon(x = c(pne1$data$date, rev(pne1$data$date)),
        y = c(pne1$data$nelb,   rev(pne1$data$neub)),
        col =  yarrr::transparent("blue", trans.val = .9), border=NA)

lines(pne2$data$date, pne2$data$nemed, col="red", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size",)


polygon(x = c(pne2$data$date, rev(pne2$data$date)),
        y = c(pne2$data$nelb,   rev(pne2$data$neub)),
        col =  yarrr::transparent("red", trans.val = .9), border=NA)

lines(pne3$data$date, pne3$data$nemed, col="darkgreen", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")
polygon(x = c(pne3$data$date, rev(pne3$data$date)),
        y = c(pne3$data$nelb,   rev(pne3$data$neub)),
        col =  yarrr::transparent("darkgreen", trans.val = .9), border=NA)

legend("left", legend = c("Partition 1", "Partition 2", "Partition 3"),
       fill=c("blue", "red", "darkgreen"), bty="n")




setwd("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Analysis for Gag Final/Mlesky_20")

# third plot, gag3

f1 <- readRDS("seq_partition_1.rds")
f2 <- readRDS("seq_partition_2.rds")
f3 <- readRDS("seq_partition_3.rds")
for (i in 1:3){
  ftemp <- get(paste0("f",i))
  pne = plot(ftemp, logy=F, ggplot=T) + theme_classic() + xlab('')
  pne$data$date <- as.Date( date_decimal( pne$data$t ))
  assign(paste0("pne",i), pne)
}
plot(pne1$data$date, pne1$data$nemed, col="blue", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size", sub= "Precision parameter=20")
#lines(pne1$data$date, pne1$data$nelb, col="lightblue")
#lines(pne1$data$date, pne1$data$neub, col="lightblue")

polygon(x = c(pne1$data$date, rev(pne1$data$date)),
        y = c(pne1$data$nelb,   rev(pne1$data$neub)),
        col =  yarrr::transparent("blue", trans.val = .9), border=NA)

lines(pne2$data$date, pne2$data$nemed, col="red", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")


polygon(x = c(pne2$data$date, rev(pne2$data$date)),
        y = c(pne2$data$nelb,   rev(pne2$data$neub)),
        col =  yarrr::transparent("red", trans.val = .9), border=NA)

lines(pne3$data$date, pne3$data$nemed, col="darkgreen", type="l", ylim=c(0,1000000),xlab="Time",ylab="Effective population size")
polygon(x = c(pne3$data$date, rev(pne3$data$date)),
        y = c(pne3$data$nelb,   rev(pne3$data$neub)),
        col =  yarrr::transparent("darkgreen", trans.val = .9), border=NA)

legend("left", legend = c("Partition 1", "Partition 2", "Partition 3"),
       fill=c("blue", "red", "darkgreen"), bty="n")




setwd("/Users/peggy/Downloads/Treestructure/Treestructure Final Analysis_01122021/Treestructure Final Analysis_01122021/Mlesky Final_20")
# Third plot, pol3  

f1 <- readRDS("seq_partition_1.rds")
f2 <- readRDS("seq_partition_2.rds")
f3 <- readRDS("seq_partition_3.rds")
for (i in 1:3){
  ftemp <- get(paste0("f",i))
  pne = plot(ftemp, logy=F, ggplot=T) + theme_classic() + xlab('')
  pne$data$date <- as.Date( date_decimal( pne$data$t ))
  assign(paste0("pne",i), pne)
}
plot(pne1$data$date, pne1$data$nemed, col="blue", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size", sub= "Precision parameter=20")
#lines(pne1$data$date, pne1$data$nelb, col="lightblue")
#lines(pne1$data$date, pne1$data$neub, col="lightblue")

polygon(x = c(pne1$data$date, rev(pne1$data$date)),
        y = c(pne1$data$nelb,   rev(pne1$data$neub)),
        col =  yarrr::transparent("blue", trans.val = .9), border=NA)

lines(pne2$data$date, pne2$data$nemed, col="red", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")


polygon(x = c(pne2$data$date, rev(pne2$data$date)),
        y = c(pne2$data$nelb,   rev(pne2$data$neub)),
        col =  yarrr::transparent("red", trans.val = .9), border=NA)

lines(pne3$data$date, pne3$data$nemed, col="darkgreen", type="l", ylim=c(0,1000000),xlab="Time", ylab="Effective population size")
polygon(x = c(pne3$data$date, rev(pne3$data$date)),
        y = c(pne3$data$nelb,   rev(pne3$data$neub)),
        col =  yarrr::transparent("darkgreen", trans.val = .9), border=NA)

legend("left", legend = c("Partition 1", "Partition 2", "Partition 3"),
       fill=c("blue", "red", "darkgreen"), bty="n")
dev.off()










###CI Calculations (TMRCA)####
install.packages("Rmisc") 
library(Rmisc)

gfinal <-read.csv(paste(geneRegion,".seq_partition_20.csv", sep=""), stringsAsFactors = T)
final_CI <-CI(gfinal$tMRCA, ci=0.95)

