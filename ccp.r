library(plyr)
library(edgeR)
library(ConsensusClusterPlus)

args <- commandArgs()

data <- read.table(args[[6]],row.name=1,header=T,sep="\t")
data <- data[, colSums(data != 0) > 0]

#myCounting[1,2] == quantidad de levels de controle
#myCounting[2,2] == quantidad de levels de tumor
myCounting <- count(substr(gsub("_[0-9]*","",colnames(data[,])), 0, 5))

if (nrow(myCounting) == 1){
        groupData <- c(rep(1, times = myCounting[1,2]))
} else if (nrow(myCounting) == 2){
        groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]))
} else if (nrow(myCounting) == 3){
        groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]))
} else if (nrow(myCounting) == 4){
        groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]))
} else if (nrow(myCounting) == 5){
        groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]))
} else if (nrow(myCounting) == 6){
        groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]), rep(6, times = myCounting[6,2]))
} else if (nrow(myCounting) == 7){
        groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]), rep(6, times = myCounting[6,2]), rep(7, times = myCounting[7,2]))
}

group <- factor(groupData)

d = DGEList(counts=data, group=group)
d = calcNormFactors(d)
d = cpm(d)

dataProteinCodingFilter = sweep(d,1, apply(d,1,median,na.rm=T))
dataProteinCodingMatrix <- as.matrix(dataProteinCodingFilter)

pdf(paste(args[[7]], args[[8]], sep=""))
results = ConsensusClusterPlus(dataProteinCodingMatrix, maxK=10, reps=1000, pItem=0.80, pFeature=1, title=paste("Cluster Consensus: ", args[[9]], sep=""), clusterAlg="hc", distance="pearson", seed=1262118388.71279)

iclProteinCoding = calcICL(results,title="ICL : ProteinCoding")
dev.off()


write.table(results[[2]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kTwo.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(results[[3]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kThree.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(results[[4]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kFour.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(results[[5]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kFive.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(results[[6]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kSix.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(results[[7]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kSeven.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(results[[8]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kEigth.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(results[[9]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kNine.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
write.table(results[[10]]$consensusClass,file=paste(args[[7]],"outputCCP.",args[[10]],".kTen.CClass.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)

write.table(iclProteinCoding$clusterConsensus,file=paste(args[[7]],"outputCCP.",args[[10]],".ClusterC.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)

write.table(iclProteinCoding$itemConsensus,file=paste(args[[7]],"outputCCP.",args[[10]],".ItemC.txt", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)

##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5680778/ -> Usou esses parametros

#Transcription levels quantified by RSEM were filtered to remove genes whose expression was quantified as zero by RSEM in more than 75% of the tumor samples, reducing the set of genes from 20,531 to 15,951. Gene quantifications were subsequently log2 transformed, with zero values set to missing. To identify genes whose expression was variable, the gene set was filtered to remove genes that demonstrated a standard deviation below 2.0 across all tumor samples, resulting in a set of 1,868 genes with high variability in expression.
#The log2 transformed expression values were then median centered prior to clustering analysis.
#Cluster analysis was performed using ConsensusClusterPlus (Wilkerson and Hayes, 2010), using agglomerative hierarchical clustering with a 1-Pearson correlation distances and resampling 80% of the samples for 1000 repetitions. The optimal number of clusters was determined using the empirical cumulative distribution function plot.

##Session_info()