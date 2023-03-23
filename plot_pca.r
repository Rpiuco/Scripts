library(plyr)
library(edgeR)
library(ggplot2)

data <- read.table("BRCA/expression_ajcc/matrix.pirna.norxtum.10p",row.name=1,header=T,sep="\t")
pirnaNames <- read.table("BRCA/expression_ajcc/heatmap_005_10.norxtum.all",row.name=1,header=FALSE,sep="\t")
data <- data[, colSums(data != 0) > 0]

dataSelected = t(subset(t(data), select=c(rownames(pirnaNames))))

colorSamples <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")

#myCounting[1,2] == quantidad de levels de controle
#myCounting[2,2] == quantidad de levels de tumor
myCounting <- count(substr(gsub("_[0-9]*","",colnames(dataSelected[,])), 0, 5))

if (nrow(myCounting) == 1){
	groupData <- c(rep(1, times = myCounting[1,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]))
} else if (nrow(myCounting) == 2){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]))
} else if (nrow(myCounting) == 3){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]))
} else if (nrow(myCounting) == 4){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]))
} else if (nrow(myCounting) == 5){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]))
} else if (nrow(myCounting) == 6){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]), rep(6, times = myCounting[6,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]), rep(colorSamples[6], times = myCounting[6,2]))
} else if (nrow(myCounting) == 7){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]), rep(6, times = myCounting[6,2]), rep(7, times = myCounting[7,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]), rep(colorSamples[6], times = myCounting[6,2]), rep(colorSamples[7], times = myCounting[7,2]))
}

group <- factor(groupData)

all_stacked_dge<-DGEList(counts=dataSelected, group=group)
d = calcNormFactors(all_stacked_dge)
dataHeatMap<-cpm(d)
dataHeatMapLog <- log2(dataHeatMap)

dataHeatMapLog[mapply(is.infinite, dataHeatMapLog)] <- 0

df_pca <- prcomp(t(dataHeatMapLog))

df_out <- as.data.frame(df_pca$x)
df_out$group <- group

write.table(df_out,file="output.PCA.BRCA.norxtum.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)

percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )

pdf("plotPCA.BRCA.norxtum.pdf")
p<-ggplot(df_out,aes(x=PC1,y=PC2,color=group ))
p<-p+geom_point() + xlab(percentage[1]) + ylab(percentage[2])
p
dev.off()


##Session_info()