library(plyr)
library(edgeR)
library(ggplot2)
library(gplots)

data <- read.table("BRCA/expression_ajcc/matrix.pirna.norxtum.10p",row.name=1,header=T,sep="\t")
pirnaNames <- read.table("BRCA/expression_ajcc/heatmap_005_10.norxtum.all",row.name=1,header=FALSE,sep="\t")
data <- data[, colSums(data != 0) > 0]

dataSelected = t(subset(t(data), select=c(rownames(pirnaNames))))

colorSamples <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")

#myCounting[1,2] == quantidad de levels de controle
#myCounting[2,2] == quantidad de levels de tumor
myCounting <- count(substr(gsub("_[0-9]*","",colnames(data[,])), 0, 4))

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

cormat2 = cor(dataHeatMap, method = "pearson")
cormat3 = cor(dataHeatMapLog, method = "pearson")



scaleRYG <- colorRampPalette(c("green", "darkgoldenrod1", "brown1"), space = "rgb")(256)

pdf("correlation.BRCA.norxtum.pdf")
heatmap.2(x = cormat2, col=scaleRYG, symm = TRUE, main = "Correlation: ALL", cexCol=0.5, cexRow=0.5, margins = c(6, 9), trace="none", density.info="none", dendrogram='none', Rowv=TRUE, Colv=TRUE, labRow=FALSE, labCol=FALSE, ColSideColors = groupDataColor)
legend("topright", legend=c(levels(myCounting$x)), col=colorSamples, fill=colorSamples, ncol=nrow(myCounting))
heatmap.2(x = cormat3, col=scaleRYG, symm = TRUE, main = "Correlation: ALL", cexCol=0.5, cexRow=0.5, margins = c(6, 9), trace="none", density.info="none", dendrogram='none', Rowv=TRUE, Colv=TRUE, labRow=FALSE, labCol=FALSE, ColSideColors = groupDataColor)
legend("topright", legend=c(levels(myCounting$x)), col=colorSamples, fill=colorSamples, ncol=nrow(myCounting))
dev.off()


##Session_info()