library(plyr)
library(edgeR)
library(gplots)

args <- commandArgs()

data <- read.table(args[[6]],row.name=1,header=T,sep="\t")
pirnaNames <- read.table(args[[7]],row.name=1,header=FALSE,sep="\t")

dataSelected = t(subset(t(data), select=c(rownames(pirnaNames))))
dataSelected <- dataSelected[, colSums(dataSelected != 0) > 0]

colorSamples <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628")

#myCounting[1,2] == quantidad de levels de controle
#myCounting[2,2] == quantidad de levels de tumor
myCounting <- count(substr(gsub("_[0-9]*","",colnames(dataSelected[,])), 0, 5))

if (nrow(myCounting) == 1){
	groupData <- c(rep(1, times = myCounting[1,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]))

	newGroup = factor(c(1))
	colorSampleSeparator<-c(rep("red", times = myCounting[1,2]))

} else if (nrow(myCounting) == 2){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]))

	newGroup = factor(c(1,2))
	colorSampleSeparator<-c(rep("red", times = myCounting[1,2]),rep("blue", times = myCounting[2,2]))

} else if (nrow(myCounting) == 3){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]))

	newGroup = factor(c(1,2,3))
	colorSampleSeparator<-c(rep("red", times = myCounting[1,2]),rep("blue", times = myCounting[2,2]),rep("black", times = myCounting[3,2]))

} else if (nrow(myCounting) == 4){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]))

	newGroup = factor(c(1,2,3,4))
	colorSampleSeparator<-c(rep("red", times = myCounting[1,2]),rep("blue", times = myCounting[2,2]),rep("black", times = myCounting[3,2]),rep("brown", times = myCounting[4,2]))

} else if (nrow(myCounting) == 5){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]))

	newGroup = factor(c(1,2,3,4,5))
	colorSampleSeparator<-c(rep("red", times = myCounting[1,2]),rep("blue", times = myCounting[2,2]),rep("black", times = myCounting[3,2]),rep("brown", times = myCounting[4,2]),rep("green", times = myCounting[5,2]))

} else if (nrow(myCounting) == 6){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]), rep(6, times = myCounting[6,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]), rep(colorSamples[6], times = myCounting[6,2]))

	newGroup = factor(c(1,2,3,4,5,6))
	colorSampleSeparator<-c(rep("red", times = myCounting[1,2]),rep("blue", times = myCounting[2,2]),rep("black", times = myCounting[3,2]),rep("brown", times = myCounting[4,2]),rep("green", times = myCounting[5,2]),rep("yellow", times = myCounting[6,2]),rep("white", times = myCounting[7,2]))

} else if (nrow(myCounting) == 7){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]), rep(6, times = myCounting[6,2]), rep(7, times = myCounting[7,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]), rep(colorSamples[6], times = myCounting[6,2]), rep(colorSamples[7], times = myCounting[7,2]))

	newGroup = factor(c(1,2,3,4,5,6,7))
	colorSampleSeparator<-c(rep("red", times = myCounting[1,2]),rep("blue", times = myCounting[2,2]),rep("black", times = myCounting[3,2]),rep("brown", times = myCounting[4,2]),rep("green", times = myCounting[5,2]),rep("yellow", times = myCounting[6,2]),rep("white", times = myCounting[7,2]))

}

group <- factor(groupData)

y<-DGEList(counts=dataSelected, group=group)


if (nrow(myCounting) == 1){
	SampleType1<-rowSums(y$counts[,1:myCounting[1,2]])
	all_stacked<-cbind(SampleType1)
	colnames(all_stacked)<-c(args[[10]])

} else if (nrow(myCounting) == 2){
	SampleType1<-rowSums(y$counts[,1:myCounting[1,2]])
	SampleType2<-rowSums(y$counts[,(myCounting[1,2]+1):(myCounting[1,2]+myCounting[2,2])])
	all_stacked<-cbind(SampleType1, SampleType2)
	colnames(all_stacked)<-c(args[[10]],args[[11]])

} else if (nrow(myCounting) == 3){
	SampleType1<-rowSums(y$counts[,1:myCounting[1,2]])
	SampleType2<-rowSums(y$counts[,(myCounting[1,2]+1):(myCounting[1,2]+myCounting[2,2])])
	SampleType3<-rowSums(y$counts[,(myCounting[2,2]+1):(myCounting[2,2]+myCounting[3,2])])
	all_stacked<-cbind(SampleType1, SampleType2, SampleType3)
	colnames(all_stacked)<-c(args[[10]],args[[11]],args[[12]])

} else if (nrow(myCounting) == 4){
	SampleType1<-rowSums(y$counts[,1:myCounting[1,2]])
	SampleType2<-rowSums(y$counts[,(myCounting[1,2]+1):(myCounting[1,2]+myCounting[2,2])])
	SampleType3<-rowSums(y$counts[,(myCounting[2,2]+1):(myCounting[2,2]+myCounting[3,2])])
	SampleType4<-rowSums(y$counts[,(myCounting[3,2]+1):(myCounting[3,2]+myCounting[4,2])])
	all_stacked<-cbind(SampleType1, SampleType2, SampleType3, SampleType4)
	colnames(all_stacked)<-c(args[[10]],args[[11]],args[[12]],args[[13]])

} else if (nrow(myCounting) == 5){
	SampleType1<-rowSums(y$counts[,1:myCounting[1,2]])
	SampleType2<-rowSums(y$counts[,(myCounting[1,2]+1):(myCounting[1,2]+myCounting[2,2])])
	SampleType3<-rowSums(y$counts[,(myCounting[2,2]+1):(myCounting[2,2]+myCounting[3,2])])
	SampleType4<-rowSums(y$counts[,(myCounting[3,2]+1):(myCounting[3,2]+myCounting[4,2])])
	SampleType5<-rowSums(y$counts[,(myCounting[4,2]+1):(myCounting[4,2]+myCounting[5,2])])
	all_stacked<-cbind(SampleType1, SampleType2, SampleType3, SampleType4, SampleType5)
	colnames(all_stacked)<-c(args[[10]],args[[11]],args[[12]],args[[13]],args[[14]])

} else if (nrow(myCounting) == 6){
	SampleType1<-rowSums(y$counts[,1:myCounting[1,2]])
	SampleType2<-rowSums(y$counts[,(myCounting[1,2]+1):(myCounting[1,2]+myCounting[2,2])])
	SampleType3<-rowSums(y$counts[,(myCounting[2,2]+1):(myCounting[2,2]+myCounting[3,2])])
	SampleType4<-rowSums(y$counts[,(myCounting[3,2]+1):(myCounting[3,2]+myCounting[4,2])])
	SampleType5<-rowSums(y$counts[,(myCounting[4,2]+1):(myCounting[4,2]+myCounting[5,2])])
	SampleType6<-rowSums(y$counts[,(myCounting[5,2]+1):(myCounting[5,2]+myCounting[6,2])])
	all_stacked<-cbind(SampleType1, SampleType2, SampleType3, SampleType4, SampleType5, SampleType6)
	colnames(all_stacked)<-c(args[[10]],args[[11]],args[[12]],args[[13]],args[[14]],args[[15]])

} else if (nrow(myCounting) == 7){
	SampleType1<-rowSums(y$counts[,1:myCounting[1,2]])
	SampleType2<-rowSums(y$counts[,(myCounting[1,2]+1):(myCounting[1,2]+myCounting[2,2])])
	SampleType3<-rowSums(y$counts[,(myCounting[2,2]+1):(myCounting[2,2]+myCounting[3,2])])
	SampleType4<-rowSums(y$counts[,(myCounting[3,2]+1):(myCounting[3,2]+myCounting[4,2])])
	SampleType5<-rowSums(y$counts[,(myCounting[4,2]+1):(myCounting[4,2]+myCounting[5,2])])
	SampleType6<-rowSums(y$counts[,(myCounting[5,2]+1):(myCounting[5,2]+myCounting[6,2])])
	SampleType7<-rowSums(y$counts[,(myCounting[6,2]+1):(myCounting[6,2]+myCounting[7,2])])
	all_stacked<-cbind(SampleType1, SampleType2, SampleType3, SampleType4, SampleType5, SampleType6, SampleType7)
	colnames(all_stacked)<-c(args[[10]],args[[11]],args[[12]],args[[13]],args[[14]],args[[15]],args[[16]])

}


all_stacked_dge<-DGEList(counts=all_stacked, group=newGroup)

all_stacked_norm<-calcNormFactors(all_stacked_dge, method="TMM")
dataHeatMapStacked<-cpm(all_stacked_norm)

all_norm<-calcNormFactors(y, method="TMM")
dataHeatMapAll<-cpm(all_norm)
dataHeatMapAllLog <- log2(dataHeatMapAll)
dataHeatMapStackedLog <- log2(dataHeatMapStacked)

dataHeatMapAllLog[mapply(is.infinite, dataHeatMapAllLog)] <- 0
dataHeatMapStackedLog[mapply(is.infinite, dataHeatMapStackedLog)] <- 0

scaleRYG <- colorRampPalette(c("green", "darkgoldenrod1","brown1"), space = "rgb")(256)

pdf(args[[8]])
heatmap.2(
	dataHeatMapStackedLog,
	main=paste(args[[9]], ": Stacked", sep=""),
	trace="none", 
	density.info="none", 
	cexCol=0.5, 
	cexRow=0.5, 
	col=scaleRYG, 
	margins = c(6, 9),
	Rowv=TRUE, 
	Colv=TRUE
)
heatmap.2(
	dataHeatMapAllLog, 
	main=args[[9]],
	trace="none", 
	density.info="none", 
	cexCol=0.5, 
	cexRow=0.5, 
	col=scaleRYG, 
	margins = c(6, 9),
	labCol = FALSE, 
	Rowv=TRUE, 
	Colv=TRUE, 
	ColSideColors = colorSampleSeparator
)

dev.off()