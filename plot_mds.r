library(plyr)
library(edgeR)

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
	pch <- c(0)
} else if (nrow(myCounting) == 2){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]))
	pch <- c(0,1)
} else if (nrow(myCounting) == 3){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]))
	pch <- c(0,1,2)
} else if (nrow(myCounting) == 4){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]))
	pch <- c(0,1,2,3)
} else if (nrow(myCounting) == 5){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]))
	pch <- c(0,1,2,3,15)
} else if (nrow(myCounting) == 6){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]), rep(6, times = myCounting[6,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]), rep(colorSamples[6], times = myCounting[6,2]))
	pch <- c(0,1,2,3,15,16)
} else if (nrow(myCounting) == 7){
	groupData <- c(rep(1, times = myCounting[1,2]), rep(2, times = myCounting[2,2]), rep(3, times = myCounting[3,2]), rep(4, times = myCounting[4,2]), rep(5, times = myCounting[5,2]), rep(6, times = myCounting[6,2]), rep(7, times = myCounting[7,2]))
	groupDataColor <- c(rep(colorSamples[1], times = myCounting[1,2]), rep(colorSamples[2], times = myCounting[2,2]), rep(colorSamples[3], times = myCounting[3,2]), rep(colorSamples[4], times = myCounting[4,2]), rep(colorSamples[5], times = myCounting[5,2]), rep(colorSamples[6], times = myCounting[6,2]), rep(colorSamples[7], times = myCounting[7,2]))
	pch <- c(0,1,2,3,15,16,17)
}

group <- factor(groupData)

all_stacked_dge = DGEList(counts=dataSelected, group=group)
d = calcNormFactors(all_stacked_dge)
dataHeatMap<-cpm(d)
dataHeatMapLog <- log2(dataHeatMap)

dataHeatMapLog[mapply(is.infinite, dataHeatMapLog)] <- 0

pdf(args[[8]])
plotMDS(dataHeatMapLog, col=colorSamples[group], pch=pch[group])
legend("topleft", legend=c(levels(myCounting$x)), pch=pch, col=colorSamples, ncol=nrow(myCounting))

dev.off()

##Session_info()