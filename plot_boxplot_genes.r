library(plyr)
library(ggplot2)
library(edgeR)

args <- commandArgs()

data <- read.table(args[[6]],row.name=1,header=T,sep="\t")
listGenes <- read.table(args[[7]],row.name=1,sep="\t")
data <- data[, colSums(data != 0) > 0]

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

y<-DGEList(counts=data, group=group)
all_norm<-calcNormFactors(y, method="TMM")
dataHeatMapAll<-cpm(all_norm)

rownames(listGenes) <- substr(rownames(listGenes[,]), 0, 15)
rownames(dataHeatMapAll) <- substr(rownames(dataHeatMapAll), 0, 15)
listGenes <- listGenes[(rownames(listGenes) %in% rownames(dataHeatMapAll)),]
dataSelected = t(subset(t(dataHeatMapAll), select=c(rownames(listGenes))))
rownames(dataSelected) <- listGenes[,1]

colnames(dataSelected) <- substr(gsub("_[0-9]*","",colnames(dataSelected[,])), 0, 5)
dataSelected <- t(dataSelected)

pdf(paste(args[[8]],"/plots/plotBoxplot.",args[[9]],".pdf", sep = ""))
for (eachGene in colnames(dataSelected)) {

	p <- ggplot(as.data.frame(dataSelected[,eachGene]), aes(x=rownames(dataSelected), y=dataSelected[,eachGene], fill=rownames(dataSelected))) + 
  	geom_boxplot(outlier.shape = NA) +
  	labs(title = "Gene expression by tissue classification", subtitle = eachGene, y = "Normalized Expression (TMM)", x = " ") + 
  	scale_y_continuous(trans='log2') + 
  	guides(fill = guide_legend(title = "Classification", title.position = "top")) + 
  	geom_jitter(shape=16, position=position_jitter(0.2))

  	print(p)
		
}
dev.off()