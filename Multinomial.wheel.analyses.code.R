### Multinomial wheel analyses and plot using the original count matrix "ori_matrix" and the clustering result from Seurat "Seurat_matrix"

library(Seurat)

DimPlot(Seurat_matrix.a, reduction = "umap", pt.size=0.1, label=TRUE, label.size=5) 
nclusters <- nlevels(Seurat_matrix.a@meta.data$seurat_clusters) ## nclusters==11 in this dataset

temp <- CellsByIdentities(Seurat_matrix.a, idents = 0)
names(temp) <- "bc"
clusters <- data.frame("bc"=temp, "g"=0, stringsAsFactors=FALSE)
for (i in 1:(nclusters - 1)) {
	temp <- CellsByIdentities(Seurat_matrix.a, idents = i)
	names(temp) <- "bc"
	clusters <- rbind.data.frame(clusters, data.frame("bc"=temp, "g"=i, stringsAsFactors=FALSE))
}

grouped_data <- aggregate(clusters$bc, by=list(clusters$g), FUN=length)
names(grouped_data) <- c("g", "cell.count")

### make aggregated vector for each Seurat cluster

dna.aggregate <- matrix(0, nrow=300, ncol=nrow(grouped_data))
for (i in 1:nrow(grouped_data)) {
	theseCells <- clusters$bc[which(clusters$g == grouped_data$g[i])]
	dna.aggregate[, i] <- rowSums(ori_matrix[, theseCells])
}

df <- ori_matrix[, which(colnames(ori_matrix) %in% colnames(Seurat_matrix))] 

### make multinomial vector dictionary

multistates <- matrix(0, nrow=9*nclusters*(nclusters-1)/2+nclusters, ncol=4)
multimatrix <- matrix(0, nrow=300, ncol=nrow(multistates)) 
rowpointer <- 1
for (i in 0:(nclusters-1)) {
	multistates[rowpointer,1] <- rowpointer
	multistates[rowpointer,2] <- i
	multistates[rowpointer,3] <- i
	multistates[rowpointer,4] <- 0
	temp <- dna.aggregate[,i+1]
	multimatrix[,rowpointer] <- temp / sum(temp)
	rowpointer <- rowpointer + 1
}
for (i in 0:(nclusters-2)) {
	temp1 <- dna.aggregate[,i+1]
	temp1 <- temp1 / sum(temp1)
	for (j in (i+1):(nclusters-1)) {
		temp2 <- dna.aggregate[,j+1]
		temp2 <- temp2 / sum(temp2)
		for (k in 1:9) {
			temp3 <- ((10 - k)*temp1 + k*temp2) / 10
			multistates[rowpointer,1] <- rowpointer
			multistates[rowpointer,2] <- i
			multistates[rowpointer,3] <- j
			multistates[rowpointer,4] <- k
			multimatrix[,rowpointer] <- temp3
			rowpointer <-  rowpointer + 1
		}				
	}
}
dimnames(multistates)[[2]] <- c("pointer", "start.cluster", "end.cluster", "distance")

### END of making multinomial vector dictionary

DNA.dist.results <- matrix(0, nrow=ncol(df), ncol=9)

for (i in 1:ncol(df)) {
	thisBc <- names(df)[i]
	cat(i, thisBc, "\t")
	thisCountVector <- df[, i]
	temp.best.pval <- -10000
	best.pointer <- -1
	random2 <- sample.int(nclusters, 2)
	thisCluster <- clusters$g[which(clusters$bc == thisBc)]

	for (k in 1:nrow(multistates)) {
		if ((multistates[k,2] == thisCluster) | (multistates[k,3] == thisCluster)) {
			testProb <- multimatrix[,multistates[k,1]]
			testPval <- dmultinom(x=thisCountVector, prob=testProb, log=T)
			if (testPval > temp.best.pval) {
				temp.best.pval <- testPval
				best.pointer <- multistates[k,1]
			}
		}
	}
	
	DNA.dist.results[i, 1] <- i
	DNA.dist.results[i, 2] <- best.pointer
	DNA.dist.results[i, 3] <- temp.best.pval
	DNA.dist.results[i, 4] <- thisCluster

	if (thisCluster == multistates[best.pointer,2]) {
		DNA.dist.results[i, 5] <- thisCluster
		DNA.dist.results[i, 6] <- multistates[best.pointer,3]
		DNA.dist.results[i, 7] <- multistates[best.pointer,4]
		if (DNA.dist.results[i, 7] <= 5) {
			DNA.dist.results[i, 8] <- 1
		} else {
			DNA.dist.results[i, 8] <- -1
		}
	} else if (thisCluster == multistates[best.pointer,3]) {
		DNA.dist.results[i, 5] <- thisCluster
		DNA.dist.results[i, 6] <- multistates[best.pointer,2]
		DNA.dist.results[i, 7] <- 10 - multistates[best.pointer,4]
		if (DNA.dist.results[i, 7] <= 5) {
			DNA.dist.results[i, 8] <- 1
		} else {
			DNA.dist.results[i, 8] <- -1
		}
	} 
	if (DNA.dist.results[i, 5] == DNA.dist.results[i, 6]) {
		if ((random2[1]-1) == DNA.dist.results[i,5]) {
			DNA.dist.results[i,6] <- random2[2]-1
		} else {
			DNA.dist.results[i,6] <- random2[1]-1
		} 
	}
}

dimnames(DNA.dist.results)[[2]] <- c("cell.num", "pointer.in.multistates", "best.pval", "assigned.cluster", "home.cluster", "other.cluster", "mix.tenths", "correct.by.seurat","patient")

DNAwheel <- cbind.data.frame("bc"=colnames(df), data.frame(DNA.dist.results), stringsAsFactors=F)

### plot multinomial wheel

x11 <- cos(2 * pi * seq(0, 10, 1) / 11)
y11 <- sin(2 * pi * seq(0, 10, 1) / 11)
# plot(x11, y11) shows this is correct

xyCoords <- matrix(0, ncol=2, nrow=nrow(DNAwheel))
thisColors <- rep("black", nrow(DNAwheel))
thisRatios <- rep(0, nrow(DNAwheel))

for (i in 1:nrow(DNAwheel)) {
	x1 <- x11[DNAwheel$home.cluster[i] + 1]
	y1 <- y11[DNAwheel$home.cluster[i] + 1]
	x2 <- x11[DNAwheel$other.cluster[i] + 1]
	y2 <- y11[DNAwheel$other.cluster[i] + 1]

	thisRatio <- DNAwheel$mix.tenths[i] / 10 + runif(1, -0.04, 0.04)
	thisRatios[i] <- thisRatio
	
	thisX <- x1 + thisRatio * (x2 - x1)
	thisY <- y1 + thisRatio * (y2 - y1)
	xyCoords[i, 1] <- thisX
	xyCoords[i, 2] <- thisY
	thisColors[i] <- "#3333CC33" # some alpha blue (transparent)
	if ((DNAwheel$mix.tenths[i] > 5) | (DNAwheel$correct.by.seurat[i] == -1)) {
		thisColors[i] <- "#FF330033" # some alpha red
	}
}

lineSegments <- matrix(0, nrow=11 * 11, ncol=4)
for (i in 1:11) {
	for (j in 1:11) {
		thisRow <- (i - 1) * 11 + j
		lineSegments[thisRow, 1] <- x11[i]
		lineSegments[thisRow, 2] <- y11[i]
		lineSegments[thisRow, 3] <- x11[j]
		lineSegments[thisRow, 4] <- y11[j]
	}
}

randomOrder <- sample(1:nrow(xyCoords))

plot(1.3 * x11, 1.3 * y11, type="n")
text(1.2 * x11, 1.2 * y11, labels=clusterNewNames, cex=1.2)
segments(x0=lineSegments[, 1], y0=lineSegments[, 2], x1=lineSegments[, 3], y1=lineSegments[, 4], lwd=0.2, col = "grey50")
points(xyCoords[randomOrder, ], col=thisColors[randomOrder], cex = 0.7, pch=16)

### The END