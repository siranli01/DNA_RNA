######
### Alluvial Diagram plotting, using Tumor2 as an example. 
### Using DNA Seurat object "Tumor2_DNA" and RNA Seurat object "Tumor2_RNA" as input files.
### "Tumor2_DNA" and "Tumor2_RNA" are from the hybrid protocol, so they contain the same group of nuclei barcodes

library(alluvial)

DimPlot(Tumor2_DNA, reduction = "umap", label = TRUE, label.size = 4, pt.size=0.5)
DNA_N <- CellsByIdentities(Tumor2_DNA, idents = 'N')$'N' 
DNA_Nx <- CellsByIdentities(Tumor2_DNA, idents = 'Nx')$'Nx' 
DNA_T2A <- CellsByIdentities(Tumor2_DNA, idents = 'T2A')$'T2A'
DNA_T2B <- CellsByIdentities(Tumor2_DNA, idents = 'T2B')$'T2B' 
DNA_T2C <- CellsByIdentities(Tumor2_DNA, idents = 'T2C')$'T2C' 
DNA_T2D <- CellsByIdentities(Tumor2_DNA, idents = 'T2D')$'T2D'


DimPlot(Tumor2_RNA, reduction = "umap", label = TRUE, label.size = 4, pt.size=0.3)
nclusters <- nlevels(Tumor2_RNA@meta.data$seurat_clusters)

temp <- CellsByIdentities(Tumor2_RNA, idents = 'RNA2a')
names(temp) <- "bc"
RNA.clusters <- data.frame("bc"=temp, "g"="0", stringsAsFactors=FALSE)
temp <- CellsByIdentities(Tumor2_RNA, idents = 'RNA2b')
names(temp) <- "bc"
RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="1", stringsAsFactors=FALSE))
if (nclusters > 2) {
	temp <- CellsByIdentities(Tumor2_RNA, idents = 'Mono')
	names(temp) <- "bc"
	RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="3", stringsAsFactors=FALSE))
}
if (nclusters > 3) {
	temp <- CellsByIdentities(Tumor2_RNA, idents = 'Lymphs')
	names(temp) <- "bc"
	RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="2", stringsAsFactors=FALSE))
}
if (nclusters > 4) {
	temp <- CellsByIdentities(Tumor2_RNA, idents = 'Fibroblast')
	names(temp) <- "bc"
	RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="4", stringsAsFactors=FALSE))
}
if (nclusters > 5) {
	temp <- CellsByIdentities(Tumor2_RNA, idents = 'EC')
	names(temp) <- "bc"
	RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="5", stringsAsFactors=FALSE))
}
if (nclusters > 6) {
	temp <- CellsByIdentities(Tumor2_RNA, idents = 'EP')
	names(temp) <- "bc"
	RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="6", stringsAsFactors=FALSE))
}
if (nclusters > 7) {
	temp <- CellsByIdentities(Tumor2_RNA, idents = 'Plasma cells')
	names(temp) <- "bc"
	RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="7", stringsAsFactors=FALSE))
}
if (nclusters > 8) {
	temp <- CellsByIdentities(Tumor2_RNA, idents = 8)
	names(temp) <- "bc"
	RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="8", stringsAsFactors=FALSE))
}
if (nclusters > 9) {
	temp <- CellsByIdentities(Tumor2_RNA, idents = 9)
	names(temp) <- "bc"
	RNA.clusters <- rbind.data.frame(RNA.clusters, data.frame("bc"=temp, "g"="9", stringsAsFactors=FALSE))
}

Cluster.info <- data.frame("bc"=intersect(DNA_T2A,colnames(Tumor2_RNA)), "DNA.cluster"="A", stringsAsFactors=FALSE)
Cluster.info <- rbind.data.frame(Cluster.info, data.frame("bc"=intersect(DNA_T2B,colnames(Tumor2_RNA)), "DNA.cluster"="B", stringsAsFactors=FALSE))
Cluster.info <- rbind.data.frame(Cluster.info, data.frame("bc"=intersect(DNA_T2C,colnames(Tumor2_RNA)), "DNA.cluster"="C", stringsAsFactors=FALSE))
Cluster.info <- rbind.data.frame(Cluster.info, data.frame("bc"=intersect(DNA_T2D,colnames(Tumor2_RNA)), "DNA.cluster"="D", stringsAsFactors=FALSE))
Cluster.info <- rbind.data.frame(Cluster.info, data.frame("bc"=intersect(DNA_Nx,colnames(Tumor2_RNA)), "DNA.cluster"="E", stringsAsFactors=FALSE))
Cluster.info <- rbind.data.frame(Cluster.info, data.frame("bc"=intersect(DNA_N,colnames(Tumor2_RNA)), "DNA.cluster"="F", stringsAsFactors=FALSE))

Cluster.info$RNA.cluster <- RNA.clusters$g[match(Cluster.info$bc, RNA.clusters$bc, nomatch=NA)]

### Build connections between DNA and RNA identities

x <- Cluster.info[, c("DNA.cluster", "RNA.cluster")]
table(x)

xmat <- matrix(0, nrow=(length(unique(x$DNA.cluster)) * length(unique(x$RNA.cluster))), ncol=3)
xmat[, 1] <- rep(0:5, each=8)
xmat[, 2] <- rep(0:7, 6)
xmat[, 3] <- as.vector(t(table(x)))
xdf <- data.frame(xmat)
names(xdf) <- c("DNA.cluster", "RNA.cluster", "freq")

x2d <- aggregate( freq ~ DNA.cluster + RNA.cluster, data=xdf, sum)

x2d$DNA.cluster <- factor(x2d$DNA.cluster) 
levels(x2d$DNA.cluster) <- c("T2A", "T2B","T2C","T2D","Nx","N")
x2d$RNA.cluster <- factor(x2d$RNA.cluster) 

levels(x2d$RNA.cluster) <- c("RNA2a", "RNA2b", "Lymphs","Mono","Fibroblast","EC","EP","Plasma cells")

alluvial( x2d[,1:2], freq=x2d$freq, xw=0.2, alpha=0.8, gap.width=0.1, col=c("#00BA38", "#00BFC4","#F564E3","#4155FF","#B79F00", "#F8766D"), border="grey90", cex.axis=1.2, cex=0.9)

### END


