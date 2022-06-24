### DNA clustering based on copy number with sliding windows
### input file is the normalized bin counts matrix "norm_bin_counts"
### Genomic bin information is in the main folder named "no_gene.avoid50closeTN.bins.txt"

norm_bin_counts <- read.table("/.../norm_bin_counts.txt", header=T, as.is=T, stringsAsFactors=F)

suffix <- seq(1:300)
bin.names <- paste("b", suffix, sep="")
data_matrix <- as.matrix(norm_bin_counts)
rownames(data_matrix) <- bin.names
data_matrix_norm <- normalize.matrix(data_matrix,2)
colnames(data_matrix_norm) <- paste0(colnames(data_matrix_norm),'_2')

c.varbin <- read.table("/.../no_gene.avoid50closeTN.bins.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
chrom.starts <- c(which(c.varbin[, 1] != c(c.varbin[, 1][-1], c.varbin[, 1][300])) + 1, 301)

rollsum.median <- function (x, i) {
  xm <- matrix(0, nrow=length(x), ncol=i)
  xm[, 1] <- x
  for (j in 1:(i-1)) {
    xm[, (j+1)] <- c(x[-(1:j)], rep(x[length(x)], j))
  }
  return(apply(xm, 1, median))
}

combined_matrix <- data_matrix_norm
for (j in 2:3) {
	dfi <- apply(data_matrix_norm, 2, rollsum.median, i=j)
	remove.indexes <- c()
	for (k in 1:(j - 1)) {
		remove.indexes <- c(remove.indexes, chrom.starts - k)
	}
	cat(nrow(dfi), length(remove.indexes), length(unique(remove.indexes)), "\n")
	combined_matrix <- rbind(combined_matrix, dfi[-remove.indexes, ])
}

### Seurat clustering based on copy number

rownames(combined_matrix) <- paste("b", 1:nrow(combined_matrix), sep="")

combined_matrix.a <- CreateSeuratObject(counts = combined_matrix, project = "test", min.cells = 30, min.features = 200) 
combined_matrix.a <- FindVariableFeatures(combined_matrix.a, selection.method = "vst", nfeatures = 500)
combined_matrix.a <- ScaleData(combined_matrix.a, verbose = FALSE)
combined_matrix.a <- RunPCA(combined_matrix.a, npcs = 50, verbose = FALSE)
ElbowPlot(combined_matrix.a, ndims = 50)
combined_matrix.a <- RunUMAP(combined_matrix.a, reduction = "pca", dims = 1:30)
combined_matrix.a <- FindNeighbors(combined_matrix.a, reduction = "pca", dims = 1:30)
combined_matrix.a <- FindClusters(combined_matrix.a, resolution = 0.4)
DimPlot(combined_matrix.a, reduction = "umap", pt.size=0.5)

### END

