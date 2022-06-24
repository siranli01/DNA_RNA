# DNA_RNA
# R code to create figures for the manuscript "High-throughput single-nucleus hybrid sequencing reveals genome-transcriptome correlations in cancer"

### Clusteing.based.on.copy.number.with.sliding.window.R
It performs DNA clustering based on copy number. 
Input file is the bin-count matrix.
### Standard.copy.number.profile.plotting.R
It plots copy number profile with segmentation information using DNAcopy (in the folder).
Input file is a 300-element bin-count vector.
### Alluvial.diagiam.hybrid.protocol.R
It plots the alluvial plots connecting DNA clustering and RNA clustering results from the new hybrid protocol.
Input files are the DNA and RNA Seurat objects with clustering information.
### Multinomial.wheel.analyses.code.R
It performs multinomial wheel analyses which quantify the similarity of every cell to the major clusters.
Input files are the original bin-counts or gene-count matrix and the Seurat object with clustering information.
### bins.boundaries.txt
The detailed genomic bin information for the hybrid protocol. It is required for copy-number profile plotting.
