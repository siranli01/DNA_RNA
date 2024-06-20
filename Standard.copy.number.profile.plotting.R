######
### Standard 300-bin copy-number-profile plotting 
### input file is a 300-element vector "vector" representing the 300 template counts in each bin 
### bin-information file "bins.boundaries.txt" is in the main folder

vector <- as.matrix(vector) ### length(vector) == 300
suffix <- seq(1:300)
bin.names <- paste("b", suffix, sep="")
rownames(vector) <- bin.names

c.varbin <- read.table("/.../bins.boundaries.txt", sep="\t", header=T, as.is=T, stringsAsFactors=F)
chrom.starts <- c(which(c.varbin[, 1] != c(c.varbin[, 1][-1], c.varbin[, 1][300])) + 1, 301)

CN.results <- cbs.segment.uber01.0_3k(outdir="/ any directory", indir="", varbin.gc="/.../bins.boundaries.txt", thisUber=vector, abspos.col=c.varbin$abspos, alpha=0.1, nperm=1000, undo.SD=0.25, min.width=2)

CN.ratio.mat <- CN.results[[1]]
CN.seg.mat <- CN.results[[2]]

### plot standard copy-number profile

asdf <- CN.seg.mat
asdf.ratio <- CN.ratio.mat

chr <- c.varbin$chr
chr.shift <- c(chr[-1], chr[length(chr)])
vlines <- c(1, c.varbin$abspos[which(chr != chr.shift) + 1], c.varbin$abspos[nrow(asdf)])
hlines <- c(0.5, 1, 2, 3, 4, 5, 6, 7)
chr.text <- c(1:10,"",12,"",14,"",16,"",18,"","","","", "X")
vlines.shift <- c(vlines[-1], 4*10^9)
chr.at <- vlines + (vlines.shift - vlines) / 2
x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

qRatio <- matrix(, nrow=nrow(asdf), ncol=ncol(asdf))
qSeg <- matrix(, nrow=nrow(asdf), ncol=ncol(asdf))

qSeg <- asdf

for (i in 1:ncol(asdf)) {
	work.seg <- asdf[,i]
	thisGrid <- seq(1.5, 4.0, by=0.05)
	thisOuter <- work.seg %o% thisGrid
	thisOuterRound <- round(thisOuter)
	thisOuterDiff <- (thisOuter - thisOuterRound) ^ 2
	thisOuterColsums <- colSums(thisOuterDiff, na.rm = FALSE, dims = 1)
	thisMultiplier <- thisGrid[which.min(thisOuterColsums)]
	thisError <- min(thisOuterColsums)
	qRatio[,i] <- asdf.ratio[,i] * thisMultiplier
	qSeg[,i] <- asdf[,i] * thisMultiplier
}


xmax <- c.varbin$abspos[nrow(c.varbin)]

plot(x=c.varbin$abspos, y=qRatio[, 1], log="y", main=thisCell, cex.main = 0.8, xaxt="n", xlab="Genome Position (Gb)", ylab="Copy Number", col="grey70", ylim=c(0.2, 12.0), xlim=c(0, xmax),  mgp=c(2.0,0.5,0), cex=0.4, pch=16)
axis(1, at=x.at, labels=x.labels)
lines(x=c.varbin$abspos, y=qRatio[, 1], col="grey85", cex=0.1)

points(x=c.varbin$abspos, y=qSeg[, 1], col="darkred", cex=1, pch=16)
lines(x=c.varbin$abspos, y=qSeg[, 1], col="darkred", cex=0.2)

abline(h=hlines, col="#666600", lwd=0.2, lty=3, xlim=c(0, xmax))
abline(v=vlines, col="grey40", lwd=0.2)
mtext(chr.text, at = chr.at, line=-1, cex=0.8)


### DNA segmentation functions to run first

library("DNAcopy", lib.loc="/.../DNAcopy_1.50.1")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

cbs.segment.uber01.0_3k <- function(outdir, indir, varbin.gc, thisUber, abspos.col, alpha, nperm, undo.SD, min.width) {

	gc <- read.table(varbin.gc, header=T)

	chrom.numeric <- substring(gc$chr, 4)
	chrom.numeric[which(gc$chr == "chrX")] <- "23"
	chrom.numeric[which(gc$chr == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisUberRatio <- thisUber
#	thisUberRatioQuantal <- thisUber
	thisUberSeg <- thisUber
#	thisUberSegQuantal <- thisUber

	for (i in 1:ncol(thisUber)) {
		sample.name <- dimnames(thisUber)[[2]][i]
		cat(i, dimnames(thisUber)[[2]][i])
		
		thisBincount <- thisUber[, i] + 1
		thisRatio <- thisBincount / mean(thisBincount)
		thisLowratio <- lowess.gc(gc$gc, thisRatio)
		
		thisUberRatio[, i] <- thisLowratio
		
		set.seed(25) 
		CNA.object <- CNA(log(thisLowratio), chrom.numeric, gc$start, data.type="logratio", sampleid=dimnames(thisUber)[[2]][i]) 
		smoothed.CNA.object <- smooth.CNA(CNA.object) 
		segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=2) 
		thisShort <- segment.smoothed.CNA.object[[2]]

		m <- matrix(data=0, nrow=nrow(thisUber), ncol=1)	
		prevEnd <- 0
		for (j in 1:nrow(thisShort)) {
			thisStart <- prevEnd + 1
			thisEnd <- prevEnd + thisShort$num.mark[j]
			m[thisStart:thisEnd, 1] <- exp(thisShort$seg.mean[j])
			prevEnd = thisEnd
		}

		thisUberSeg[, i] <- m[, 1]
		write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.0_3k.k50.varbin.short.txt", sep=""), quote=F, row.names=F) 

		chr <- chrom.numeric
		chr.shift <- c(chr[-1], chr[length(chr)])
		vlines <- c(1, abspos.col[which(chr != chr.shift) + 1], abspos.col[nrow(thisUber)])
		hlines <- c(0.5, 1.0, 1.5, 2.0)
		chr.text <- c(1:22, "X", "Y")
		vlines.shift <- c(vlines[-1], 4*10^9)
		chr.at <- vlines + (vlines.shift - vlines) / 2
		x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
		x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")

		png(paste(outdir, "/", sample.name, ".0_3k.wg.png", sep=""), height=800, width=1200)
		plot(x=abspos.col, y=thisLowratio, log="y", main=paste(sample.name, ""), xaxt="n", xlab="Genome Position Gb", ylab="Ratio", col="#CCCCCC", ylim=c(0.05, 10.0))
		axis(1, at=x.at, labels=x.labels)
		lines(x=abspos.col, y=thisLowratio, col="#CCCCCC")
		points(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		lines(x=abspos.col, y=thisUberSeg[, i], col="#0000AA")
		abline(h=hlines)
		abline(v=vlines)
		mtext(chr.text, at = chr.at)
		dev.off()
		
		thisRatioOut <- data.frame(gc[, 1:3], "bincount"=thisBincount, "ratio" = thisRatio, "gc" = gc$gc, "lowratio" = thisLowratio, "seg.mean.LOWESS" = thisUberSeg[, i])
		write.table(thisRatioOut, sep="\t", file=paste(outdir, "/", sample.name, ".hg19.0_3k.k50.varbin.data.txt", sep=""), quote=F, row.names=F)
		
	}
	
	return(list(thisUberRatio, thisUberSeg))
}

### END









