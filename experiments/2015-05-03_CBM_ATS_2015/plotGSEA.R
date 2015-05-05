plot.gsea <- function(gene.rank.file, gene.set.files, low="blue",
			    mid="white", high="red", num.bins=1000, gene.set.lwd=2,
			    max.enrich=1, num.rows=NULL, gene.set.scale=1,
			    gene.rank.scale=1, gene.set.cex=1, gene.rank.cex=1,
                      plot.labels=T, lab.rank=NULL, lab.gene.set=NULL) {


######## Arguments ############################################################
#
#	gene.rank.file:	Ranked gene list in the edb folder from the GSEA results
#	
#	gene.set.files:	A list of files corresponding to the enrichment result
#				tables.  Each file should have the gene names in the 
#				second column and the enrichment score in the seventh
#				column.  This table can be found by opening the
#				HTML file called "index.html" in the GSEA results
#				folder and then opening the detailed results for each
#				gene set.  Copy and paste the table into notepad without
#				the column headers.
#
#	          low:	The lowest color for the negative values in the ranked
#				gene list.
#
#	          mid:	The middle color for the ranked gene list.
#
#	         high:	The highest color for positive values in the ranked
#				gene list.
#
#	     num.bins:	The number of bins in the color bar to draw.
#
#	 gene.set.lwd:	The width of the gene set lines.
#
#	   max.enrich:	The maximum enrichement score to draw.  If an enrichment
#				score is above max.enrich, the line will be cut off at
#				max.enrich.  If the enrichment scores are far below
#				max.enrich, then there may be a lot of extra white space
#				that can be eliminated by lowering max.enrich.
#
#	     num.rows:	Set the number of items (e.g. rows) to be plotted. The
#				default is the number of gene sets plus one for the
#				ranked gene	list.  If you want the enrichment score
#				lengths to be consistent between GSEAs with different
#				numbers of gene sets, you may need to set this number
#				manually to the maximum number of gene sets between all
#				your GSEAs.
#
#    gene.set.scale:	Scales the gene set enrichment scores.  This allows you
#				to tweek the heights of the bars.
#
#   gene.rank.scale:	Scales the color bar for the ranked gene list.  This
#				allows you to tweek the heights of the color bar.
#
#      gene.set.cex:	Scales the gene set labels.
#
#	gene.rank.cex:	Scales the ranked gene list label.
#
#	  plot.labels:	logical. Whether or not to plot the labels.
#
#	     lab.rank:	The label for the color bar.  Defaults to the file name
#				of the ranked gene list.
#
#      lab.gene.set:	A vector of labels for the gene sets.  Defaults to the
#				vector of file names for the gene sets.
#
###############################################################################

	## Load Gplots library
	library(gplots, verbose=F)

	## Read in the ranked gene list from GSEA
	genes <- read.table(gene.rank.file, stringsAsFactors=F, sep="\t", quote="")
	gene.rank <- genes[,2]
	names(gene.rank) <- genes[,1]
	
	## Sort ranked genes
	gene.rank <- gene.rank[order(gene.rank, decreasing=F)]

	## Set up the labels if none were provided
	if(is.null(lab.rank) & plot.labels) {
		lab.rank <- gene.rank.file
	}
	if(is.null(lab.gene.set) & plot.labels) {
		lab.gene.set <- gene.set.files
	}

	## Set up the color bar
	num.genes <- length(gene.rank)
	pos.col <- colorpanel(sum(gene.rank > 0), low=mid, high=high)
	neg.col <- colorpanel(sum(gene.rank <= 0), low=low, high=mid)
	col <- c(neg.col, pos.col)

	reference.seq.pos <- seq(0, max(gene.rank), length.out=sum(gene.rank > 0))
	reference.seq.neg <- seq(min(gene.rank), 0, length.out=sum(gene.rank <= 0))
	reference.seq <- c(reference.seq.neg, reference.seq.pos)

	## Set the number of items (e.g. rows) to be plotted
	## The default is the number of gene sets plus the ranked
	## gene list.  If you want the enrichment score lengths to be
	## consistent between GSEAs with different numbers of gene sets,
	## you may need to set this number manually
	if(is.null(num.rows)) {
		num.rows <- length(gene.set.files) + 1
	}
	
	## Set up layout
	layout(matrix(num.rows:1, ncol=1, nrow=num.rows))

	## Set up margins
	bottom.mar <- 0
	if(plot.labels) {
		bottom.mar <- 2
	}

	## Set up blank plot for the ranked gene list
	par(mar=c(bottom.mar,0,0,0), oma=c(0,0,0,0))
	plot(0, xlim=c(0,1), ylim=c(0,1), type="n", axes=F, ann=F, frame.plot=F)

	## Plot gene rank label
	if(plot.labels) {
		mtext(lab.rank, side=1, cex=gene.rank.cex)
	}

	## Calculate bins
	my.bins <- round(seq(1, length(gene.rank), length.out=num.bins)) 

	## Plot the color bar using a rectangle for each bin
	for (i in my.bins) {
		col.ind <- which.min(abs(reference.seq - gene.rank[i]))
		percent <- i/num.genes
		rect(percent, gene.rank.scale*1, percent + (1/num.bins), 0, col=col[col.ind], border=col[col.ind]) 
	}
 
	## Plot the gene sets
	for (i in length(gene.set.files):1) {

		## Set up blank plot for gene set
		par(mar=c(bottom.mar,0,0,0))
		plot(0, xlim=c(0,1), ylim=c(0,max.enrich), type="n", ann=F, axes=F, frame.plot=F)

		## Plot gene set label
		if(plot.labels) {
			mtext(lab.gene.set[i], side=1, line=0, cex=gene.set.cex)		
		}

		## Read in enrichment scores from file
		enrich <- read.table(gene.set.files[i], stringsAsFactors=F, sep="\t")

		## Plot each enrichment score
		for (j in 1:nrow(enrich)) {
			index <- which(names(gene.rank) == enrich[j,2]) / num.genes
			enrich.score <- gene.set.scale*abs(enrich[j,7])
			lines(c(index, index), c(0, enrich.score), lwd=gene.set.lwd)
		}
	}
}
