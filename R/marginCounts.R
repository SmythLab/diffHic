marginCounts <- function(files, param, width=50000)
# Gets the marginal counts i.e. sum of counts for each bin or region.
# This is useful to determine the `genomic coverage' of each region,
# based on the number of Hi-C read pairs involving that region.
#
# written by Aaron Lun
# some time ago.
# last modified 22 November 2015
{
	nlibs <- length(files)
	width <- as.integer(width)
	fragments <- param$fragments
	frag.by.chr <- .splitByChr(fragments)

	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	if (width < 0) { stop("width must be a non-negative integer") }
	new.pts <- .getBinID(fragments, width)
	bin.by.chr <- .splitByChr(new.pts$region)

	total.counts <- matrix(0L, length(new.pts$region), nlibs)
	full.sizes <- integer(nlibs)
	chrs <- seqlevels(fragments)

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files, chrs, restrict)
	for (anchor1 in names(overall)) {
		current <- overall[[anchor1]]
		first.anchor1 <- bin.by.chr$first[[anchor1]]
		last.anchor1 <- bin.by.chr$last[[anchor1]]
		nabins <- last.anchor1 - first.anchor1 + 1L
		keep.a <- first.anchor1:last.anchor1

		for (anchor2 in names(current)) {
			first.anchor2 <- bin.by.chr$first[[anchor2]]
			last.anchor2 <- bin.by.chr$last[[anchor2]]
			ntbins <- last.anchor2 - first.anchor2 + 1L
			keep.t <- first.anchor2:last.anchor2

			pairs <- .baseHiCParser(current[[anchor2]], files, anchor1, anchor2,
				chr.limits=frag.by.chr, discard=discard, cap=cap)

			# Aggregating them for each library.
			for (lib in seq_len(nlibs)) {
				a.counts <- tabulate(new.pts$id[pairs[[lib]]$anchor1.id]-first.anchor1+1L, nbins=nabins)
				total.counts[keep.a,lib] <- total.counts[keep.a,lib] + a.counts
				t.counts <- tabulate(new.pts$id[pairs[[lib]]$anchor2.id]-first.anchor2+1L, nbins=ntbins)
				total.counts[keep.t,lib] <- total.counts[keep.t,lib] + t.counts
				full.sizes[lib] <- full.sizes[lib] + nrow(pairs[[lib]])
			}
		}	
	}
	
	# Aggregating all elements.
	return(SummarizedExperiment(list(counts=total.counts), colData=DataFrame(totals=full.sizes), 
			rowRanges=new.pts$region, metadata=List(param=param)))
}

