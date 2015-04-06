# This generates a BAM file tailored to the 'cut.rds' object.

require(GenomicRanges)
cuts <- readRDS("cuts.rds")
ofile <- "hic.sam"

addread <- function(name, fragId, isForward, offset=alen, isFirst=TRUE, isPrimary=TRUE, isDup=FALSE, isUnmapped=FALSE,
		alen=10, hanging=NULL, mapq=200, overrun=FALSE) {
	if (!isPrimary) { offset <- width(cuts[fragId])-offset + alen } # as 3' segments point away from the binding site.
	if (isForward) { 
		actual.pos <- end(cuts[fragId]) - offset + 1L
		if (!overrun) { stopifnot(actual.pos >= start(cuts[fragId])) }
	} else {
		actual.pos <- start(cuts[fragId]) + offset - alen
		if (!overrun) { stopifnot(actual.pos + alen - 1L <= end(cuts[fragId])) }
	}

	# Multiple hangings to allow for tripartite reads.
	cigared <- paste0(alen, "M")
	if (length(hanging)) {
		for (i in 1:length(hanging)) { 
			if ((i==1)==(isPrimary==isForward)) { 
				cigared <- paste0(cigared, hanging[i], "H")
			} else {
				cigared <- paste0(hanging[i], "H", cigared)
			}
		}
	}
	write(file=ofile, append=TRUE, paste(name, 1+ifelse(isUnmapped, 4, 0)+ifelse(isForward, 0, 16)+ifelse(isFirst, 64, 128)+ifelse(isPrimary, 0, 256)+ifelse(isDup, 1024, 0),
		as.character(seqnames(cuts[fragId])), actual.pos, mapq, cigared, "*", 0, 0, paste(rep("N", alen), collapse=""), 
		paste(rep("h", alen), collapse=""), sep="\t"), ncol=1)
	invisible(NULL)
}

####################################################################################################
# All right, generating the file to specifications.

write(file=ofile, paste("@SQ", paste0("SN:", seqlevels(cuts)), paste0("LN:", seqlengths(cuts)), sep="\t"), ncol=1)

# Adding some good read pairs.

addread("good.1", 2, isForward=FALSE, offset=40, isFirst=TRUE)
addread("good.1", 1, isForward=TRUE, offset=40, isFirst=FALSE, mapq=20)

addread("good.2", 2, isForward=TRUE, offset=30, isFirst=TRUE)
addread("good.2", 1, isForward=FALSE, offset=30, isFirst=FALSE)

addread("good.3", 4, isForward=FALSE, offset=25, isFirst=TRUE)
addread("good.3", 1, isForward=FALSE, offset=15, isFirst=FALSE, mapq=20)

addread("good.4", 4, isForward=TRUE, offset=35, isFirst=TRUE)
addread("good.4", 2, isForward=TRUE, offset=45, isFirst=FALSE)

addread("good.5", 5, isForward=TRUE, offset=10, isFirst=TRUE)
addread("good.5", 3, isForward=FALSE, offset=20, isFirst=FALSE)

addread("good.6", 7, isForward=FALSE, offset=20, isFirst=TRUE)
addread("good.6", 5, isForward=FALSE, offset=15, isFirst=FALSE)

# ... and a self-circle.

addread("self.1", 6, isForward=TRUE, offset=15, isFirst=TRUE, isDup=TRUE)
addread("self.1", 6, isForward=FALSE, offset=15, isFirst=FALSE, isDup=TRUE)

# Adding some singleton reads.

addread("singleton.1", 6, isForward=TRUE, offset=20, isFirst=TRUE)
addread("singleton.2", 6, isForward=TRUE, offset=20, isFirst=FALSE)

# Adding instances where the two reads form non-dangling end, non-self-circle pairs on the same fragment.

addread("other.1", 4, isForward=TRUE, offset=20, isFirst=TRUE, isDup=TRUE)
addread("other.1", 4, isForward=TRUE, offset=45, isFirst=FALSE, isDup=TRUE)

addread("other.2", 7, isForward=TRUE, offset=15, isFirst=TRUE)
addread("other.2", 7, isForward=FALSE, offset=15, isFirst=FALSE)

# Adding some dangling ends.

addread("dangling.1", 3, isForward=TRUE, offset=20, isFirst=TRUE, isDup=TRUE)
addread("dangling.1", 3, isForward=FALSE, offset=40, isFirst=FALSE, isDup=TRUE)

addread("dangling.2", 5, isForward=TRUE, offset=20, isFirst=TRUE, mapq=20)
addread("dangling.2", 5, isForward=FALSE, offset=20, isFirst=FALSE)

addread("dangling.3", 7, isForward=TRUE, offset=20, isFirst=TRUE, mapq=20)
addread("dangling.3", 7, isForward=FALSE, offset=20, isFirst=FALSE)

# Checking what happens if the restriction site is overrun.

addread("good.7", 6, isForward=TRUE, offset=5, isFirst=TRUE, overrun=TRUE)
addread("good.7", 5, isForward=FALSE, offset=20, isFirst=FALSE)

addread("good.8", 5, isForward=TRUE, offset=5, isFirst=TRUE, overrun=TRUE)
addread("good.8", 2, isForward=FALSE, offset=5, isFirst=FALSE, overrun=TRUE)

addread("self.2", 4, isForward=FALSE, offset=30, isFirst=TRUE)
addread("self.2", 4, isForward=TRUE, offset=5, isFirst=FALSE, overrun=TRUE)

# Generating chimeras. Alignment length needs to be 5, to get past the overhang.

addread("chimeric.good.1", 6, isForward=TRUE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=5)
addread("chimeric.good.1", 1, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=FALSE, hanging=5, mapq=20)
addread("chimeric.good.1", 1, isForward=TRUE, offset=30, isFirst=FALSE)

addread("chimeric.good.2", 7, isForward=FALSE, offset=15, isFirst=TRUE, isDup=TRUE)
addread("chimeric.good.2", 7, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=FALSE, hanging=5)
addread("chimeric.good.2", 3, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=TRUE, hanging=5, isDup=TRUE)

addread("chimeric.good.3", 7, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=5)
addread("chimeric.good.3", 6, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=FALSE, hanging=5)
addread("chimeric.good.3", 6, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=TRUE, hanging=5)
addread("chimeric.good.3", 7, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=FALSE, hanging=5)

# Invalid because the 3' segment extends past the mate 5' segment.

addread("chimeric.invalid.1", 6, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=10, isDup=TRUE)
addread("chimeric.invalid.1", 3, isForward=FALSE, alen=10, isFirst=TRUE, isPrimary=FALSE, hanging=5)
addread("chimeric.invalid.1", 3, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=TRUE, hanging=5, isDup=TRUE)
addread("chimeric.invalid.1", 6, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=FALSE, hanging=5)

# Invalid because both reads map to multiple locations.

addread("chimeric.invalid.2", 6, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=5)
addread("chimeric.invalid.2", 3, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=FALSE, hanging=5)
addread("chimeric.invalid.2", 2, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=TRUE, hanging=5, mapq=20)
addread("chimeric.invalid.2", 7, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=FALSE, hanging=5, mapq=20)

# Invalid because one read maps to inconsistent locations.

addread("chimeric.invalid.3", 6, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=5)
addread("chimeric.invalid.3", 4, isForward=TRUE, alen=5, isFirst=TRUE, isPrimary=FALSE, hanging=5, mapq=20)
addread("chimeric.invalid.3", 3, isForward=FALSE, offset=15, isFirst=FALSE)

# Invalid because the 3' and mate 5' segments don't form an inward pair.

addread("chimeric.invalid.4", 4, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=5)
addread("chimeric.invalid.4", 2, isForward=TRUE, alen=5, isFirst=TRUE, isPrimary=FALSE, hanging=5, isDup=TRUE)
addread("chimeric.invalid.4", 2, isForward=TRUE, offset=15, isFirst=FALSE)

# Invalid because one read has 3 mapping locations.

addread("chimeric.invalid.5", 2, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=28)
addread("chimeric.invalid.5", 5, isForward=FALSE, alen=23, isFirst=TRUE, isPrimary=FALSE, hanging=c(5,5))
addread("chimeric.invalid.5", 1, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=FALSE, hanging=28)
addread("chimeric.invalid.5", 1, isForward=TRUE, offset=35, isFirst=FALSE)

# Invalid because there's non-specific ligation to a shared fragment.

addread("chimeric.invalid.6", 4, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=5)
addread("chimeric.invalid.6", 3, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=FALSE, hanging=5)
addread("chimeric.invalid.6", 3, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=TRUE, hanging=5, mapq=20)
addread("chimeric.invalid.6", 1, isForward=TRUE, alen=5, isFirst=FALSE, isPrimary=FALSE, hanging=5, mapq=20)

# Invalid because it extends past its mate, again.

addread("chimeric.invalid.7", 6, isForward=FALSE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=20)
addread("chimeric.invalid.7", 3, isForward=TRUE, alen=20, isFirst=TRUE, isPrimary=FALSE, hanging=5)
addread("chimeric.invalid.7", 3, isForward=FALSE, offset=10, isFirst=FALSE)

# Offsite, but not necessarily invalid.

addread("chimeric.good.4", 5, isForward=TRUE, alen=5, offset=10, isFirst=TRUE, isPrimary=TRUE, hanging=5, isDup=TRUE)
addread("chimeric.good.4", 4, isForward=TRUE, alen=5, isFirst=TRUE, isPrimary=FALSE, hanging=5)
addread("chimeric.good.4", 4, isForward=FALSE, offset=20, isFirst=FALSE, isDup=TRUE)

addread("chimeric.good.5", 5, isForward=TRUE, alen=5, isFirst=TRUE, isPrimary=TRUE, hanging=5)
addread("chimeric.good.5", 4, isForward=FALSE, alen=5, offset=10, isFirst=TRUE, isPrimary=FALSE, hanging=5)
addread("chimeric.good.5", 4, isForward=TRUE, offset=20, isFirst=FALSE, isDup=TRUE)

# Read pairs with one unmapped component.

addread("unmap.1", 4, isForward=FALSE, offset=20, isFirst=TRUE, isUnmapped=TRUE)
addread("unmap.1", 1, isForward=TRUE, offset=30, isFirst=FALSE)

addread("unmap.2", 5, isForward=TRUE, offset=10, isFirst=TRUE)
addread("unmap.2", 3, isForward=FALSE, offset=25, isFirst=FALSE, isUnmapped=TRUE)

addread("unmap.3", 6, isForward=FALSE, offset=20, isFirst=TRUE, isUnmapped=TRUE)
addread("unmap.3", 5, isForward=FALSE, offset=15, isFirst=FALSE, isUnmapped=TRUE)

####################################################################################################
## Converting it to a BAM file.

require(Rsamtools)
asBam(ofile, destination="hic", overwrite=TRUE, indexDestination=FALSE)
sortBam("hic.bam", destination="hic_sort", byQname=TRUE)
unlink(c("hic.bam", "hic.sam"))

# require(diffHic)
# preparePairs(bam="hic_sort.bam", dir="stuff", fragments=cuts, dedup=FALSE)
# getFragmentData("stuff/")

