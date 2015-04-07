# This script converts the BAM files to index files for processing within diffHic.
# Note, if a custom FASTA file was used to build the genome indices for alignment, then 
# the path to that FASTA file should be used as the argument to cutGenome.

require(diffHic)
require(BSgenome.Mmusculus.UCSC.mm10)
require(BSgenome.Hsapiens.UCSC.hg19)

for (i in 1:2) {
	if (i==1L) {
		# Rickman.
		bam <- c("SRR493818.bam", "SRR493819.bam", "SRR493820.bam", "SRR493821.bam", "SRR493822.bam", "SRR493823.bam", "SRR493824.bam", "SRR493825.bam", "SRR493826.bam", "SRR493827.bam")
		fragments <- cutGenome(BSgenome.Hsapiens.UCSC.hg19, pattern="AAGCTT", overhang=4)
	} else if (i==2L) { 
		# Sofueva.
		bam <- c("SRR941267.bam", "SRR941268.bam", "SRR941269.bam", "SRR941270.bam", "SRR941271.bam", "SRR941272.bam", "SRR941273.bam", "SRR941274.bam", "SRR941275.bam", "SRR941276.bam", "SRR941277.bam", "SRR941278.bam", "SRR941279.bam", "SRR941280.bam", "SRR941281.bam", "SRR941282.bam")
		fragments <- cutGenome(BSgenome.Mmusculus.UCSC.mm10, pattern="AAGCTT", overhang=4)
	}

	param <- pairParam(fragments=fragments)
	for (bf in bam) {
		outfile <- sub("\\.bam$", ".h5", bf)
		preparePairs(bf, param=param, file=outfile, dedup=TRUE, minq=10)
		prunePairs(outfile, param=param, min.inward=1000, min.outward=25000, max.frag=600)
	}

	# Merging technical replicates.
	if (i==1L) { 
		mergePairs(c("SRR493820.h5", "SRR493821.h5", "SRR493822.h5", "SRR493823.h5"), "merged_erg.h5")
		mergePairs(c("SRR493824.h5", "SRR493825.h5", "SRR493826.h5", "SRR493827.h5"), "merged_gfp.h5")
	} else if (i==2L) {
		mergePairs(c("SRR941267.h5", "SRR941268.h5", "SRR941269.h5", "SRR941270.h5"), "merged_flox_1.h5")
		mergePairs(c("SRR941271.h5", "SRR941272.h5", "SRR941273.h5", "SRR941274.h5"), "merged_flox_2.h5")
		mergePairs(c("SRR941275.h5", "SRR941276.h5", "SRR941277.h5", "SRR941278.h5"), "merged_ko_1.h5")
		mergePairs(c("SRR941279.h5", "SRR941280.h5", "SRR941281.h5", "SRR941282.h5"), "merged_ko_2.h5")
	}
}

