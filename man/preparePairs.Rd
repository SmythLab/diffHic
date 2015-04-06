\name{preparePairs}
\alias{preparePairs}

\title{Prepare Hi-C pairs}

\description{Identifies the interacting pair of restriction fragments corresponding to each read pair in a Hi-C library.}

\usage{
preparePairs(bam, param, file, dedup=TRUE, yield=1e7, ichim=TRUE, minq=NA)
}

\arguments{
	\item{bam}{a character string containing the path to a name-sorted BAM file}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{file}{a character string specifying the path to an output index file}
	\item{dedup}{a logical scalar indicating whether marked duplicate reads should be removed}
	\item{yield}{a numeric scalar specifying the number of reads to extract at every iteration}
	\item{ichim}{a logical scalar indicating whether invalid chimeras should be counted}
	\item{minq}{an integer scalar specifying the minimum mapping quality for each read}
}

\section{Converting to restriction fragment indices}{
The resolution of a Hi-C experiment is defined by the distribution of restriction sites across the genome. 
Thus, it makes sense to describe interactions in terms of restriction fragments. 
This function identifies the interacting fragments corresponding to each pair of reads in a Hi-C library. 
To save space, it stores the indices of the interacting fragments for each read pair, rather than the fragments themselves.

Indexing is performed by matching up the mapping coordinates for each read with the restriction fragment boundaries in \code{param$fragments}.
Needless to say, the boundary coordinates in \code{param$fragments} must correspond to the reference genome being used. 
In most cases, these can be generated using the \code{\link{cutGenome}} function from any given \code{BSgenome} object. 
If, for any reason, a modified genome is used for alignment, then the coordinates of the restriction fragments on the modified genome are required.

Each read pair subsequently becomes associated with a pair of restriction fragments.
The anchor fragment is that with the higher genomic coordinate, i.e., the larger index in \code{param$fragments}.
The target fragment is that with the smaller coordinate/index. 
This definition avoids the need to consider both permutations of indices in a pair during downstream processing.
}

\section{Handling of read pairs}{
For pairs with chimeric reads, the alignment of the 5' end of the read is used to identify the interacting fragments.  
Invalid chimeras arise when the position of the 3' segment of a chimeric read is not consistent with that of the mate read.  
These are generally indicative of mapping errors but can also form due to non-specific ligation events. 
While they can be explicitly removed, setting \code{ichim=TRUE} is recommended to avoid excessive filtering of reads when alignment of short chimeric segments is inaccurate.

Self-circles are outward-facing read pairs mapped to the same restriction fragment.
These are formed from inefficient cross-linking and are generally uninformative.
Dangling ends are inward-facing read pairs mapped to the same fragment, and are generated from incomplete ligation of blunt ends.
Both constructs are detected and discarded within the function.

% Also note that overlapping pairs of self-circles are retained, as are dangling ends where one end extends past the other.
% This because they are neither valid self-circles or dangling ends.

In all cases, if the 5' end of either read is unavailable (e.g. unmapped, mapping quality score below \code{minq}, marked as a duplicate), the read pair is recorded as unpaired and discarded.
By default, no MAPQ filtering is performed when \code{minq} is set to \code{NA}.
Any duplicate read must be marked in the bit field of the BAM file using a tool like Picard's \code{MarkDuplicates} if it is to be removed with \code{dedup=TRUE}. 
For chimeric reads, the recommended approach is to designate the 5' end as the primary alignment. 
This ensures that the duplicate computations are performed on the most relevant alignments for each read pair.
}

\section{Miscellaneous information}{
The value of \code{yield} simply controls the amount of memory used by regulating the number of read pairs that are loaded in at any given time. 
Large values will require more memory but will reduce the computation time. 
Input files must be sorted by read name for proper execution of the function.
Read pair information need not be synchronised so long as the pairing, first read and second read flags are appropriately set.

Users should note that the use of a \code{pairParam} object for input is strictly for convenience.
Any non-empty values of \code{param$discard} and \code{param$restrict} will be ignored here.
Reads will not be discarded if they lie outside the specified chromosomes, or if they lie within blacklisted regions.
}

\value{
Multiple dataframe objects are stored within the specified \code{file} using the HDF5 format.
Each object corresponds to a pair of anchor/target chromosomes.
Each row of the dataframe contains information for a read pair, with one read mapped to each chromosome. 
The dataframe contains several integer fields:
\describe{
  \item{\code{anchor.id}:}{index of the anchor restriction fragment}
  \item{\code{target.id}:}{index of the target restriction fragment}
  \item{\code{anchor.pos}:}{mapping coordinate of the read (or for chimeras, the 5' segment thereof) on the anchor fragment}
  \item{\code{target.pos}:}{mapping coordinate of read on the target fragment}
  \item{\code{anchor.len}:}{length of the alignment on the anchor fragment, set to a negative value if the alignment is on the reverse strand}
  \item{\code{target.len}:}{length of the alignment on the target fragment, set to a negative value if reverse stranded}
}
See \code{\link{getPairData}} for more details on \code{length}, \code{orientation} and \code{gap}.  

An integer vector is also returned from the function, containing various diagnostics:
\describe{
	\item{\code{pairs}:}{an integer vector containing \code{total}, the total number of read pairs; \code{marked}, read pairs with at least one marked 5' end; \code{filtered}, read pairs lost due to low MAPQ scores for the 5' end(s); \code{mapped}, read pairs considered as successfully mapped (i.e., not marked or filtered)}
	\item{\code{same.id}:}{an integer vector containing \code{dangling}, the number of read pairs that are dangling ends; and \code{self.circles}, the number of read pairs forming self-circles}
	\item{\code{singles}:}{an integer scalar specifying the number of reads without a mate}
	\item{\code{chimeras}:}{an integer vector containing \code{total}, the total number of read pairs with one chimeric read; \code{mapped}, chimeric read pairs with both 5' ends mapped; \code{multi}, mapped chimeric pairs with at least one successfully mapped 3' segment; and \code{invalid}, read pairs where the 3' location of one read disagrees with the 5' location of the mate}
}
}

\author{Aaron Lun}

\seealso{
\code{\link{cutGenome}}, 
\code{\link{prunePairs}}, 
\code{\link{mergePairs}}, 
\code{\link{getPairData}}
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts) 

tmpf <- "gunk.h5"
preparePairs(hic.file, param, tmpf)
preparePairs(hic.file, param, tmpf, minq=50)
preparePairs(hic.file, param, tmpf, ichim=TRUE)
preparePairs(hic.file, param, tmpf, dedup=FALSE)

\dontshow{
unlink(tmpf)
}
}

\references{
Imakaev M et al. (2012). Iterative correction of Hi-C data reveals hallmarks of chromosome organization. \emph{Nat. Methods} 9, 999-1003.

Belton, RP et al. (2012). Hi-C: a comprehensive technique to capture the conformation of genomes. \emph{Methods} 58, 268-276.
}

\keyword{preprocessing}
