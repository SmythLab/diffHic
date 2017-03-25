\name{preparePairs}
\alias{preparePairs}

\title{Prepare Hi-C pairs}

\description{Identifies the interacting pair of restriction fragments corresponding to each read pair in a Hi-C library.}

\usage{
preparePairs(bam, param, file, dedup=TRUE, minq=NA, ichim=TRUE, 
    chim.dist=NA, output.dir=NULL, storage=5000L)
}

\arguments{
	\item{bam}{a character string containing the path to a name-sorted BAM file}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{file}{a character string specifying the path to an output index file}
	\item{dedup}{a logical scalar indicating whether marked duplicate reads should be removed}
	\item{minq}{an integer scalar specifying the minimum mapping quality for each read}
	\item{ichim}{a logical scalar indicating whether invalid chimeras should be counted}
	\item{chim.dist}{an integer scalar specifying the maximum distance between segments for a valid chimeric read pair}
	\item{output.dir}{a character string specifying a directory for temporary files}
    \item{storage}{an integer scalar specifying the maximum number of pairs to store in memory before writing to file}
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
The anchor1 fragment is that with the higher genomic coordinate, i.e., the larger index in \code{param$fragments}.
The anchor2 fragment is that with the smaller coordinate/index. 
This definition avoids the need to consider both permutations of indices in a pair during downstream processing.
}

\section{Details of read pair processing}{
For pairs with chimeric reads, the alignment of the 5' end of the read is used to identify the interacting fragments.  
Invalid chimeras arise when the position of the 3' segment of a chimeric read is not consistent with that of the mate read. 
These are generally indicative of mapping errors but can also form due to non-specific ligation events. 
While they can be explicitly removed, setting \code{ichim=TRUE} is recommended to avoid excessive filtering of reads when alignment of short chimeric segments is inaccurate.

Invalid chimeras can be defined in two ways:
\itemize{
\item If \code{chim.dist=NA}, a chimeric pair is considered to be invalid if the 3' segment and the mate do not map onto the same restriction fragment in an inward-facing orientation.
This reflects the resolution limits of the Hi-C protocol.
\item If \code{chim.dist} is not \code{NA}, chimeras are defined based on distance.
A pair is considered invalid if the distance between the segment and mate is greater than \code{chim.dist}, or if the alignments are not inward-facing.
}
The second approach is more relevant in situations involving inefficient cleavage, where the mapping locations are broadly consistent but do not fall in the same restriction fragment.
The maximum size of the ligation products can be used as a reasonable value for \code{chim.dist}, e.g., 1000 bp. 

Self-circles are outward-facing read pairs mapped to the same restriction fragment.
These are formed from inefficient cross-linking and are generally uninformative.
Dangling ends are inward-facing read pairs mapped to the same fragment, and are generated from incomplete ligation of blunt ends.
Both constructs are detected and discarded within the function.
Note that this does not consider dangling ends or self-circles formed from incompletely digested fragments, which must be removed with \code{\link{prunePairs}}.

% Also note that overlapping pairs of self-circles are retained, as are dangling ends where one end extends past the other.
% This because they are neither valid self-circles or dangling ends.

If the 5' end of either read is unavailable (e.g. unmapped, mapping quality score below \code{minq}, marked as a duplicate), the read pair is recorded as unpaired and discarded.
No MAPQ filtering is performed when \code{minq} is set to \code{NA}.
Any duplicate read must be marked in the bit field of the BAM file using a tool like Picard's \code{MarkDuplicates} if it is to be removed with \code{dedup=TRUE}. 
For chimeric reads, the recommended approach is to designate the 5' end as the only primary or non-supplementary alignment. 
This ensures that the duplicate computations are performed on the most relevant alignments for each read pair.
}

\section{Processing DNase Hi-C experiments}{ 
DNase Hi-C involves random fragmentation with DNase instead of restriction enzymes.
To indicate that the data are generated by DNase Hi-C, an empty GRanges object should be supplied as the \code{fragments} in \code{\link{pairParam}}.
Genome information will instead be extracted from the \code{seqlengths} of the GRanges object.
This empty object can be generated by \code{\link{emptyGenome}}.

The BAM file can be processed with \code{preparePairs} in a manner that is almost identical to that of normal Hi-C experiments.
However, there are some key differences:
\itemize{
\item No reporting or removal of self-circles or dangling ends is performed, as these have no meaning when restriction fragments are not involved.
\item Invalidity of chimeras is determined by checking whether the 3' end is more than \code{chim.dist} away from the mate read.
If \code{chim.dist=NA}, it will be set to a default of 1000 bp.
\item All fragment IDs in the output HDF5 file will be set to zero.
}
}

\section{Miscellaneous information}{
If \code{output.dir} is not specified, a directory name is constructed using the \code{\link{tempfile}} command.
If it is specified, users should make sure that no file already exists with that name.
Otherwise, an error will be raised.
This directory is used to store intermediate files that will be eventually processed into the HDF5 output file.

For low-memory systems or in cases where there are many chromosome pairs, users may need to reduce the value of \code{storage}.
This will write data to file more frequently, which reduces memory usage at the cost of speed.

Users should note that the use of a \code{pairParam} object for input is strictly for convenience.
Only the value of \code{param$fragments} will be used.
Any non-empty values of \code{param$discard} and \code{param$restrict} will be ignored here.
Reads will not be discarded if they lie outside the specified chromosomes, or if they lie within blacklisted regions.
}

\value{
Multiple dataframe objects are stored within the specified \code{file} using the HDF5 format.
Each object corresponds to a pair of chromosomes, designated as the anchor1 (later) or anchor2 (earlier) chromosome based on the order of their names.
Each row of the dataframe contains information for a read pair, with one read mapped to each chromosome. 
The dataframe contains several integer fields:
\describe{
  \item{\code{anchor1.id}, \code{anchor2.id}:}{index of the anchor1 or anchor2 restriction fragment}
  \item{\code{anchor1.pos}, \code{anchor2.pos}:}{mapping coordinate of the read (or for chimeras, the 5' segment thereof) on the anchor1 or anchor2 fragment}
  \item{\code{anchor1.len}, \code{anchor2.len}:}{length of the alignment on the anchor1 or anchor2 fragment, set to a negative value if the alignment is on the reverse strand}
}
See \code{\link{getPairData}} for more details on \code{length}, \code{orientation} and \code{gap}.  

A list is also returned from the function, containing various diagnostics:
\describe{
	\item{\code{pairs}:}{an integer vector containing \code{total}, the total number of read pairs; \code{marked}, read pairs with at least one marked 5' end; \code{filtered}, read pairs where the MAPQ score for either 5' end is below \code{minq}; \code{mapped}, read pairs considered as successfully mapped (i.e., not filtered, and also not marked if \code{dedup=TRUE})}
	\item{\code{same.id}:}{an integer vector containing \code{dangling}, the number of read pairs that are dangling ends; and \code{self.circles}, the number of read pairs forming self-circles}
	\item{\code{singles}:}{an integer scalar specifying the number of reads without a mate}
	\item{\code{chimeras}:}{an integer vector containing \code{total}, the total number of read pairs with one chimeric read; \code{mapped}, chimeric read pairs with both 5' ends mapped; \code{multi}, mapped chimeric pairs with at least one successfully mapped 3' segment; and \code{invalid}, read pairs where the 3' location of one read disagrees with the 5' location of the mate}
}

For DNase Hi-C data, the \code{anchor1.id} and \code{anchor2.id} fields are set to zero, and the \code{same.id} field in the output list is removed.
}

\author{Aaron Lun,
with use of Rhtslib based on comments from Alessandro Mammana.
}

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

# Pretending it's DNase Hi-C:
eparam <- pairParam(cuts[0]) 
preparePairs(hic.file, eparam, tmpf)
preparePairs(hic.file, eparam, tmpf, dedup=FALSE)
preparePairs(hic.file, eparam, tmpf, minq=50)
preparePairs(hic.file, eparam, tmpf, chim.dist=20)

\dontshow{
unlink(tmpf)
}
}

\references{
Imakaev M et al. (2012). Iterative correction of Hi-C data reveals hallmarks of chromosome organization. \emph{Nat. Methods} 9, 999-1003.

Belton, RP et al. (2012). Hi-C: a comprehensive technique to capture the conformation of genomes. \emph{Methods} 58, 268-276.
}

\keyword{preprocessing}
