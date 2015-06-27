\name{DNaseHiC}
\alias{prepPseudoPairs}
\alias{segmentGenome}

\title{Methods for processing DNase Hi-C data}
\description{Processing of BAM files for DNase Hi-C into index files}

\usage{
segmentGenome(bs, size=500) 
prepPseudoPairs(bam, param, file, dedup=TRUE, yield=1e7, 
    ichim=TRUE, chim.span=1000, minq=NA)
}

\arguments{
	\item{bs}{a \code{BSgenome} object, or a character string pointing to a FASTA file, or a named integer vector of chromosome lengths}
	\item{size}{an integer scalar indicating the size of the pseudo-fragments}
	\item{bam}{a character string containing the path to a name-sorted BAM file}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{file}{a character string specifying the path to an output index file}
	\item{dedup}{a logical scalar indicating whether marked duplicate reads should be removed}
	\item{yield}{a numeric scalar specifying the number of reads to extract at every iteration}
	\item{ichim}{a logical scalar indicating whether invalid chimeras should be counted}
	\item{chim.span}{an integer scalar specifying the maximum span between a chimeric 3' end and a mate read}
	\item{minq}{an integer scalar specifying the minimum mapping quality for each read}
}

\details{
DNase Hi-C involves random fragmentation with DNase instead of restriction enzymes.
This is accommodated in \code{diffHic} by partitioning the genome into small pseudo-fragments, using \code{segmentGenome}.
Reads are then assigned into these pseudo-fragments using \code{prepPseudoPairs}.
The rest of the analysis pipeline can then be used in the same manner as that for standard Hi-C.

The behaviour of \code{prepPseudoPairs} is almost identical to that for \code{\link{preparePairs}}, if the latter were asked to assign reads into pseudo-fragments.
However, for \code{prepPseudoPairs}, no reporting or removal of self-circles or dangling ends is performed, as these have no meaning for artificial fragments.
Also, invalidity of chimeras is determined by checking whether the 3' end is more than \code{chim.span} away from the mate read, rather than checking for localization in different fragments.

The size of the pseudo-fragments is determined by, well, \code{size} in \code{segmentGenome}.
Smaller sizes provide better resolution but increase computational work.
Needless to say, the \code{param$fragments} field should contain the output from \code{segmentGenome}, rather than from \code{\link{cutGenome}}.
Also see \code{\link{cutGenome}} documentation for a warning about the chromosome names.

Some loss of spatial resolution is inevitable when reads are summarized into pseudo-fragments.
This is largely irrelevant, though, as counting across the interaction space will ultimately use much larger bins (usually at least 2 kbp).
}

\value{
For \code{segmentGenome}, a \code{GRanges} object is produced containing the coordinates of the pseudo-fragments in the specified genome.

For \code{prepPseudoPairs}, a HDF5-formatted index file is produced at the specified location.
A list of diagnostic vectors are also returned in the same format as that from \code{\link{preparePairs}}, without the \code{same.id} entry.
}


\seealso{
\code{\link{preparePairs}},
\code{\link{cutGenome}}
}

\author{Aaron Lun}

\examples{
require(BSgenome.Ecoli.NCBI.20080805)
segmentGenome(BSgenome.Ecoli.NCBI.20080805)
segmentGenome(BSgenome.Ecoli.NCBI.20080805, size=1000)

# Pretend that this example is DNase Hi-C.
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
pseudo <- segmentGenome(seqlengths(cuts), size=50) 
param <- pairParam(pseudo) 

tmpf <- "gunk.h5"
prepPseudoPairs(hic.file, param, tmpf)
prepPseudoPairs(hic.file, param, tmpf, dedup=FALSE)
prepPseudoPairs(hic.file, param, tmpf, minq=50)
prepPseudoPairs(hic.file, param, tmpf, chim.span=20)

\dontshow{
unlink(tmpf)
}
}

\keyword{preprocessing}
