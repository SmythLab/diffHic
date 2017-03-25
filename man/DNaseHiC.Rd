\name{DNaseHiC}
\alias{prepPseudoPairs}
\alias{segmentGenome}
\alias{emptyGenome}

\title{Methods for processing DNase Hi-C data}
\description{Processing of BAM files for DNase Hi-C into index files}

\usage{
emptyGenome(bs)

# Deprecated
segmentGenome(bs) 

# Deprecated
prepPseudoPairs(bam, param, file, dedup=TRUE, minq=NA, ichim=TRUE, 
    chim.span=1000, output.dir=NULL, storage=5000L)
}

\arguments{
	\item{bs}{a \code{BSgenome} object, or a character string pointing to a FASTA file, or a named integer vector of chromosome lengths}
	\item{bam}{a character string containing the path to a name-sorted BAM file}
	\item{param}{a \code{pairParam} object containing read extraction parameters}
	\item{file}{a character string specifying the path to an output index file}
	\item{dedup}{a logical scalar indicating whether marked duplicate reads should be removed}
	\item{minq}{an integer scalar specifying the minimum mapping quality for each read}
	\item{ichim}{a logical scalar indicating whether invalid chimeras should be counted}
	\item{chim.span}{an integer scalar specifying the maximum span between a chimeric 3' end and a mate read}
    \item{output.dir}{a character string specifying a directory for temporary files}
    \item{storage}{an integer scalar specifying the maximum number of pairs to store in memory before writing to file}
}

\details{
DNase Hi-C uses DNase to randomly fragment the genome, rather than using restriction fragments.
This requires some care to handle in \pkg{diffHic}, as most functions rely on fragment assignments in many functions.
To specify that the data are from a DNase Hi-C experiment, an empty GRanges object should be supplied as the \code{fragments} in \code{\link{pairParam}}.
Most functions will automatically recognise that the data are DNase Hi-C and behave appropriately.
This reflects the fact that no restriction fragments are involved in this analysis. 
Genome information will instead be extracted from the \code{seqlengths} of the GRanges object.

\code{prepPseudoPairs} and \code{segmentGenome} are deprecated in favour of \code{\link{preparePairs}} and \code{emptyGenome}, respectively.
In the case of \code{\link{preparePairs}}, it will automatically detect if the BAM file is to be treated as DNase Hi-C based on \code{param$fragments}.
}

\value{
For \code{emptyGenome}, an empty GRanges object is produced containing the sequence names and lengths.
}

\seealso{
\code{\link{preparePairs}},
\code{\link{cutGenome}}
}

\author{Aaron Lun}

\examples{
require(BSgenome.Ecoli.NCBI.20080805)
emptyGenome(Ecoli)
emptyGenome(c(chrA=100, chrB=200))
}

\keyword{preprocessing}
