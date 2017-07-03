\name{savePairs}
\alias{savePairs}

\title{Save Hi-C read pairs}

\description{Save a dataframe of read pairs into a directory structure for rapid chromosomal access.}

\usage{
savePairs(x, file, param)
}

\arguments{
	\item{x}{A dataframe with integer fields \code{anchor1.id} and \code{anchor2.id}.
Each row corresponds to a single read pair.}
	\item{file}{A character string specifying the path for the output index file.}
	\item{param}{A \code{pairParam} object containing read extraction parameters.
In particular, \code{param$fragments} should contain genomic regions corresponding to the \code{anchor*.id} values.}
}

\value{
An index file is produced at the specified \code{file} location, containing the interaction data.
A \code{NULL} value is invisibly returned.
}

\details{
This function facilitates the input of processed Hi-C data from other sources into the current pipeline. 
Each row of \code{x} corresponds to a read pair, and each entry in \code{x$anchor1.id} and \code{x$anchor2.id} contains an index for \code{param$fragments}.
Thus, the pair of indices for each row denotes the the interacting regions for each read pair.
These regions are generally expected to be restriction fragments in conventional Hi-C experiments.

Obviously, the coordinates of the restriction fragment boundaries in \code{param$fragments} should correspond to the genome to which the reads were aligned.
These can be generated using the \code{\link{cutGenome}} function from any given \code{BSgenome} object or the FASTA files used for alignment.
Values of \code{param$discard} and \code{param$restrict} will not be used here and can be ignored.

Any additional fields in \code{x} will also be saved to file.
Users are recommended to put in \code{anchor1.pos}, \code{anchor1.len}, \code{anchor2.pos} and \code{anchor2.len} fields.
These should mimic the output of \code{\link{preparePairs}}:
\describe{
\item{\code{anchorY.pos}:}{Integer field, containing the 1-based genomic position of the left-most aligned base of read Y.}
\item{\code{anchorY.len}:}{Integer field, containing the length of the alignment of read Y on the reference sequence. 
This should be multiplied by -1 if the alignment was on the negative strand.}
}
These fields enable the use of more \pkg{diffHic} functions, e.g., removal of reads in \code{param$discard} during counting with \code{\link{squareCounts}},
correct calculation of statistics with \code{\link{getPairData}}, quality control with \code{\link{prunePairs}}.

For storing DNase Hi-C data, \code{param$fragments} should be empty but the \code{seqinfo} should contain the lengths and names of all chromosomes.
Here, the input \code{anchor1.id} and \code{anchor2.id} should contain indices of the \code{seqlengths}.
This specifies the chromosome to which each read is aligned, e.g., an \code{anchor1.id} of 2 means that read 1 is aligned to the second chromosome in \code{seqlengths}.
Note that, for this type of data, it is essential to store the position and length fields mentioned above.

When constructing the output file, \code{x} will be resorted by \code{anchor1.id}, then \code{anchor2.id}.
If necessary, anchor1 and anchor2 indices will be switched such that the former is never less than the latter. 
For DNase Hi-C data, both of these fields will ultimately be set to zero - see \code{\link{prepPseudoPairs}} for more details.
}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <-readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts)

n <- 1000
all.a <- as.integer(runif(n, 1L, length(cuts)))
all.t <- as.integer(runif(n, 1L, length(cuts)))
x <- data.frame(anchor1.id=all.a, anchor2.id=all.t,
	anchor1.pos=runif(1:100), anchor1.len=10, 
	anchor2.pos=runif(1:100), anchor2.len=-10)

fout <- "temp2.h5"
savePairs(x, fout, param)
require(rhdf5)
head(h5read(fout, "chrA/chrA"))

\dontshow{
unlink(fout)
}
}

\author{Aaron Lun}

\seealso{
	\code{\link{preparePairs}}, 
	\code{\link{cutGenome}}
}

\keyword{preprocessing}
