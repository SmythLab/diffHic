\name{cutGenome}
\alias{cutGenome}

\title{Cut up the genome}

\description{Perform an in silico restriction digest of a target genome.}

\usage{
cutGenome(bs, pattern, overhang=4L)
}

\arguments{
	\item{bs}{a \code{BSgenome} object or a character string pointing to a FASTA file}
	\item{pattern}{character string describing the recognition site}
	\item{overhang}{integer scalar specifying the length of the 5' overhang}
}

\value{
A \code{GRanges} object containing the boundaries of each restriction fragment in the genome.
}

\details{
This function simulates a restriction digestion of a specified genome, given the recognition site and 5' overhang of the cutter
The total sequence spanned by each fragment is recorded, including the two sticky ends. 
No support is currently provided for searching the reverse strand, so the recognition site should be an inverse palindrome.

The genome should be specified as a \code{BSgenome} object.
However, a character string can also be provided, specifying a FASTA file containing all the reference sequences in a genome.
The latter may be necessary to synchronise the fragments with the genome used for alignment.

Note that some of the reported fragments may be physically impossible to form, e.g., for overlapping sites or consecutive sites when \code{overhang==nchar(pattern)}.
Nonetheless, they are still reported to maintain the correspondence between fragments and cut sites.
Cleavage sites on the forward strand can be obtained as the \code{start} locations of all fragments (excepting the first fragment on each chromosome).
}

\section{Warning}{
If \code{bs} is a FASTQ file, the chromosome names in the FASTQ headers will be loaded faithfully by \code{cutGenome}.
However, many mapping pipelines will drop the rest of the name past the first whitespace when constructing the alignment index.
To be safe, users should ensure that the chromosome names in the FASTQ headers consist of one word.
Otherwise, there will be a discrepancy between the chromosome names in the output \code{GRanges}, and those in the BAM files after alignment.
}

% Interpretations of consecutive sites is generally tricky.
% For starters, the 'remainder' is so low that the strands are unlikely to stay stuck together until the fill-in step.
% This becomes an impossibility if remainder is zero, such that ssDNA is formed after cleavage of consecutive sites.
% It's also unclear whether cleavage is even possible when the recognition site is at the very end of the fragment (e.g., after one cleavage).
% That's not even considering the grief that's possible when a site overlaps with itself.
% In short, the fragments that will be reported by cutGenome might be a bit silly in such cases; but, in many respects, it doesn't matter.
% Fragment-level resolution is never used, and those formed between consecutive sites will be so small that they'll have no effect on read assignment.

\author{Aaron Lun}

\examples{
require(BSgenome.Ecoli.NCBI.20080805)

cutGenome(Ecoli, "AAGCTT", overhang=4L) # HindIII
cutGenome(Ecoli, "CCGCGG", overhang=2L) # SacII
cutGenome(Ecoli, "AGCT", overhang=0L) # AluI

# Trying with FastA files.
x <- system.file("extdata", "fastaEx.fa", package="Biostrings")
cutGenome(x, "AGCT", overhang=2)
cutGenome(x, "AGCT", overhang=4)
}

\seealso{
\code{\link{matchPattern}}
}

\keyword{preprocessing}
