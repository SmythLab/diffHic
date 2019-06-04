\name{cutGenome}
\alias{cutGenome}

\title{Cut up the genome}

\description{Perform an in silico restriction digest of a target genome.}

\usage{
cutGenome(bs, pattern, overhang=4L)
}

\arguments{
\item{bs}{A \linkS4class{BSgenome} object or a character string containing a path to a FASTA file.}
\item{pattern}{A character vector containing one or more recognition sites.}
\item{overhang}{An integer vector specifying the length of the 5' overhang for each element in \code{pattern}.}
}

\value{
A \linkS4class{GRanges} object containing the boundaries of each restriction fragment in the genome.
}

\details{
This function simulates a restriction digestion of a specified genome, 
given the recognition site and 5' overhang of the restriction enzyme.
The total sequence spanned by each fragment is recorded including the two sticky ends after they are filled in.
No support is currently provided for searching the reverse strand, so the recognition site should be an inverse palindrome.

The genome should be specified as a \linkS4class{BSgenome} object.
However, a character string can also be provided, specifying a FASTA file containing all the reference sequences in a genome.
The latter may be necessary to synchronise the fragments with the genome used for alignment.

Multiple restriction enzymes can be specified by passing vectors to \code{pattern} and \code{overhang}.
All recognition sites are expected to be inverse palindromes,
and for any given enzyme, the width of \code{pattern} and \code{overhang} must be both odd or even.

No attempt is made to remove fragments that cannot be physically formed,
e.g., from recognition sites that overlap with themselves or each other.
This generally is not problematic for downstream analysis, as they are short and will not have many assigned reads.
If they are a concern, most of them can be removed by simply applying a suitable threshold (e.g., 10 bp) on the fragment width.
However, the best solution is to simply choose (combinations of) restriction enzymes that do not overlap.
}

\section{Note on FASTA sequence names}{
If \code{bs} is a FASTQ file, the chromosome names in the FASTQ headers will be loaded faithfully by \code{cutGenome}.
However, many mapping pipelines will drop the rest of the name past the first whitespace when constructing the alignment index.
To be safe, users should ensure that the chromosome names in the FASTQ headers consist of one word.
Otherwise, there will be a discrepancy between the chromosome names in the output \code{GRanges}
and those in the BAM files after alignment.
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

# Trying with FASTA files.
x <- system.file("extdata", "fastaEx.fa", package="Biostrings")
cutGenome(x, "AGCT", overhang=2)
cutGenome(x, "AGCT", overhang=4)

# Multiple sites with different overhangs are supported.
cutGenome(x, c("AGCT", "AGNCT"), overhang=c(4, 3))
}

\seealso{
\code{\link{matchPattern}}
}

\keyword{preprocessing}
