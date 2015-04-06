\name{pairParam}
\docType{class}

\alias{pairParam}
\alias{pairParam-class}
\alias{$,pairParam-method}
\alias{show,pairParam-method}
\alias{reform}
\alias{reform,pairParam-method}

\title{pairParam class and methods}
\description{Class to specify read pair loading parameters}

\details{
Each pairParam object stores a number of parameters to extract reads from a BAM file.
Slots are defined as:
\describe{
	\item{\code{fragments}:}{a \code{GRanges} object containing the coordinates of the restriction fragments}
	\item{\code{restrict}:}{a character vector or 1-by-2 matrix containing the names of allowable chromosomes from which reads will be extracted}
	\item{\code{discard}:}{a \code{GRanges} object containing intervals in which any alignments will be discarded}
	\item{\code{cap}:}{an integer scalar, specifying the maximum number of read pairs per pair of restriction fragments}
}

The \code{fragments} object defines the genomic interval spanned by each restriction fragment.
All reads are generated around restriction sites, so the spatial resolution of the experiment depends on such sites.
The object can be can be obtained by applying \code{\link{cutGenome}} on an appropriate \code{BSgenome} object.

If \code{restrict} is supplied, reads will only be extracted for the specified chromosomes. 
This is useful to restrict the analysis to interesting chromosomes, e.g., no contigs/scaffolds or mitochondria. 
\code{restrict} can also be a 1-by-2 matrix, specifying a single pair of chromosomes over which read pairs are to be counted.

If \code{discard} is set, a read will be removed if the corresponding alignment is wholly contained within the supplied ranges. 
Any pairs involving reads discarded in this manner will be ignored.
This is useful for removing unreliable alignments in repeat regions.

If \code{cap} is set to a non-\code{NA} value, an upper bound will be placed on the number of read pairs that are counted for each fragment pair (after any removal due to \code{discard}).
This protects against spikes in the read pair density throughout the interaction space.
Such spikes may be caused by technical artifacts like PCR duplication or repeats, which were not successfully removed in prior processing steps.

% Nominally, cap should be sample-specific, as you should use different caps for different library sizes.
% However, the expected count per fragment pair is quite low at your typical sequencing depth. If you define
% the cap as a proportional function of the expected count, you'd end up with a value less than 1 for most
% library sizes. This would just get rounded to a constant value of 1 anyway.
}

\section{Constructor}{
\describe{
	\item{}{\code{pairParam(fragments, discard=GRanges(), restrict=NULL, cap=NA)}:
		Creates a pairParam object. 
		Each argument is placed in the corresponding slot, with coercion into the appropriate type.
	}
}
}

\section{Subsetting}{
In the code snippets below, \code{x} is a pairParam object.
\describe{
	\item{}{\code{x$name}: 
	Returns the value in slot \code{name}.
	}
}
}

\section{Other methods}{
In the code snippets below, \code{x} is a pairParam object.
\describe{
    \item{}{\code{show(x)}:
	Describes the parameter settings in plain English.
  	}
	\item{}{\code{reform(x, ...)}:
		Creates a new pairParam object, based on the existing \code{x}.
		Any named arguments in \code{...} are used to modify the values of the slots in the new object, with type coercion as necessary.
	}
}	
}

\author{
Aaron Lun
}

\seealso{
\code{\link{cutGenome}},
\code{\link{squareCounts}}
}

\examples{
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
		
blah <- pairParam(cuts)
blah <- pairParam(cuts, discard=GRanges("chrA", IRanges(1, 10)))
blah <- pairParam(cuts, restrict='chr2')
blah$fragments
blah$restrict
blah$cap

# Use 'reform' if only some arguments need to be changed.
blah
reform(blah, restrict='chr3')
reform(blah, discard=GRanges())
reform(blah, cap=10)

# Different restrict options.
pairParam(cuts, restrict=c('chr2', 'chr3'))
pairParam(cuts, restrict=cbind('chr2', 'chr3'))
}

\keyword{counting}
