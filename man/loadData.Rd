\name{loadData}
\alias{loadData}
\alias{loadChromos}

\title{Load data from an index file}
\description{Load read pair data and chromosome names from a HDF5 index file.}

\usage{
loadChromos(file)
loadData(file, anchor, target)
}

\arguments{
\item{file}{a character string containing a path to a index file}
\item{anchor}{a character string, specifying the name of one chromosome in a pair}
\item{target}{a character string, specifying the name of the other chromosome in the pair}
}

\value{
The \code{loadChromos} function will return a dataframe with character fields \code{anchors} and \code{targets}.
Each row represents a pair of chromosomes, the names of which are stored in the fields. 
The presence of a row indicates that the data for the corresponding pair exists in the \code{file}.

The \code{loadData} function will return a dataframe where each row contains information for one read pair.
Refer to \code{\link{preparePairs}} for more details on the type of fields that are included.
}

\details{
The purpose of these function is to allow users to perform custom analyses by extracting the data manually from each index file.
This may be desirable, e.g., when preparing data for input into other tools.
To extract all data, users are advised to run \code{loadData} iteratively on each pair of chromosomes as obtained with \code{loadChromos}.

Note that \code{loadData} will successfully operate even if the anchor/target specification is mixed up.
In this case, it will return a warning to inform the user that the names should be switched.
}

\author{Aaron Lun}

\examples{
hic.file <- system.file("exdata", "hic_sort.bam", package="diffHic")
cuts <- readRDS(system.file("exdata", "cuts.rds", package="diffHic"))
param <- pairParam(cuts) 

tmpf <- "gunk.h5"
preparePairs(hic.file, param, tmpf)

loadChromos(tmpf)
loadData(tmpf, "chrA", "chrA")
loadData(tmpf, "chrB", "chrA")
loadData(tmpf, "chrA", "chrB")
try(loadData(tmpf, "chrA2", "chrB2"))

\dontshow{
unlink(tmpf)
}
}

\keyword{preprocessing}
