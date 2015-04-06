\name{diffHicUsersGuide}
\alias{diffHicUsersGuide}
\alias{diffHic}

\title{View diffHic user's guide}
\description{Finds the location of the user's guide and opens it for viewing.}

\usage{
diffHicUsersGuide(view=TRUE)
}

\arguments{
  \item{view}{logical scalar specifying whether the document should be opened}
}

\value{
A character string giving the file location. 
If \code{view=TRUE}, the system's default PDF document reader is started and the user's guide is opened.
}

\details{
The diffHic package is designed for the detection of differential interactions from Hi-C data.
It provides methods for read pair counting, normalization, filtering and statistical analysis via edgeR.
As the name suggests, the diffHic user's guide for can be obtained by running this function.

For non-Windows operating systems, the PDF viewer is taken from \code{Sys.getenv("R_PDFVIEWER")}. 
This can be changed to \code{x} by using \code{Sys.putenv(R_PDFVIEWER=x)}. 
For Windows, the default viewer will be selected to open the file. 

Note that this guide is not a true vignette as it is not generated using \code{\link{Sweave}} when the package is built. 
This is due to the time-consuming nature of the code when run on realistic case studies.
}

\seealso{
\code{\link{system}}
}

\author{Aaron Lun}

\examples{
# To get the location:
diffHicUsersGuide(view=FALSE)
# To open in pdf viewer:
\dontrun{diffHicUsersGuide()}
}

\keyword{documentation}
