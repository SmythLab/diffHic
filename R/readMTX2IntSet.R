#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom utils read.table
#' @importFrom InteractionSet InteractionSet GInteractions
#' @importFrom SummarizedExperiment assayNames
#' @importFrom S4Vectors DataFrame mcols<- metadata<-
#' @importFrom BiocGenerics match unique sort colnames width
#' @importMethodsFrom InteractionSet match c unique sort
#' @importFrom stats median
readMTX2IntSet <- function(mtx, bed, as.integer=TRUE)
# Read contact matrix in Matrix Market Exchange Format into an Interaction Set
#
# written by Gordon Smyth
# with modifications by Aaron Lun
# created 22 June 2018
{
    GR <- import.bed(bed)
    GR <- sort(GR)
    mcols(GR) <- NULL

    collected.gi <- collected.counts <- vector("list", length(mtx))
    for (mdx in seq_along(mtx)) {
        current <- read.table(mtx[mdx], comment.char="%", quote="", 
            colClasses=c("integer", "integer", if (as.integer) "integer" else "numeric"))

        if (current[1,1]!=length(GR) || current[1,2]!=length(GR)) {
            stop("dimensions in 'mtx' are not equal to number of regions in 'bed'")
        }

        collected.gi[[mdx]] <- GInteractions(current[-1,1], current[-1,2], GR, mode="reverse") 
        collected.counts[[mdx]] <- current[-1,3] 
    }

    # Defining the common set of interactions.
    if (length(mtx)) { 
        all.gi <- do.call(c, collected.gi)
        all.gi <- unique(sort(all.gi))
    } else {
        all.gi <- GInteractions(integer(0), integer(0), GR)
    }

    output <- matrix(if (as.integer) 0L else 0, nrow=length(all.gi), ncol=length(collected.gi))
    for (mdx in seq_along(mtx)) {
        location <- match(collected.gi[[mdx]], all.gi)
        output[location,mdx] <- collected.counts[[mdx]]
    }

    output <- InteractionSet(output, all.gi, colData=DataFrame(totals=colSums(output)))
    assayNames(output) <- "counts"
    colnames(output) <- names(mtx)
    metadata(output)$width <- median(width(regions(output)))
    return(output)
}
