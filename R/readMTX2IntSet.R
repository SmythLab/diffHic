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
    mcols(GR) <- NULL

    nsamples <- length(mtx)
    collected.first <- collected.second <- collected.counts <- vector("list", nsamples)
    for (mdx in seq_len(nsamples)) {
        current <- read.table(mtx[mdx], comment.char="%", quote="", 
            colClasses=c("integer", "integer", if (as.integer) "integer" else "numeric"))

        if (current[1,1]!=length(GR) || current[1,2]!=length(GR)) {
            stop("dimensions in 'mtx' are not equal to number of regions in 'bed'")
        }

        collected.first[[mdx]] <- current[-1,1]
        collected.second[[mdx]] <- current[-1,2]
        collected.counts[[mdx]] <- current[-1,3] 
    }

    # Defining the common set of interactions.
    sample.ids <- rep(seq_len(nsamples), lengths(collected.first))
    if (nsamples) {
        collected.first <- unlist(collected.first)
        collected.second <- unlist(collected.second)
        collected.counts <- unlist(collected.counts)
    } else {
        collected.first <- collected.second <- collected.counts <- integer(0)
    }
    all.gi <- GInteractions(collected.first, collected.second, GR, mode="reverse")
    all.gi$samples <- sample.ids
    all.gi$counts <- collected.counts

    all.gi <- sort(all.gi)
    first.of.type <- !duplicated(all.gi)
    row.num <- cumsum(first.of.type)

    # Creating the output matrix.
    output <- matrix(if (as.integer) 0L else 0, nrow=sum(first.of.type), ncol=nsamples)
    for (mdx in seq_len(nsamples)) {
        chosen <- all.gi$samples==mdx
        output[row.num[chosen],mdx] <- all.gi$counts[chosen]
    }

    # Creating the output ISet.
    uniq.gi <- all.gi[first.of.type]
    mcols(uniq.gi) <- NULL
    output <- InteractionSet(output, uniq.gi, colData=DataFrame(totals=colSums(output)))
    assayNames(output) <- "counts"
    colnames(output) <- names(mtx)
    metadata(output)$width <- median(width(regions(output)))
    return(output)
}
