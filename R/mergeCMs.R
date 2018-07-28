#' @export
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom InteractionSet deflate regions anchors
#' @importMethodsFrom InteractionSet as.matrix cbind
#' @importFrom methods as
#' @importClassesFrom InteractionSet ReverseStrictGInteractions
#' @importFrom BiocGenerics colnames width
#' @importFrom S4Vectors metadata<-
#' @importFrom stats median
mergeCMs <- function(..., deflate.args=list())
# This function merges contact matrices into a single InteractionSet object.
# It effectively mimics a squareCounts call, but with contact matrix inputs.
#
# written by Aaron Lun
# created 21 September 2017
{
    inputs <- list(...)
    if (length(inputs)==0) {
        stop("at least one ContactMatrix must be supplied")
    }

    # Choosing which interactions to keep.
    total.sum <- as.matrix(inputs[[1]])
    ref.regions <- regions(inputs[[1]])
    ref.anchors <- anchors(inputs[[1]], id=TRUE)

    for (cm in inputs[-1]) {
        if (length(ref.regions)!=length(regions(cm)) || !all(regions(cm)==ref.regions)) {
            stop("'regions' should be the same across ContactMatrix objects")
        } else if (!identical(ref.anchors, anchors(cm, id=TRUE))) {
            stop("anchor indices should be the same across ContactMatrix objects")
        }
        total.sum <- total.sum + as.matrix(cm)
    }

    to.keep <- total.sum >= 1L 

    # Constructing standardized ISet objects for cbinding.
    isets <- mapply(deflate, inputs, MoreArgs=c(list(extract=to.keep), deflate.args), SIMPLIFY=FALSE)
    data <- do.call(cbind, isets)
    data$totals <- unname(colSums(assay(data, withDimnames=FALSE)))

    # Cleaning up the output.
    interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")
    assayNames(data) <- "counts"
    colnames(data) <- names(inputs)
    metadata(data)$width <- median(width(regions(data)))
    return(data)
}
