mergeCMs <- function(..., deflate.args=list(), filter=1, binned=TRUE) 
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
    to.keep <- total.sum >= filter

    # Constructing standardized ISet objects.
    isets <- vector("list", length(inputs))
    totals <- numeric(length(inputs))
    for (i in seq_along(inputs)) {
        cm <- inputs[[i]]
        isets[[i]] <- do.call(deflate, c(list(cm, extract=to.keep), deflate.args))
        totals[i] <- sum(as.matrix(cm))
    }

    # Merging together into a single output object.
    data <- do.call(cbind, isets)
    interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")
    colnames(data) <- names(inputs)
    data$totals <- totals

    if (binned) {
        metadata(data)$width <- median(width(regions(data)))
    }
    return(data)
}
