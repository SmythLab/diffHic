# Constructor object.
DIList <- function(counts, totals=colSums(counts), anchors, targets, regions, exptData=List(), ...) {
	.Defunct(new="InteractionSet", package="InteractionSet")
}

DI2IS <- function(x) {
    InteractionSet(list(counts=x@counts), 
                   GInteractions(anchor1=x@anchors, anchor2=x@targets, regions=x@regions, mode="reverse"),
                   colData=x@colData, metadata=x@exptData)
}
