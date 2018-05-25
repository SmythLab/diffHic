filterPeaks <- function(data, enrichment, assay.bp=1, assay.neighbors=NULL, get.enrich=FALSE,
                        min.enrich=log2(1.5), min.count=5, min.diag=2L, ...)
# This is a wrapper function that takes the enrichment values from filterPeaks and 
# identifies those bin pairs that satisfy the enrichment threshold, and other cut-offs.
# 
# written by Aaron Lun
# created 23 March 2015
{
    .check_StrictGI(data)
    ab <- scaledAverage(data, assay.id=assay.bp, ...)

    # Computing the enrichment value, if it is missing.
    if (missing(enrichment)) { 
        assay.neighbors <- .neighbor_locales(assay.neighbors)
        n.names <- .neighbor_numbers(assay.neighbors)

        ab.q <- scaledAverage(data, assay.id=assay.neighbors[1], scale=mcols(data)[,n.names[1]], ...) # Quadrant
        ab.v <- scaledAverage(data, assay.id=assay.neighbors[2], scale=mcols(data)[,n.names[2]], ...) # Vertical
        ab.h <- scaledAverage(data, assay.id=assay.neighbors[3], scale=mcols(data)[,n.names[3]], ...) # Horizontal 
        ab.s <- scaledAverage(data, assay.id=assay.neighbors[4], scale=mcols(data)[,n.names[4]], ...) # Surrounding

        enrichment <- ab - pmax(ab.q, ab.v, ab.h, ab.s)
    }
    if (get.enrich) return(enrichment)
   
    # Otherwise filtering based on the enrichment values. 
    keep <- enrichment > min.enrich 
    if (!is.null(min.count)) { 
        keep <- keep & ab > aveLogCPM(min.count, lib.size=mean(data$totals), ...)
    } 
    if (!is.null(min.diag)) {
        cur.dist <- pairdist(data, type="diag")  
        keep <- keep & (is.na(cur.dist) | cur.dist >= min.diag)
    }
    return(keep)
}

