# Copied from csaw package

.weightedFDR <- function(p, w) {
    if (length(p)!=length(w)) { stop("weights and p-value vector are not of same length") }
    o <- order(p)
    p <- p[o]
    w <- w[o]
    adjp <- numeric(length(o))
    adjp[o] <- rev(cummin(rev(sum(w)*p/cumsum(w))))
    pmin(adjp, 1)
}
