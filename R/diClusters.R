diClusters <- function(data.list, result.list, target, equiweight=TRUE, cluster.args=list(), pval.col="PValue", fc.col=NA, grid.length=21, iterations=4)
# Performs post-hoc clustering of significant bin pairs, to control the 
# cluster-level FDR at 'target'. Adjusts the weights so that the contribution
# from each set of bin pairs is the same.
#
# written by Aaron Lun
# created 12 January 2016
# last modified 13 January 2016
{
    # Setting initial parameters.    
    if (missing(target)) {
        target <- 0.05
        warning("unspecified 'target' for the cluster-level FDR set to 0.05")
    }
    expanded.args <- as.list(match.call(clusterPairs, do.call(call, c("clusterPairs", cluster.args))))
    if (is.null(expanded.args$tol)) {
        cluster.args$tol <- 1
        warning("'tol' for 'clusterPairs' set to a default of 1 bp")
    }

    # Checking inputs.
    if (.notList(data.list)) {
        data.list <- list(data.list) 
    }
    if (.notList(result.list) || is.data.frame(result.list)) { 
        result.list <- list(result.list) 
    }
	nset <- length(data.list)
	if (nset!=length(result.list)) { stop("data list must have same length as result list") }
	for (x in seq_len(nset)) {
		if (length(data.list[[x]])!=nrow(result.list[[x]])) {
 		   	stop("corresponding entries of data and result lists must have same number of entries") 
        }
	}

    # Computing the adjusted FDR for all of these samples (with equi-weighting).
    nentries <- lengths(data.list)
    if (equiweight) {
        weights <- rep(1/nentries, nentries)
    } else {
        weights <- rep(1, sum(nentries))
    }

    all.ps <- in.each.group <- vector("list", nset)
    last <- 0L
    for (x in seq_len(nset)) { 
        all.ps[[x]] <- result.list[[x]][,pval.col]
        in.each.group[[x]] <- last + seq_len(nrow(result.list[[x]]))
        last <- last + nrow(result.list[[x]])
    }
    all.ps <- unlist(all.ps)
    adjp <- csaw:::.weightedFDR(all.ps, weights)

    # Getting the sign.
    if (is.na(fc.col)) { 
        all.signs <- lapply(nentries, logical)
    } else {
        all.signs <- lapply(result.list, function(x) { x[,fc.col] > 0 })
    }

	# Controlling the cluster-level FDR.
    FUN <- function(sig, index.only=TRUE) {
        pos.data.list <- neg.data.list <- vector("list", nset)
        for (x in seq_len(nset)) { 
            cur.sig <- sig[in.each.group[[x]]]
            pos.data.list[[x]] <- data.list[[x]][cur.sig & all.signs[[x]],]
            neg.data.list[[x]] <- data.list[[x]][cur.sig & !all.signs[[x]],]
        }
        pos.clust.all <- do.call(clusterPairs, c(pos.data.list, cluster.args, index.only=index.only))
        neg.clust.all <- do.call(clusterPairs, c(neg.data.list, cluster.args, index.only=index.only))

        # Extracting the indices.
        if (index.only) { 
            pos.clust <- pos.clust.all
            neg.clust <- neg.clust.all
            additional <- max(sapply(pos.clust, function(x) { 
                keep <- !is.na(x)
                if (!any(keep)) { return(0) } 
                return(max(x[keep])) 
            }))
        } else {
            pos.clust <- pos.clust.all$indices
            neg.clust <- neg.clust.all$indices
            additional <- length(pos.clust.all$interactions)
        }
    
        # Assembling it back into a single return value.
        clust.indices <- vector("list", nset)
        for (x in seq_len(nset)) {
            cur.sig <- sig[in.each.group[[x]]]
            cur.signs <- all.signs[[x]][cur.sig] 
            full.ids <- integer(sum(cur.sig))
            full.ids[cur.signs] <- pos.clust[[x]]
            full.ids[!cur.signs] <- neg.clust[[x]] + additional
            clust.indices[[x]] <- full.ids
        }
        names(clust.indices) <- names(data.list)
        if (index.only) { return(clust.indices) }
        list(indices=clust.indices, interactions=c(pos.clust.all$interactions, neg.clust.all$interactions))
    }
    out <- controlClusterFDR(target=target, adjp=adjp, FUN=function(sig) { 
        unlist(FUN(sig)) 
    }, weights=weights, grid.length=grid.length, iterations=iterations)
    sig <- adjp <= out$threshold
    clusters <- FUN(sig, index.only=FALSE)

    # Cleaning up the output.
    for (x in seq_len(nset)) {
        full.ids <- rep(NA_integer_, nrow(data.list[[x]]))
        full.ids[sig[in.each.group[[x]]]] <- clusters$indices[[x]]
        clusters$indices[[x]] <- full.ids
    }
    clusters$FDR <- out$FDR    
    return(clusters)
}

.notList <- function(x) { (!is.list(x) && !is(x, "List")) }
