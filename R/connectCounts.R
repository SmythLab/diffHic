connectCounts <- function(files, param, regions, filter=1L, type="any", second.regions=NULL, restrict.regions=FALSE)
# This counts the number of connections between specified regions in the genome (i.e. between regions
# in 'anchor' and regions in 'target'). This is designed to make it easier to analyze results in terms
# of genes. Note that everything is rounded up to the nearest outside restriction site (or to the
# nearest inside restriction site, depending on the overlap parameters).
#
# written by Aaron Lun
# a long time ago.
{
    if (!.isDNaseC(param)) { 
        frag.out <- .processRegionsFrag(param, regions=regions, type=type, second.regions=second.regions, restricted=restrict.regions)
        out <- .generalLoopFUN(files, param, frag.out$regions, retain=c("anchor1.id", "anchor2.id"), filter=filter, linkFUN=.linkFrag, 
                               fragments=param$fragments, olap1=frag.out$overlap1, olap2=frag.out$overlap2)
    } else {
        raw.out <- .processRegionsRaw(param, regions=regions, second.regions=second.regions, restricted=restrict.regions)
        out <- .generalLoopFUN(files, param, raw.out$regions, retain=c("anchor1.pos", "anchor1.len", "anchor2.pos", "anchor2.len"),
                               filter=filter, linkFUN=.linkRaw, regions=raw.out$regions, region1=raw.out$region1, region2=raw.out$region2, 
                               index1=raw.out$index1, index2=raw.out$index2, type=type)
    }

    # Returning an ISet object.
    obj <- InteractionSet(list(counts=do.call(rbind, out$counts)), 
                          colData=DataFrame(totals=out$totals),
                          interactions=GInteractions(anchor1=unlist(out$left),
                                                     anchor2=unlist(out$right), 
                                                     regions=out$regions, mode="reverse"),
                          metadata=List(param=param))
    sort(obj)
}

##############################################################################################

.processRegionsFrag <- function(param, regions, type, second.regions, restricted=FALSE)
# Processes the regions into a common format for further use. Namely, we 
# identify the restriction fragments overlapping each region, and expand
# the ranges to include the overlapped fragments. 
{
    # Eliminating irrelevant strand information and metadata.
    fragments <- param$fragments
    strand(fragments) <- "*"
    strand(regions) <- "*"
    mcols(regions) <- NULL

    # Performing restriction.
    rest.out <- .restrictRegions(regions, param, restricted)
    regions <- rest.out$regions
    first.original <- rest.out$original

    # Checking out which regions overlap with each fragment.
    olap <- suppressWarnings(findOverlaps(fragments, regions, type=type))
    regions <- .redefineRegions(olap, fragments, regions)

    if (!is.null(second.regions)) {
        if (is(second.regions, "GRanges")) {
            strand(second.regions) <- "*"
            mcols(second.regions) <- NULL

            # Again, subsetting the regions now, if requested.
            sec.rest.out <- .restrictRegions(second.regions, param, restricted)
            second.regions <- sec.rest.out$regions
            second.original <- sec.rest.out$original

            olap2 <- suppressWarnings(findOverlaps(fragments, second.regions, type=type))
            second.regions <- .redefineRegions(olap2, fragments, second.regions)

        } else {
            # Forming bins if we have a scalar 'second.regions'.
            second.regions <- as.integer(second.regions)
            if (second.regions < 0) {
                stop("bin size must be a positive integer")
            }
            binned <- .assignBins(param, second.regions, restricted=restricted)
            second.regions <- binned$region
            second.original <- rep(NA_integer_, length(second.regions))

            keep <- !is.na(binned$id) # As NA's are possible when restricted=TRUE.
            olap2 <- Hits(which(keep), binned$id[keep], nLnode=length(fragments), 
                          nRnode=length(second.regions), sort.by.query=TRUE)
            olap2 <- sort(olap2)
        }

        n.first <- length(regions) # needs to be done here, as we modify it in the next line!
        regions <- suppressWarnings(c(regions, second.regions))
        regions$is.second <- rep(c(FALSE, TRUE), c(n.first, length(second.original)))
        regions$original <- c(first.original, second.original)

        # Remapping Hits due to the concatenation of second.regions onto regions.
        olap2 <- remapHits(olap2, new.nRnode=length(regions), Rnodes.remapping=seq_along(second.regions) + n.first)
     } else {
        regions$original <- first.original
        olap2 <- NULL
    }

    return(list(regions=regions, overlap1=olap, overlap2=olap2))
}

.redefineRegions <- function(olaps, fragments, regions)
# Stretches out each region to encompass the fragments it overlaps
# (regions that don't overlap any fragments are not modified,
# which allows safe use for DNase-C data.
{
    so <- subjectHits(olaps)
    qo <- queryHits(olaps)
    reo <- order(so, qo)
    so <- so[reo]
    qo <- qo[reo]

    s.rle <- rle(so)
    r.fin <- cumsum(s.rle$length)
    r.beg <- r.fin - s.rle$length + 1L
    
    # This step is valid because fragments are sorted and non-nested.
    ranges(regions)[s.rle$value] <- IRanges(start(fragments)[qo[r.beg]],
                                            end(fragments)[qo[r.fin]])

    nfrags <- integer(length(regions))
    nfrags[s.rle$value] <- s.rle$length
    regions$nfrags <- nfrags
    return(regions)
}

.linkFrag <- function(anchor, target, curpairs, fragments, olap1, olap2) {
    gi <- GInteractions(curpairs$anchor1.id, curpairs$anchor2.id, fragments) 
    if (is.null(olap2)){ 
        linked <- linkOverlaps(gi, olap1) 
    } else {
        linked <- linkOverlaps(gi, olap1, olap2) 
    }
    linked$query <- NULL
    return(linked)
}

# Note that anchor and target don't get used by .linkFrag, but are nontheless necessary for 
# compatibility with the general linkFUN template.

##############################################################################################

.processRegionsRaw <- function(param, regions, second.regions, restricted=FALSE) {
    strand(regions) <- "*"
    mcols(regions) <- NULL
    regions$nfrags <- 0L

    # Applying restriction if requested.
    rest.out <- .restrictRegions(regions, param, restricted)
    regions <- rest.out$regions
    first.original <- rest.out$original

    if (!is.null(second.regions)) {
        if (is.numeric(second.regions)) {
            # Creating bins if requested.
            region2 <- .createBins(param, second.regions, restricted=restricted)$region
            original2 <- rep(NA_integer_, length(region2))

        } else {
            # Cleaning up the second set of ranges.
            region2 <- second.regions
            mcols(region2) <- NULL
            strand(region2) <- "*"
            region2$nfrags <- 0L

            # Applying restriction if necessary.
            res.out2 <- .restrictRegions(region2, param, restricted=restricted)
            region2 <- res.out2$region
            original2 <- res.out2$original
        }

        region1 <- regions
        n.first <- length(region1)
        n.second <- length(region2)

        # Redefining the full set of regions.
        regions <- suppressWarnings(c(region1, region2))
        regions$is.second <- rep(c(FALSE, TRUE), c(n.first, n.second))
        regions$original <- c(first.original, original2)

        # Creating indices for entry in the concatenated object.
        d1 <- seq_len(n.first)
        d2 <- seq_len(n.second) + n.first
    } else {
        regions$original <- first.original
        region1 <- region2 <- NULL
        d1 <- d2 <- NULL
    }

    return(list(regions=regions, region1=region1, region2=region2, index1=d1, index2=d2))
}

.linkRaw <- function(anchor, target, curpairs, regions, region1, region2, index1, index2, type) { 
    sgi <- suppressWarnings(GInteractions(GRanges(anchor, IRanges(curpairs$anchor1.pos, width=abs(curpairs$anchor1.len))),
        GRanges(target, IRanges(curpairs$anchor2.pos, width=abs(curpairs$anchor2.len))), mode="strict"))
     if (is.null(region2)) {
        li <- linkOverlaps(sgi, regions, type=type)
    } else {
        li <- linkOverlaps(sgi, region1, region2, type=type)
        li$subject1 <- index1[li$subject1]
        li$subject2 <- index2[li$subject2]
    }
    li$query <- NULL
    return(li)
}

##############################################################################################

.restrictRegions <- function(regions, param, restricted) {
    # Subsetting the regions now, if requested.
    if (restricted && length(param$restrict)) {
        original <- which(seqnames(regions) %in% param$restrict)
        regions <- regions[original]
    } else {
        original <- seq_along(regions)
    }
    return(list(regions=regions, original=original))
}

.generalLoopFUN <- function(files, param, regions, retain, filter, linkFUN, ...) 
# This is a generic looping function that can be applied to both 
# restriction Hi-C and DNase-C data, to avoid rewriting it between
# the two modes. The ellipsis arguments are passed to linkFUN.
{
    nlibs <- length(files)
    if (nlibs==0L) {
        stop("number of libraries must be positive")
    }
    filter <- as.integer(filter)

    # Setting up output containers.
    full.sizes <- integer(nlibs)
    out.counts <- list(matrix(0L, 0, nlibs))
    out.right <- out.left <- list(integer(0))
    idex <- 1L

    my.chrs <- unique(runValue(seqnames(regions)))
    loadfuns <- preloader(files, param=param, retain=retain)
    for (anchor in names(loadfuns)) {
        current <- loadfuns[[anchor]]
        for (target in names(current)) {
            curfuns <- current[[target]]

            pairs <- vector("list", nlibs)
            for (lib in seq_len(nlibs)) {
                pairs[[lib]] <- curfuns[[lib]]()
            }
            full.sizes <- full.sizes + vapply(pairs, FUN=nrow, FUN.VALUE=0L)

            # This check needs to be done after loading to obtain a valid 'full.sizes'.
            if (! (target %in% my.chrs) || ! (anchor %in% my.chrs)) {
                next
            }
               
            # Calling the link function.
            collected <- vector("list", length(pairs))
            for (lib in seq_along(pairs)) {
                linked <- linkFUN(anchor=anchor, target=target, curpairs=pairs[[lib]], ...)
                o <- order(linked$subject1, linked$subject2)
                collected[[lib]] <- as.list(linked[o,])
            }

            # Compiling the counts together.
            out <- .Call(cxx_count_connect, collected, filter)
            out.counts[[idex]] <- out[[3]]
            out.left[[idex]] <- out[[1]]
            out.right[[idex]] <- out[[2]]
            idex <-  idex + 1L
        }
    }

    return(list(counts=out.counts, left=out.left, right=out.right, regions=regions, totals=full.sizes))
}

