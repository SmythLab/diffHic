connectCounts <- function(files, param, regions, filter=1L, type="any", second.regions=NULL)
# This counts the number of connections between specified regions in the genome (i.e. between regions
# in 'anchor' and regions in 'target'). This is designed to make it easier to analyze results in terms
# of genes. Note that everything is rounded up to the nearest outside restriction site (or to the
# nearest inside restriction site, depending on the overlap parameters).
#
# written by Aaron Lun
# a long time ago.
# last modified 18 November 2017
{
    if (.isDNaseC(param)) {
        return(.connectCountsRaw(files, param, regions, filter=filter, type=type, second.regions=second.regions))
    }    
    
    nlibs <- length(files)
    if (nlibs==0L) { 
        stop("number of libraries must be positive") 
    }
    filter <- as.integer(filter)

    # Processing regions.
    fragments <- param$fragments
    reg.out <- .processRegions(param, regions, type, second.regions)
    regions <- reg.out$regions
    frag.ids <- reg.out$frag.ids
    reg.ids <- reg.out$reg.ids

    # Ordering regions, consistent with the previous definitions of anchor/targets.
    # Stable sort preserves order, if expanded intervals are identical (NA's get sorted towards the end and can be ignored).
    ordered.chrs <- as.character(runValue(seqnames(fragments)))
    matched <- match(as.character(seqnames(regions)), ordered.chrs)
    o <- order(matched, start(regions), end(regions)) 

    nfrags <- length(fragments)
    regions <- regions[o]
    ranked <- integer(length(regions))
    ranked[o] <- seq_along(o)
    reg.ids <- ranked[reg.ids]

    use.second <- !is.null(second.regions)
    if (!use.second) {
        by.frag1 <- .retrieveHits(frag.ids, nfrags)
        by.frag2 <- by.frag1
        reg.id1 <- reg.ids
        reg.id2 <- NULL
    } else {
        is2 <- regions$is.second[reg.ids]
        reg.id1 <- reg.ids[!is2]
        by.frag1 <- .retrieveHits(frag.ids[!is2], nfrags) 
        reg.id2 <- reg.ids[is2]
        by.frag2 <- .retrieveHits(frag.ids[is2], nfrags) 
    }

    # Setting up output containers.
    full.sizes <- integer(nlibs)
    out.counts <- list(matrix(0L, 0, nlibs))
    out.right <- out.left <- list(integer(0))
    idex <- 1L

    my.chrs <- unique(runValue(seqnames(regions)))
    loadfuns <- preloader(files, param=param, retain=c("anchor1.id", "anchor2.id"))
    for (anchor in names(loadfuns)) {
        current <- loadfuns[[anchor]]
        for (target in names(current)) {
            curfuns <- current[[target]]

            pairs <- vector("list", nlibs)
            for (lib in seq_len(nlibs)) { 
                pairs[[lib]] <- curfuns[[lib]]()
            }
            full.sizes <- full.sizes + sapply(pairs, FUN=nrow)

            # This check needs to be done after loading to obtain a valid 'full.sizes'.
            if (! (target %in% my.chrs) || ! (anchor %in% my.chrs)) { 
                next 
            }
            
            # Extracting counts. Running through the fragments and figuring out what matches where.
            out <- .Call(cxx_count_connect, pairs, by.frag1$start, by.frag1$end, reg.id1, 
                         by.frag2$start, by.frag2$end, reg.id2, filter)
            out.counts[[idex]] <- out[[3]]
            out.left[[idex]] <- out[[1]]
            out.right[[idex]] <- out[[2]]
            idex <-  idex + 1L
        }
    }

    .generateOutput(out.counts, out.left, out.right, full.sizes, regions, param)
}

##############################################################################################

.processRegions <- function(param, regions, type, second.regions)
# Processes the regions into a common format for further use. Namely, we do
# pre-screening to remove entries with chromosomes beyond those in 'fragments',
# and identify the restriction fragments overlapping each region..
{
    # Eliminating irrelevant strand information and metadata.
    fragments <- param$fragments
    strand(fragments) <- "*"
    strand(regions) <- "*"
    mcols(regions) <- NULL

    # Checking out which regions overlap with each fragment.
    olaps <- suppressWarnings(findOverlaps(fragments, regions, type=type))
    frag.ids <- queryHits(olaps)
    reg.ids <- subjectHits(olaps)
    regions <- .redefineRegions(olaps, fragments, regions)

    if (!is.null(second.regions)) { 
        if (is(second.regions, "GRanges")) {
            strand(second.regions) <- "*"
            mcols(second.regions) <- NULL

            lap2 <- suppressWarnings(findOverlaps(fragments, second.regions, type=type))
            to.add.query <- queryHits(lap2)
            to.add.subject <- subjectHits(lap2)
            second.regions <- .redefineRegions(lap2, fragments, second.regions)
            second.original <- seq_along(second.regions)

        } else {
            second.regions <- as.integer(second.regions)
            if (second.regions < 0) { stop("bin size must be a positive integer") }
            binned <- .assignBins(param, second.regions, restricted=FALSE)
            to.add.query <- seq_along(fragments)
            to.add.subject <- binned$id 
            second.regions <- binned$region
            second.original <- rep(NA_integer_, length(second.regions))
        }

        n.first <- length(regions) # need here, as 'regions' is modified in the next line!
        regions <- suppressWarnings(c(regions, second.regions))
        regions$is.second <- rep(c(FALSE, TRUE), c(n.first, length(second.regions)))
        regions$original <- c(seq_len(n.first), second.original)

        frag.ids <- c(frag.ids, to.add.query)
        reg.ids <- c(reg.ids, to.add.subject + n.first)
        o <- order(frag.ids, reg.ids)
        frag.ids <- frag.ids[o]
        reg.ids <- reg.ids[o]
     } else {
        regions$original <- seq_along(regions)
    }

    return(list(regions=regions, frag.ids=frag.ids, reg.ids=reg.ids))
}

.retrieveHits <- function(frag.id, nfrags) 
# Figures out the start and end index of the queryHits vector for
# each fragment, allowing it to rapidly index the subjectHits vector.
{ 
    start <- end <- integer(nfrags)
    is.first <- c(TRUE, diff(frag.id)!=0L)
    start[frag.id[is.first]] <- which(is.first)
    end[frag.id[is.first]] <- c(which(is.first)[-1], length(frag.id)+1L)
    return(list(start=start, end=end))
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
    ranges(regions)[s.rle$value] <- IRanges(start(fragments)[qo[r.beg]], 
                                            end(fragments)[qo[r.fin]])
    # The preceding step is valid because fragments are sorted and non-nested.

    nfrags <- integer(length(regions))
    nfrags[s.rle$value] <- s.rle$length
    regions$nfrags <- nfrags
    return(regions)        
}

.generateOutput <- function(out.counts, out.left, out.right, full.sizes, regions, param) 
# A function to generate output.
{
    out.counts <- do.call(rbind, out.counts)
    anchors <- unlist(out.left)
    targets <- unlist(out.right)
    out <- InteractionSet(list(counts=out.counts), colData=DataFrame(totals=full.sizes), 
        interactions=GInteractions(anchor1=anchors, anchor2=targets, regions=regions, mode="reverse"), 
        metadata=List(param=param))
    sort(out)
}

##############################################################################################

.connectCountsRaw <- function(files, param, regions, filter=1L, type="any", second.regions=NULL)
# An equivalent function for DNase-C data. This uses the linkOverlaps() machinery
# to do the heavy lifting, as the C++ code written above is designed for restriction fragments.
#
# written by Aaron Lun
# created 17 March 2017
# last modified 18 November 2017
{
    nlibs <- length(files)
    if (nlibs==0L) { stop("number of libraries must be positive") } 
    filter <- as.integer(filter)
    
    # Processing regions.
    strand(regions) <- "*"
    mcols(regions) <- NULL
    if (!is.null(second.regions)) { 
        if (is.numeric(second.regions)) { 
            region2 <- .createBins(param, second.regions)$region
        } else {
            region2 <- second.regions
            mcols(region2) <- NULL
            strand(region2) <- "*"
            region2$nfrags <- 0L
        }
        region1 <- regions
        region1$nfrags <- 0L
        n.first <- length(region1)
        n.second <- length(region2)

        # Redefining the full set of regions.
        regions <- suppressWarnings(c(region1, region2))
        regions$is.second <- rep(c(FALSE, TRUE), c(n.first, n.second))
        d1 <- seq_len(n.first)
        d2 <- seq_len(n.second) 
        regions$original <- c(d1, d2)
        d2 <- d2 + n.first
    } else {
        regions$original <- seq_along(regions)
    }

    # Setting up output vectors.
    full.sizes <- integer(nlibs)
    out.counts <- list(matrix(0L, 0, nlibs))
    out.right <- out.left <- list(integer(0))
    idex <- 1L

    my.chrs <- unique(runValue(seqnames(regions)))
    loadfuns <- preloader(files, param=param, retain=c("anchor1.pos", "anchor1.len", "anchor2.pos", "anchor2.len"))
    for (anchor in names(loadfuns)) {
        current <- loadfuns[[anchor]]
        for (target in names(current)) {
            curfuns <- current[[target]]

            # Extracting counts.
            pairs <- vector("list", nlibs)
            for (lib in seq_len(nlibs)) { 
                pairs[[lib]] <- curfuns[[lib]]()
            }
            full.sizes <- full.sizes + sapply(pairs, FUN=nrow)

            # Again, this needs to be done after full.sizes is collated.
            if (! (target %in% my.chrs) || ! (anchor %in% my.chrs)) { 
                next 
            }

            # Forming a GInteractions object.  
            collected <- vector("list", length(pairs))
            for (lib in seq_along(pairs)) {
                cpair <- pairs[[lib]] 
                sgi <- suppressWarnings(GInteractions(GRanges(anchor, IRanges(cpair$anchor1.pos, width=abs(cpair$anchor1.len))),
                                                      GRanges(target, IRanges(cpair$anchor2.pos, width=abs(cpair$anchor2.len))),
                                                      mode="strict"))
                if (is.null(second.regions)) { 
                    li <- linkOverlaps(sgi, regions, type=type)
                } else {
                    li <- linkOverlaps(sgi, region1, region2, type=type)
                    li$subject1 <- d1[li$subject1]
                    li$subject2 <- d2[li$subject2]
                }
                li <- li[,c("subject1", "subject2")]
                o <- do.call(order, li)
                collected[[lib]] <- li[o,]
            }
        
            # Compiling the counts together.
            out <- .Call(cxx_count_reconnect, collected, filter)
            out.counts[[idex]] <- out[[3]]
            out.left[[idex]] <- out[[1]]
            out.right[[idex]] <- out[[2]]
            idex <-  idex + 1L
        }
    }

    .generateOutput(out.counts, out.left, out.right, full.sizes, regions, param)
}
