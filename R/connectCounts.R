connectCounts <- function(files, param, regions, filter=1L, type="any", second.regions=NULL)
# This counts the number of connections between specified regions in the genome (i.e. between regions
# in 'anchor' and regions in 'target'). This is designed to make it easier to analyze results in terms
# of genes. Note that everything is rounded up to the nearest outside restriction site (or to the
# nearest inside restriction site, depending).
#
# written by Aaron Lun
# a long time ago.
# last modified 27 March 2015
{
	nlibs <- length(files)
	if (nlibs==0L) { stop("number of libraries must be positive") } 
	filter<-as.integer(filter)
	fragments <- param$fragments

	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	# Checking out which regions overlap with each fragment.
	if (any(strand(regions)!="*")) { 
		warning("stranded region ranges have no interpretation, coercing unstrandedness") 
		strand(regions) <- "*"
	}
	if (any(strand(fragments)!="*")) { 
		warning("stranded fragment ranges have no interpretation, coercing unstrandedness") 
		strand(fragments) <- "*"
	}
	olaps <- suppressWarnings(findOverlaps(fragments, regions, type=type))
	frag.ids <- queryHits(olaps)
	reg.ids <- subjectHits(olaps)
	regions <- .redefineRegions(olaps, fragments, regions)

	# Including second region information.
	if (!is.null(second.regions)) { 
		if (is(second.regions, "GRanges")) {
			if (any(strand(second.regions)!="*")) { 
				warning("stranded region ranges have no interpretation, coercing unstrandedness") 
				strand(second.regions) <- "*"
			}
			lap2 <- suppressWarnings(findOverlaps(fragments, second.regions, type=type))
			to.add.query <- queryHits(lap2)
			to.add.subject <- subjectHits(lap2)
			second.regions <- .redefineRegions(lap2, fragments, second.regions)
		} else {
			second.regions <- as.integer(second.regions)
			if (second.regions < 0) { stop("bin size must be a positive integer") }
			binned <- .getBinID(fragments, second.regions)
			to.add.query <- 1:length(fragments)
			to.add.subject <- binned$id
			second.regions <- binned$region
		}

		n.first <- length(regions)
		n.second <- length(second.regions)
		regions <- suppressWarnings(c(regions, second.regions))
		regions$is.second <- rep(c(FALSE, TRUE), c(n.first, n.second))

		frag.ids <- c(frag.ids, to.add.query)
		reg.ids <- c(reg.ids, to.add.subject + n.first)
		o <- order(frag.ids, reg.ids)
		frag.ids <- frag.ids[o]
		reg.ids <- reg.ids[o]

		regions$original <- c(1:n.first, 1:n.second)
	} else {
		is.second <- NULL
		regions$original <- 1:length(regions)
	}

	# Figuring out which regions are anchor or targets.
	fdata <- .delimitFragments(fragments)
	matched <- match(as.character(seqnames(regions)), fdata$chr)
	if (any(is.na(matched))) { stop("chromosome present in regions and not in fragments") }

	nregs <- length(regions)
	o <- order(matched, start(regions), end(regions), 1:nregs) # Preserve order, if expanded intervals are identical.
	regions <- regions[o]

	ranked <- integer(nregs)
	ranked[o] <- 1:length(o)
	reg.ids <- ranked[reg.ids]
	by.frag <- .retrieveHits(frag.ids, length(fragments))

	# Setting up output containers.
    full.sizes <- integer(nlibs)
	out.counts <- list(matrix(0L, 0, nlibs))
	out.right <- out.left <- list(integer(0))
	idex<-1L

	chrs <- seqlevels(fragments)
	my.chrs <- unique(runValue(seqnames(regions)))
    overall <- .loadIndices(files, chrs, restrict)

	for (anchor in names(overall)) {
		current<-overall[[anchor]]
		for (target in names(current)) {
			if (!.checkIfPairOK(restrict, anchor, target)) { next }

           	pairs <- .baseHiCParser(current[[target]], files, anchor, target, discard=discard, cap=cap)
			for (lib in 1:length(pairs)) { 
				.checkIndexOK(fragments, anchor, pairs[[lib]]$anchor.id)
				if (anchor!=target) { .checkIndexOK(fragments, target, pairs[[lib]]$target.id) }
	            full.sizes[lib] <- full.sizes[lib] + nrow(pairs[[lib]])
			}
			if (! (target %in% my.chrs) || ! (anchor %in% my.chrs)) { next }	

			# Extracting counts. Running through the fragments and figuring out what matches where.
			out <- .Call(cxx_count_connect, pairs, by.frag$start, by.frag$end, reg.ids, filter, regions$is.second)
			if (is.character(out)) { stop(out) }
			out.counts[[idex]] <- out[[3]]
			out.left[[idex]] <- out[[1]]
			out.right[[idex]] <- out[[2]]
			idex <-  idex + 1L
		}
	}

	out.counts <- do.call(rbind, out.counts)
	anchors <- unlist(out.left)
	targets <- unlist(out.right)
	o.all <- order(anchors, targets)

	return(DIList(counts=out.counts[o.all,,drop=FALSE], totals=full.sizes, 
		anchors=anchors[o.all], targets=targets[o.all], regions=regions,
		exptData=List(param=param)))
}

.retrieveHits <- function(frag.id, nfrags) { 
	start <- end <- integer(nfrags)
	is.first <- c(TRUE, diff(frag.id)!=0L)
	start[frag.id[is.first]] <- which(is.first)
	end[frag.id[is.first]] <- c(which(is.first)[-1], length(frag.id)+1L)
	return(list(start=start, end=end))
}

.redefineRegions <- function(olaps, fragments, regions) {
	so <- subjectHits(olaps)
	qo <- queryHits(olaps)
	reo <- order(so, qo)
	so <- so[reo]
	qo <- qo[reo]
	
	s.rle <- rle(so)
	r.fin <- cumsum(s.rle$length)
	r.beg <- r.fin - s.rle$length + 1L
	ranges(regions)[s.rle$value] <- IRanges(start(fragments[qo[r.beg]]), end(fragments[qo[r.fin]]))
	# The preceding step is valid because fragments are sorted and non-nested.

	nfrags <- integer(length(regions))
	nfrags[s.rle$value] <- s.rle$length
	regions$nfrags <- nfrags
	return(regions)		
}
