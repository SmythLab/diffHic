.loadIndices <- function(y, frag.chrs=NULL, restrict=NULL)
# A quick and dirty function for index loading with multiple libraries. This
# produces a list which describes the necessary HDF5 object corresponding to
# each chromosome combination for each library. 
{
	y <- path.expand(y)
	overall <- list()
	ni<-length(y)
	
	# Checking how we want to do restriction.
	do.restrict <- length(restrict)!=0L
	paired <- attributes(restrict)$paired
	if (is.null(paired)) { 
        paired <- FALSE 
    } else if (paired) { 
        nr <- length(restrict)/2L
    }

	for (ix in seq_len(ni)) {
		all.data <- loadChromos(y[ix])

        # Initial restriction on both entries being present in the list.
        if (do.restrict) {
            keep <- all.data$anchor1 %in% restrict & all.data$anchor2 %in% restrict
            all.data <- all.data[keep,]
        }
		current <- split(all.data$anchor2, all.data$anchor1)

        # Throwing an error if we see a chromosome not present in 'frag.chrs'.
        if (!is.null(frag.chrs)) { 
            leftovers <- setdiff(c(all.data$anchor1, all.data$anchor2), frag.chrs)
            if (length(leftovers)) {
                stop("'", leftovers[1], "' in '", y[ix], "' is not present in fragment chromosomes")
            }
        }

		for (ac in names(current)) {
			subcurrent <- current[[ac]]

			# Paired restriction if necessary.
			if (do.restrict && paired) {
                tmatch <- which(restrict==ac)
                tmatch <- tmatch + nr * (-1L)^(tmatch > nr) # Getting the chromosome in the other column.
                subcurrent <- intersect(subcurrent, restrict[tmatch])
                if (!length(subcurrent)) { next }
            }
			
            if (is.null(overall[[ac]])) { overall[[ac]]<-list() }
			for (tc in subcurrent) {
				if (is.null(overall[[ac]][[tc]])) { overall[[ac]][[tc]] <- logical(ni) }
				overall[[ac]][[tc]][ix] <- TRUE
			}
		}
	}
	return(overall)
}

.getPairs <- function(y, anchor1, anchor2) { 
	y <- path.expand(y)
	out <- h5read(y, file.path(anchor1, anchor2)) 

    # For legacy purposes:
    colnames(out) <- sub("anchor\\.", "anchor1.", colnames(out))
    colnames(out) <- sub("target\\.", "anchor2.", colnames(out))
    return(out)
}

.initializeH5 <- function(y) {
	y <- path.expand(y)
	if (file.exists(y)) { unlink(y, recursive=TRUE) } 
	if (!h5createFile(y)) { stop(sprintf("failed to create '%s'", y)) }
	return(invisible(NULL))
}

.addGroup <- function(y, anchor) {
	y <- path.expand(y)
	if (!h5createGroup(y, anchor)) { stop("failed to add '%s' group to '%s'", anchor, y) }
	return(invisible(NULL))
}

.writePairs <- function(pairs, y, anchor1, anchor2) {
	y <- path.expand(y)
	rownames(pairs) <- NULL
	if (h5write(pairs, y, file.path(anchor1, anchor2))) { stop("failed to add tag pair data to '%s'", y) }
	return(invisible(NULL))
}

loadChromos <- function(file) 
# A user-accessible function, to see what chromosomes are available in the
# file. This is designed to allow users to pull out one chromosome or another.
#
# written by Aaron Lun
# created 3 November 2014
# last modified 6 January 2016
{
	current <- h5ls(file)
	keep <- current$otype=="H5I_DATASET"
	return(data.frame(anchor1=basename(current$group[keep]), 
            anchor2=current$name[keep], stringsAsFactors=FALSE))
}

loadData <- function(file, anchor1, anchor2) 
# Friendly user-exposed handling of read pair extraction, when the user
# isn't sure of the order of the anchor chromosomes.
#
# written by Aaron Lun
# created 3 November 2014
# last modified 22 November 2015
{
	stopifnot(is.character(anchor1))
	stopifnot(is.character(anchor2))
	stopifnot(is.character(file))
	tryCatch(.getPairs(file, anchor1, anchor2), error=function(e) {
		out <- tryCatch(.getPairs(file, anchor2, anchor1), error=function(e) {
			stop("no dataset corresponding to this anchor combination")
		})
		warning("anchor definitions are reversed")
		out
	})
}	
