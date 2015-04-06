.loadIndices <- function(y, frag.chrs, restrict=NULL)
# A quick and dirty function for index loading with multiple libraries. This
# produces a list which describes the necessary HDF5 object corresponding to
# each chromosome combination for each library. 
{
	y <- path.expand(y)
	overall <- list()
	ni<-length(y)
	do.restrict <- length(restrict)!=0L
	check.chrs <- !missing(frag.chrs)
	if (!check.chrs) { warning("no protection against nonsensical chromosomes in the data") }

	for (ix in 1:ni) {
		all.data <- loadChromos(y[ix])
		current <- split(all.data$targets, all.data$anchors)

		for (ac in names(current)) {
			if (check.chrs && !ac %in% frag.chrs) { 
				warning("'", ac, "' in data is not present in fragment chromosomes") 
				next
			}
			if (do.restrict && !ac %in% restrict) { next }
			if (is.null(overall[[ac]])) { overall[[ac]]<-list() }
			subcurrent <- current[[ac]]
				
			for (tc in subcurrent) {
				if (check.chrs && !tc %in% frag.chrs) { 
					warning("'", tc, "' in data is not present in fragment chromosomes") 
					next
				}
				if (do.restrict && !tc %in% restrict) { next } 
				if (is.null(overall[[ac]][[tc]])) { overall[[ac]][[tc]] <- logical(ni) }
				overall[[ac]][[tc]][ix] <- TRUE
			}
		}
	}
	return(overall)
}

.getPairs <- function(y, anchor, target) { 
	y <- path.expand(y)
	h5read(y, file.path(anchor, target)) 
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

.writePairs <- function(pairs, y, anchor, target) {
	y <- path.expand(y)
	rownames(pairs) <- NULL
	if (h5write(pairs, y, file.path(anchor, target))) { stop("failed to add tag pair data to '%s'", y) }
	return(invisible(NULL))
}

loadChromos <- function(file) 
# A user-accessible function, to see what chromosomes are available in the
# file. This is designed to allow users to pull out one chromosome or another.
#
# written by Aaron Lun
# created 3 November 2014.
{
	current <- h5ls(file)
	keep <- current$otype=="H5I_DATASET"
	return(data.frame(anchors=basename(current$group[keep]),
		targets=current$name[keep]))
}

loadData <- function(file, anchor, target) 
# Friendly user-exposed handling of read pair extraction, when the user
# isn't sure of the order of the anchor/target chromosomes.
#
# written by Aaron Lun
# created 3 November 2014
{
	stopifnot(is.character(anchor) & is.character(target) & is.character(file))
	tryCatch(.getPairs(file, anchor, target), error=function(e) {
		out <- tryCatch(.getPairs(file, target, anchor), error=function(e) {
			stop("no dataset corresponding to this anchor/target combination")
		})
		warning("anchor and target definitions are reversed")
		out
	})
}	
