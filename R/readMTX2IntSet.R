readMTX2IntSet <- function(mtx, bed, verbose=TRUE)
# Read contact matrix in Matrix Market Exchange Format into an Interaction Set
# Gordon Smyth
# Created 22 June 2018. Last modified 22 July 2018.
{
#   Read genomic regions from BED file
#   Read first line to get number of columns
    FirstLine <- utils::read.table(bed,nrows=1L,comment.char="",quote="",stringsAsFactors=FALSE)
    nbedcols <- ncol(FirstLine)
    if(nbedcols < 3L) stop("BED file must have at least 3 columns")
    colClasses <- c("character","integer","integer")
    if(nbedcols > 3L) colClasses <- c(colClasses,rep.int("NULL",nbedcols-3L))
    Regions <- utils::read.table(bed,colClasses=colClasses,comment.char="",quote="")
    if(verbose) cat("Read",bed,"with",nrow(Regions),"regions\n")
    GR <- GRanges(Regions[,1],IRanges(Regions[,2]+1L,Regions[,3]))

#   Check for Matrix Market header (so it can be skipped)
	con <- file(mtx[1], "r")	
	skip <- 0L
#   Count number of lines starting with comment: header are these plus one line.
	repeat {
		txt <- readLines(con,n=1L)
		if(substring(txt,1,1) != "%") break
		skip <- skip + 1L
	}
	if(skip > 0L) {
	    if(verbose) cat("First MTX file header:",txt,"\n")
	    skip <- skip + 1L
	}
	close(con)

#   Read the mtx files
    nfiles <- length(mtx)
    if(nfiles==1L) {

#       If only one mtx file, make InteractionSet directly. No need to merge interactions.
        x <- utils::read.table(mtx,skip=skip,sep="",colClasses=c("integer","integer","integer"),comment.char="",quote="")
        if(verbose) cat("Read",mtx,"\n")
        lib.size <- sum(x[,3])
        GI <- GInteractions(x[,1],x[,2],GR)
        IS <- InteractionSet(as.matrix(x[,3]),GI,colData=DataFrame(lib.size=lib.size))

    } else {

#       Read mtx files into a list of data.frames
        x <- list()
        lib.size <- numeric(nfiles)
        hash <- list()
        for (i in seq_len(nfiles)) {
            x[[i]] <- utils::read.table(mtx[i],skip=skip,sep="",colClasses=c("integer","integer","integer"),comment.char="",quote="")
            if(verbose) cat("Read",mtx[i],"\n")
            hash[[i]] <- x[[i]][,1] + x[[i]][,2] / (length(GR)+1)
            x[[i]] <- x[[i]][,3] # saving memory
            lib.size[i] <- sum(x[[i]])
        }

#       Find union of interactions
        if(verbose) cat("Merging ...\n")
        hashu <- hash[[1]]
        for (i in 2:nfiles) hashu <- unique(c(hashu,hash[[i]]))
        Anchor1 <- floor(hashu)
        Anchor2 <- (hashu - Anchor1) * (length(GR) + 1)
        Anchor1 <- as.integer(Anchor1)
        Anchor2 <- as.integer(round(Anchor2))
        GI <- GInteractions(Anchor1,Anchor2,GR)

#       Merge counts into one matrix
        counts <- matrix(0L,length(hashu),nfiles)
        for (i in seq_len(nfiles)) {
            m <- match(hash[[i]],hashu)
            counts[m,i] <- x[[i]]
        }
        IS <- InteractionSet(counts,GI,colData=DataFrame(lib.size=lib.size))
    }

#   Set colnames and assayNames
    assayNames(IS) <- "counts"
    mtx <- sub("\\.gz$","",mtx)
    colnames(IS) <- limma::removeExt(mtx)

    IS
}
