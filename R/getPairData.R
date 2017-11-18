getPairData <- function(file, param)
# This retrieves the fragment sizes, relative orientations and gaps from each directory produced by 
# preparePairs. This is a convenience function which allows people to avoid loading the entire directory 
# in (or manually parsing the said directory).
#
# written by Aaron Lun
# created 20 September 2014
# last modified 18 November 2017
{
    loadfuns <- preloader(file, param=param, 
        retain=c("anchor1.id", "anchor1.pos", "anchor1.len",
                 "anchor2.id", "anchor2.pos", "anchor2.len"))
    alll <- allo <- alli <- vector("list", sum(lengths(loadfuns)))
    ix <- 1L

    for (ax in names(loadfuns)) {
        current <- loadfuns[[ax]] 
        for (tx in names(current)) { 
            extracted <- current[[tx]][[1]]()
            yielded <- .getStats(extracted, ax==tx, param$fragments)
            alll[[ix]] <- yielded$length
            allo[[ix]] <- yielded$orientation
            alli[[ix]] <- yielded$insert
            ix <- ix + 1L
        }
    }

    # Return objects.
    return(data.frame(length=unlist(alll), 
        orientation=unlist(allo), insert=unlist(alli)))
}

.getStats <- function(incoming, same.chr, fragments) 
# Gets the statistics for everything, including the fragment length,
# strand orientation and gap size.
{
    if (!nrow(incoming)) {
        return(list(length=integer(0), orientation=integer(0), gap=integer(0)))
    }
    output <- .Call(cxx_pair_stats, incoming$anchor1.id, incoming$anchor2.id, incoming$anchor1.pos, incoming$anchor2.pos,
        incoming$anchor1.len, incoming$anchor2.len, same.chr, start(fragments), end(fragments))
    names(output) <- c("length", "orientation", "insert")
    return(output)    
}

