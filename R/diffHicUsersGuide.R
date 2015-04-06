diffHicUsersGuide <- function(view=TRUE)
# Find, and optionally show, the diffHic user's guide 
# 
# written by Aaron Lun, based on equivalent from edgeR by Gordon Smyth
# 7 November 2014.
{
	ugloc <- system.file("doc", "diffHicUsersGuide.pdf", package="diffHic")
	if (view) {
		if(.Platform$OS.type == "windows") {
			shell.exec(ugloc)
		} else {
			system(paste(Sys.getenv("R_PDFVIEWER"),ugloc,"&"))
		}
	}
	return(ugloc)
}
