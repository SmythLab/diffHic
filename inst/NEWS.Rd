\name{diffHicnews}
\title{diffHic News}
\encoding{UTF-8}

\section{Version 1.5.5}{\itemize{
\item Relaxed checks in preparePairs(), prepPseudoPairs() when more chromosomes are supplied than needed.

\item Fixed an issue with preparePairs() when too many file handles are open.

\item Fixed clash between BiocGenerics and Matrix which().

\item Modified domainDirections() to return a RangedSummarizedExperiment.

\item Removed defunct DIList class and methods.
}
}

\section{Version 1.4.0}{\itemize{
\item 
Deprecated DIList objects and methods in favour of InteractionSet objects.

\item
marginCounts() now returns a RangedSummarizedExperiment for all bins.

\item
Added the max.height argument to the rotPlaid() and rotDI() functions.

\item
Added the diClusters() function for post-hoc cluster-level FDR control.

\item
Added the annotatePairs() function for convenient annotation of (clusters of) interactions.

\item
Fixed a bug in plotPlaid() when the interaction space was empty.

\item
Fixed a bug in preparePairs() where unmapped chimeric segments led to the loss of the entire pair.

\item
Updated user's guide, documentation and tests.
}
}

\section{Version 1.2.0}{\itemize{
\item
Renamed normalize() to normOffsets().

\item
Added library size specification to DIList methods normOffsets(), asDGEList().

\item 
Fixed bugs under pathological settings in plotPlaid(), plotDI(), rotPlaid(), rotDI().

\item 
Optimized C++ code for connectCounts(), squareCounts().

\item 
Streamlined various R utilities used throughout all functions.

\item 
Added iter_map.py to inst/python, for iterative mapping of DNase Hi-C data.

\item
Added the neighborCounts() function, for simultaneous read counting and enrichment calculation.

\item
Added exclude for enrichedPairs(), to provide an exclusion zone in the local neighborhood.

\item
Switched default colour in rotPlaid(), plotPlaid() to black.

\item
Added compartmentalize() function to identify genomic compartments.

\item 
Added domainDirections() function to help identify domains.

\item
Modified correctedContact() to allow distance correction and report factorized probabilities directly.

\item 
Modified marginCounts() function for proper single-end-like treatment of Hi-C data.

\item
Extended clusterPairs() to merge bin pairs from multiple DILists.

\item
Switched to reporting ranges directly from boxPairs(), added support for minimum bounding box output.

\item
Modified consolidatePairs() to accept index vectors for greater modularity.

\item 
Added reference argument for large bin pairs, in filterDirect() and filterTrended().

\item
Added filterDiag() convenience function for filtering of (near-)diagonal bin pairs.

\item 
Slight change to preparePairs() diagnostic reports when dedup=FALSE, and for unpaired reads.

\item
Added option for a distance-based threshold to define invalid chimeras in preparePairs().

\item 
Updated documentation, tests and user's guide.

\item
Added diffHic paper entry to CITATION.
}}

\section{Version 1.0.0}{\itemize{
\item 
New package \code{diffHic}, for detection of differential interactions from Hi-C data.
}}
