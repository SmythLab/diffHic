\name{diffHicnews}
\title{diffHic News}
\encoding{UTF-8}

\section{Version 1.1.17}{\itemize{
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
