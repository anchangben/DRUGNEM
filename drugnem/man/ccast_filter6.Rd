% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccast_filter6.R
\name{ccast_filter6}
\alias{ccast_filter6}
\title{A function for clustering, classification and sorting tree that identifies and isolates homogeneous cell subpopulations from a heterogeneous single cell data using a decision tree representation.}
\usage{
ccast_filter6(file, asinhp = asinhp, colid, coln = NULL, rown = NULL, k,
  patient)
}
\arguments{
\item{file}{A file path or directory path name with all fcs files of interest or a cxm flow cytometry expression matrix of c number of cells and m number of markers.}

\item{asinhp}{asinh parameter p for scaling the data if necessary. Default is 1.}

\item{colid}{index of markers of interest to be used by ccast algorithm.}

\item{coln}{column names of markers to use of equal length as colid. If not provided the default fcs column names will be used.}

\item{rown}{row names of data to use of equal length as number of cells. If not provided the default fcs file names will be used.}

\item{k}{number of expected clusters for hclust.vector algorithms.}

\item{patient}{Patient or subject name. Name must be included in the file names.}
}
\value{
Several pdf files corresponding to various diagnostic plots including plots of initial and final optimized ccast decision  tree In addition to a list 8 items:
\item{initialtree }{initial ccast tree output with less homogeneous bins in the leaf nodes}
\item{finaltree }{final ccast tree output with distinct homogeneous bins in the leaf nodes}
\item{Initialdata}{Initial unfiltered data matrix with only selected markers}
\item{Finaldata}{Final filtered data matrix with only selected markers}
\item{allfinalexprdata}{Final filtered data matrix with all markers}
\item{allorigdata}{Original data matrix with all markers}
\item{treeheight}{Final tree height}
\item{responderids}{list of row ids(cell ids) corresponding to responders and non-responders. }
}
\description{
This function uses hirerachical clustering approach to quantify underlying homogeneous subpopulations of cells, then recursive partitioning techniques to generate the decision tree comprising the gating strategy for all subpopulations of interest. It further optimizes the decision tree to produce  homogeneous cell states or subgroups.
}
\examples{
ccast_filter6(file,asinhp=asinhp,colid,coln=NULL,rown=NULL,k,patient)
}
\seealso{
\code{\link{responders}} which this function wraps
}
