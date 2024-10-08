% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccast_filter5.R
\name{ccast_filter5}
\alias{ccast_filter5}
\title{A function for clustering, classification and sorting tree that identifies and isolates homogeneous cell subpopulations from a heterogeneous single cell data using a decision tree representation.}
\usage{
ccast_filter5(file, transformlogic = FALSE, asinhp = asinhp, colid,
  coln = NULL, rown = NULL, k, param = NULL, fcsdes = NULL, ylabel,
  subanalysis = FALSE, patient, subsamplesize = subsamplesize)
}
\arguments{
\item{file}{A file path or directory path name with all fcs files of interest or a cxm flow cytometry expression matrix of c number of cells and m number of markers.}

\item{transformlogic}{An character string that defines the type of transformation. See function read.FCS in R packgae flowCore for more information. Valid values are linearize (default) or scale.The linearize transformation applies the appropriate power transform to the data while the scale transformation scales all columns to $[0,10^decades]$. defaulting to decades=0 as in the FCS4 specification. A logical can also be used: TRUE is equal to linearize and FALSE(which is the default) corresponds to no transformation.}

\item{asinhp}{asinh parameter p for scaling the data if necessary. Default is 1.}

\item{colid}{index vector of markers of interest to be used by ccast algorithm.}

\item{coln}{column names of markers to use of equal length as colid. If not provided the default fcs column names will be used.}

\item{rown}{row names of data to use of equal length as number of cells. If not provided the default fcs file names will be used.}

\item{k}{number of expected clusters for both npEM or hclust.vector algorithms.}

\item{param}{string of parameter names from an FCS file. Default is NULL.}

\item{fcsdes}{Description of parameter names from an FCS file. Default is NULL.}

\item{ylabel}{y axis marker variable used to generate biaxial plots for all inner node marker variables from the optimized ccast decision tree. We recommend using a marker that has not been used by ccast to predict the decision tree.}

\item{subanalysis}{A logical to downsample original FCS file(s). Default is FALSE. Recommended if total number of cells across all FCS files > 100000.}

\item{patient}{Patient or subject name. Name must be included in the file names.}

\item{subsamplesize}{Number of cells to subsample using density downsampling approach by Qiu et al. 2011}
}
\value{
Several pdf files corresponding to various diagnostic plots including plots of initial and final optimized ccast decision  tree In addition to a list 7 items:
\item{initialtree }{initial ccast tree output with less homogeneous bins in the leaf nodes}
\item{finaltree }{final ccast tree output with distinct homogeneous bins in the leaf nodes}
\item{Initialdata}{Initial unfiltered data matrix with only selected markers}
\item{Finaldata}{Final filtered data matrix with only selected markers}
\item{allfinalexprdata}{Final filtered data matrix with all markers}
\item{allorigdata}{Original data matrix with all markers}
\item{treeheight}{Final tree height}
}
\description{
This function uses  a non-parametric mixture model or hirerachical clustering approach and silhouette measures to quantify underlying homogeneous subpopulations of cells, then recursive partitioning techniques to generate the decision tree comprising the gating strategy for all subpopulations of interest. It further optimizes the decision tree to produce  homogeneous cell states or subgroups.
}
\examples{
ccast_filter5(file,transformlogic=FALSE,asinhp=asinhp,colid,coln=NULL,rown=NULL,k,param=NULL,fcsdes=NULL,ylabel,subanalysis=FALSE,patient,subsamplesize=subsamplesize)
}
\seealso{
\code{\link{drugmnemmain}} which this function wraps
}
