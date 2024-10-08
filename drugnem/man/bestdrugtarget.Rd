% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drugmnemscoring.R
\name{bestdrugtarget}
\alias{bestdrugtarget}
\title{Scoring functions for all drug combinations using the nested effects models}
\usage{
bestdrugtarget(G, mappos, data)
}
\arguments{
\item{G}{a graphnel object derived from nem.}

\item{mappos}{List of names for each drug node in the network corresponding to the proteins connected to each node also derived from a nem output.}

\item{data}{2D matrix of posterior probabilities for each drug-protein pair derived from a nem output object.}
}
\value{
List of the drug regimen(s) with the highest score and a ranking of all drug combinations and their corresponding scores:
\item{Ordercomb}{Ranking of all drug combinations and their corresponding scores}
}
\description{
This function takes a nem graph object on drugs, a list of Drug tagets for each node and a data matrix with posterior weights for attaching a target protein to each node in the drug network. It produes a list of drug regimen(s) with the best score and a ranking of all drug combinations and their corresponding scores.
}
\examples{
out22=bestdrugtarget(G=rr$graph,mappos=mappos,data=pos)

}
\seealso{
\code{\link{drugmnemmain}} which this function wraps
}
