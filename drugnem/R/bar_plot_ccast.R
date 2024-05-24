#' Plotting Barplots for Homogeneous subpopulations
#'
#' This function takes as input 2 matrices corresponding to means and
#' standard deviations respectively of dimension mxn each across all markers
#' of size n and produces a pdf file of barplots for all subpopulation of
#' size m corresponding to the leaf of the ccast tree.
#'
#' @param M mxn matrix of means for m subpopulations and n markers
#' @param M2 mxn matrix of standard deviations for m subpopulations and n
#' markers
#' @param file name of pdf output file to be saved in the current directory
#' @return A pdf file labeled "Barplot for homogenous cells.pdf" with
#' barplots for all homogeneous subgroups.
#' @seealso \code{\link{ccast_filter5}} which this function wraps
#' @export
#' @examples
#' bar_plot_ccast(M=M, M2=M2, "Barplot for homogenous cells.pdf")
#'
bar_plot_ccast <-
function(M,M2,file) {
  pdf(file,width=15, height=10)
  u1<-ifelse(min(M)<0,(min(M)-1),min(M))
  u2<-max(M)+1
  ###par(mfrow=c(cc1,cc2))
  for (i in 1:dim(M)[1]) {
    mp<-barplot(M[i,]+u1,names.arg=colnames(M),main=rownames(M)[i],ylim=c(u1-2,(u2+1)),yaxt="n",yaxs="r",family= "sans",space=c(1,1),xlab="",ylab="",cex.lab=2)
    axis(1,mp, labels = FALSE, tick = T)
    axis(2,M[i,]+u1,labels = FALSE, tick = T)
    arrows(mp,M[i,]+u1,mp,(M[i,]+u1)+ M2[i,],length=0.1, angle = 90)
    arrows(mp,M[i,]+u1,mp,(M[i,]+u1)- M2[i,], length=0.1,angle = 90)
    mtext(min(M[i,]),side = 2, at =min(M[i,])+u1 ,family= "sans", line = 1)
    
  }
  dev.off()
}
