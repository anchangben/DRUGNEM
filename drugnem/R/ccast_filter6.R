#' A function for clustering, classification and sorting tree that identifies and isolates homogeneous cell subpopulations from a heterogeneous single cell data using a decision tree representation.
#'
#' This function uses hirerachical clustering approach to quantify underlying homogeneous subpopulations of cells, then recursive partitioning techniques to generate the decision tree comprising the gating strategy for all subpopulations of interest. It further optimizes the decision tree to produce  homogeneous cell states or subgroups.
#'
#' @param file A file path or directory path name with all fcs files of interest or a cxm flow cytometry expression matrix of c number of cells and m number of markers.
#' @param asinhp asinh parameter p for scaling the data if necessary. Default is 1.
#' @param colid index of markers of interest to be used by ccast algorithm.
#' @param coln column names of markers to use of equal length as colid. If not provided the default fcs column names will be used.
#' @param rown row names of data to use of equal length as number of cells. If not provided the default fcs file names will be used.
#' @param k number of expected clusters for hclust.vector algorithms.
#' @param patient Patient or subject name. Name must be included in the file names.
#' @return 	Several pdf files corresponding to various diagnostic plots including plots of initial and final optimized ccast decision  tree In addition to a list 8 items:
#' \item{initialtree }{initial ccast tree output with less homogeneous bins in the leaf nodes}
#' \item{finaltree }{final ccast tree output with distinct homogeneous bins in the leaf nodes}
#' \item{Initialdata}{Initial unfiltered data matrix with only selected markers}
#' \item{Finaldata}{Final filtered data matrix with only selected markers}
#' \item{allfinalexprdata}{Final filtered data matrix with all markers}
#' \item{allorigdata}{Original data matrix with all markers}
#' \item{treeheight}{Final tree height}
#' \item{responderids}{list of row ids(cell ids) corresponding to responders and non-responders. } 
#' @seealso \code{\link{responders}} which this function wraps
#' @export
#' @examples
#' ccast_filter6(file,asinhp=asinhp,colid,coln=NULL,rown=NULL,k,patient)

ccast_filter6 <-
function(file,asinhp=asinhp,colid,coln=NULL,rown=NULL,k,patient){
    
       alldata=file
    if(is.null(asinhp)) Dallf=alldata[,colid,drop=FALSE]
      else Dallf=asinh(alldata[,colid,drop=FALSE]/asinhp)
      colnames(Dallf)=coln
      if((is.null(rownames(Dallf)))|(length(unique(rownames(Dallf)))==1)) rownames(Dallf)=paste(1:dim(Dallf)[1])            
        Dall<- Dallf
        hc = hclust.vector(Dall, method='ward')
        groups <- cutree(hc, k=k)
          if(is.numeric(Dall)) {
              Dall=matrix(Dall,ncol=1)
              colnames(Dall)=coln
              rownames(Dall)=names(Dallf)
          }
        save(Dall,file="Cleandata.rda")
        save(groups,file="groups.rda")
      tabg=table(groups)  
    if (any(tabg<2)) print(paste("Little evidence of death cells in Subpopulation"))
    
### Create a dataframe of data and predicted cluster assignments ######
    D3=cbind(Dall,groups)
    dD1=as.data.frame(D3)
   
##Determining CCAST tree and prunning level L of CCAST tree #########
    mm=max(groups)
    ccastree<-ccast_L(dD1,mm)

##Optimizing CCAST tree #########
  t1=ccastree[[1]]
  s=ccastree[[2]]
  optccastree<-ccast_optimize(t1,dD1,s)

##### Plotting CCAST diagnostic output #######
#### plot splitpoint statistics ###########
############## Biaxial plot of tree #########
##ccast_biaxialplot(optccastreeoutput=optccastree,ylabel)
#### Heatmaps and barplots for Homogenous subpopulations ###########
###############################
heatmapdata<-ccast_heatmaplot2(optccastreeoutput=optccastree)
#### Barplots for Homogenous subpopulations ###########
###############################
    MM=heatmapdata[[1]]
    MM2=heatmapdata[[2]]
    idrows=heatmapdata[[3]]
   ## bar_plot_ccast(MM, MM2, "Barplot for homogenous cells.pdf")
#### Pairwise Scatterplots for original and final homogenous data for major markers ###########
###############################
####scatterplot_ccast(data= dD1,optccastreeoutput=optccastree)

cDD=alldata[as.numeric(rownames(Dall)),]
finalD<-optccastree$puredata[[length(optccastree$puredata)]]
   matchids=match(as.numeric(rownames(finalD)),rownames(cDD))
finalfcsD<-alldata[as.numeric(matchids),,drop=FALSE]


return(list(initialtree=ccastree[[3]],finaltree=optccastree[[1]],Initialdata=Dall,Finaldata=finalD,allfinalexprdata=finalfcsD,allorigdata=alldata,treeheight=s,responderids=idrows))

}
