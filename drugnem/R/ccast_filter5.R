#' A function for clustering, classification and sorting tree that identifies and isolates homogeneous cell subpopulations from a heterogeneous single cell data using a decision tree representation.
#'
#' This function uses  a non-parametric mixture model or hirerachical clustering approach and silhouette measures to quantify underlying homogeneous subpopulations of cells, then recursive partitioning techniques to generate the decision tree comprising the gating strategy for all subpopulations of interest. It further optimizes the decision tree to produce  homogeneous cell states or subgroups.
#'
#' @param file A file path or directory path name with all fcs files of interest or a cxm flow cytometry expression matrix of c number of cells and m number of markers.
#' @param transformlogic An character string that defines the type of transformation. See function read.FCS in R packgae flowCore for more information. Valid values are linearize (default) or scale.The linearize transformation applies the appropriate power transform to the data while the scale transformation scales all columns to $[0,10^decades]$. defaulting to decades=0 as in the FCS4 specification. A logical can also be used: TRUE is equal to linearize and FALSE(which is the default) corresponds to no transformation.
#' @param asinhp asinh parameter p for scaling the data if necessary. Default is 1.
#' @param colid index vector of markers of interest to be used by ccast algorithm.
#' @param coln column names of markers to use of equal length as colid. If not provided the default fcs column names will be used.
#' @param rown row names of data to use of equal length as number of cells. If not provided the default fcs file names will be used.
#' @param k number of expected clusters for both npEM or hclust.vector algorithms.
#' @param param string of parameter names from an FCS file. Default is NULL.
#' @param fcsdes Description of parameter names from an FCS file. Default is NULL.
#' @param ylabel y axis marker variable used to generate biaxial plots for all inner node marker variables from the optimized ccast decision tree. We recommend using a marker that has not been used by ccast to predict the decision tree.
#' @param subanalysis A logical to downsample original FCS file(s). Default is FALSE. Recommended if total number of cells across all FCS files > 100000.
#' @param patient Patient or subject name. Name must be included in the file names.
#' @param subsamplesize Number of cells to subsample using density downsampling approach by Qiu et al. 2011
#' @return 	Several pdf files corresponding to various diagnostic plots including plots of initial and final optimized ccast decision  tree In addition to a list 7 items:
#' \item{initialtree }{initial ccast tree output with less homogeneous bins in the leaf nodes}
#' \item{finaltree }{final ccast tree output with distinct homogeneous bins in the leaf nodes}
#' \item{Initialdata}{Initial unfiltered data matrix with only selected markers}
#' \item{Finaldata}{Final filtered data matrix with only selected markers}
#' \item{allfinalexprdata}{Final filtered data matrix with all markers}
#' \item{allorigdata}{Original data matrix with all markers}
#' \item{treeheight}{Final tree height}
#' @seealso \code{\link{drugmnemmain}} which this function wraps
#' @export
#' @examples
#' ccast_filter5(file,transformlogic=FALSE,asinhp=asinhp,colid,coln=NULL,rown=NULL,k,param=NULL,fcsdes=NULL,ylabel,subanalysis=FALSE,patient,subsamplesize=subsamplesize)

ccast_filter5 <-
function(file,transformlogic=FALSE,asinhp=asinhp,colid,coln=NULL,rown=NULL,k,param=NULL,fcsdes=NULL,ylabel,subanalysis=FALSE,patient,subsamplesize=subsamplesize){
    
    ##require(flowCore)
    
    if(is.matrix(file)) { alldata=file
    } else {
        Dall1=fcsimport3(file=file,transformlogic=transformlogic,asinhp=asinhp,colid=colid,coln=coln,rown=rown,patient,subanalysis=subanalysis,subsamplesize=subsamplesize)
        Dall=Dall1[[1]]
        fcsfile=Dall1[[3]]
        alldata=Dall1[[2]]
        rownames(alldata)=paste(1:dim(alldata)[1])
        save(alldata,file="Allexprsdata.rda")
    }
    
    ##require(fastcluster)
    ##require(cluster)
    
    if(dim(alldata)[1]<= 25000){
        if(is.null(asinhp)) Dallf=alldata[,colid,drop=FALSE]
    else Dallf=asinh(alldata[,colid]/asinhp)
    colnames(Dallf)=coln
    if((is.null(rownames(Dallf)))|(length(unique(rownames(Dallf)))==1)) rownames(Dallf)=paste(1:dim(Dallf)[1])
        hc = hclust.vector(Dallf, method="ward")
        d <- dist(Dallf, method = "euclidean")
        assign <- cutree(hc, k=k)
        sk <- silhouette(assign, d)
       # pdf(file=paste("All data sillhouette plots.pdf"),width=12, height=15)
       # plot(sk,main="All data")
       # dev.off();
        id=which(sk[,3]<=0)
        Dall3=Dallf[-id,]
        d2 <- dist(Dall3, method = "euclidean")
        hc2 = hclust.vector(Dall3, method='ward')
        res2 <- cutree(hc2, k=k)
        sk0 = silhouette(res2, d2)
        if(all(sk0[,3] >= 0)) {
        #pdf(file=paste("Filtered sillhouette plot all data.pdf"),width=12, height=15)
            #plot(sk,main="All data filtered")
            #dev.off();
        groups=res2
        Dall=Dall3
        save(groups,file="groups.rda")
        } else {
        out1=ccast_silhouette(d2,fileout="silhouette.pdf", nCluster.init = k, assign=res2, iter.max = 10, diagnosis = TRUE)
        groups=out1$optFitCluster[["cluster"]]
        Dall=Dall3
        save(groups,file="groups.rda")
        }
    
      } else {
    if(is.null(asinhp)) Dallf=alldata[,colid,drop=FALSE]
      Dallf=asinh(alldata[,colid]/asinhp)
      colnames(Dallf)=coln
      rownames(Dallf)=paste(1:dim(Dallf)[1])               
        Dall<- Dallf
        hc = hclust.vector(Dall, method='ward')
        groups <- cutree(hc, k=k)
        save(Dall,file="Cleandata.rda")
        save(groups,file="groups.rda")
        }
    
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
ccast_biaxialplot(optccastreeoutput=optccastree,ylabel)   
#### Heatmaps and barplots for Homogenous subpopulations ###########
###############################
heatmapdata<-ccast_heatmaplot(optccastreeoutput=optccastree)
#### Barplots for Homogenous subpopulations ###########
###############################
    MM=heatmapdata[[1]]
    MM2=heatmapdata[[2]]
    bar_plot_ccast(MM, MM2, "Barplot for homogenous cells.pdf")
#### Pairwise Scatterplots for original and final homogenous data for major markers ###########
###############################
####scatterplot_ccast(data= dD1,optccastreeoutput=optccastree)
if((is.null(param))& (is.null(fcsdes))) {
param=flowCore::parameters(fcsfile)
cDD=alldata[as.numeric(rownames(Dall)),]
out_frame <- flowFrame(cDD, param, description = description(fcsfile))
finalD<-optccastree$puredata[[length(optccastree$puredata)]]
finalfcsD<-alldata[as.numeric(rownames(finalD)),]
out_frame2 <- flowFrame(finalfcsD, param, description = description(fcsfile))
} else {
cDD=alldata
out_frame <- flowFrame(cDD, param, description = fcsdes)
finalD<-optccastree$puredata[[length(optccastree$puredata)]]
finalfcsD<-alldata[as.numeric(rownames(finalD)),]
out_frame2 <- flowFrame(finalfcsD, param, description = fcsdes)
}
write.FCS(out_frame, paste(getwd(),"origfiltereddata.fcs",sep="/"))
write.FCS(out_frame2, paste(getwd(),"finalfiltereddata.fcs",sep="/"))

return(list(initialtree=ccastree[[3]],finaltree=optccastree[[1]],Initialdata=Dall,Finaldata=finalD,allfinalexprdata=finalfcsD,allorigdata=alldata,treeheight=s))
}
