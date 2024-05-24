#' Plotting Heatmaps for Homogeneous subpopulations
#'
#' This function takes as input an optimal ccast decision tree output from the ccast_optimize function and produces 2 matrices corresponding to means and standard deviations respectively of dimension mxn each across all markers of size n and produces a pdf file of heatmaps for all subpopulation of size m.
#'
#' @param optccastreeoutput An optimal ccast tree(decision tree) object produced by the function ccast_optimize
#' @return an R objects corresponding to (1) a list of row ids(cell ids) (2) Tree node ids corresponding to all Subpopulations  and a pdf file labeled "Heatmaps Homogeneous cells.pdf" of heatmaps for all leaf nodes from the ccast tree. In addition to a list of 3 items:
#' \item{ccastmeans }{mxn matrix of means for m subpopulations and n markers }
#' \item{ccastsds }{mxn matrix of standard deviations for m subpopulations and n markers}
#' \item{idcells}{list of cell ids for all lineage subpopulations}
#' @seealso \code{\link{ccast_filter5}} which this function wraps
#' @export
#' @examples
#' ccast_heatmaplot(optccastreeoutput=optccastree)
#'
ccast_heatmaplot <-
function(optccastreeoutput,file="Heatmaps Homogeneous cells.pdf"){
    pdf(file=file,width=10, height=10)
    tree3<-optccastreeoutput$tree
    tt1=party::nodes(tree3, unique(where(tree3)))
    dD3=optccastreeoutput$puredata[[length(optccastreeoutput$puredata)]]
    ##require(RColorBrewer)
    nrcolors = 100
    half = 1 + nrcolors/2
    colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
    colorpalette = colorRampPalette(colpal)(nrcolors)
    
    ##layout(matrix(1:length(tt1), ncol = round(sqrt(length(tt1)),0)))
    MM <- NULL
    MM2 <- NULL
    rname<-NULL
    rowids<-list()
    ###par(mfrow=c(4,2)))
    homodata=list()
    for( j in 1:length(tt1)) {
        mx=which.max(tt1[[j]]$prediction)
        w1=which(tt1[[j]]$weights==1)
        w2=w1[which(dD3$groups[w1]==mx)]
        print(table(dD3$groups[w2]))
        minv=min(dD3[,-dim(dD3)[2]],na.rm=TRUE)
        maxv=max(dD3[,-dim(dD3)[2]],na.rm=TRUE)
        key_color1=seq(minv,maxv,length.out=dim(dD3[,-dim(dD3)[2]])[2])
        
        rowids[[j]]=as.numeric(rownames(dD3[w2,-dim(dD3)[2]]))
        BB=as.matrix(rbind(dD3[w2,-dim(dD3)[2]],key_color1))
        rownames(BB)=rep(paste("Cell",j),dim(BB)[1])
        mm <- round(c(apply(BB[-dim(BB)[1],], 2, mean)),3)
        MM <- rbind(MM,mm)
        mm2 <- round(c(apply(BB[-dim(BB)[1],], 2, sd)),3)
        MM2 <- rbind(MM2,mm2)
        image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = paste("Celltype",mx),xaxt="n",yaxt="n",ylab = "Cells",family="sans",cex.lab=1.4,main=paste("Node",tt1[[j]]$ nodeID))
        mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=0.8)
        mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=1,line=1, family= "sans", cex=0.8)
        homodata[[j]]=BB[-dim(BB)[1],]
        rownames(homodata[[j]])=paste(rowids[[j]])
        rname=c(rname,paste("Node",tt1[[j]]$ nodeID))
        
    }
    rownames(MM)=rname
    rownames(MM2)=rname
    save(homodata,file="combinedatalist.rdata")
    save(rname,file="NodeIDs.rdata")
    save(rowids,file="rowids.rdata")
    dev.off();
   return(list(ccastmeans=MM,ccastsds=MM2,idcells=rowids))
}
