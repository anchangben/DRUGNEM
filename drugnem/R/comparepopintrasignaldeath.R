#' A function to compare intracellular or functional signaling across pre defined death cell populations.
#'
#' This function takes as input treatment response single-cell expression data labeled by cell types and returns  a list of dataframes and a matrix with relative percentage of death cells for all drugs across all subpopulations.
#'
#' @param rowids List of row ids(cell ids) for each sub population of cells stored as R object.
#' @param dataP Marker expression data matrix of pooled cells from original FCS file stored as R object.
#' @param asinhp asinh parameter p for scaling the data if necessary. Default is 1.
#' @param drugs Parameter vector corresponding to the names of all single drugs. 
#' @param sigids Vector of marker ids corresponding to the functional proteins used for building the Nested effects models.
#' @param colid column ids of lineage markers to use of equal length as colid. Must be provided by user.
#' @param antibody Vector of marker names of length equal to length of the columns in each FCS file for each channel in the FCS file.
#' @param patient Patient or subject name. Name must be included exactly as in the file names. 

#' @return 	List of dataframes with mean and fold change expression values for each drug relative to the baseline and heatmaps of fold changes with relative percentage of cells compared to basal for all drugs. The output of this functionis a list 3 items:
#'\item{initialtree }{initial ccast tree output with less homogeneous bins in the leaf nodes}
#' \item{singlecelldata}{List of dataframes for each subpopulation with drug response values for each cell across functional proteins }
#' \item{Popmeanexprsdata}{List of dataframes for each subpopulation with mean expression for proteins(rows) and drugs (columns).}
#' \item{Diffpopmeanexprsdata}{List of dataframes for each subpopulation with Fold change expression for proteins(rows) and drugs (columns) relative to the basal treatment.}
# \item{deathdist}{Matrix of percentages of cells in all death subpopulations(rows) for drugs (columns).} 
#' @seealso \code{\link{drugmnemmain}} which this function wraps
#' @export
#' @examples
#' res2=comparepopintrasignaldeath(rowids2,dataP,asinhp=asinhp,drugs=drugs,sigids=sigids,colid=colid,antibody=antibody,patient=patient)
#'

##### Comparison of intracellular signaling across death subpopulations ########
###### count percentage of death cells in sub populations#####
comparepopintrasignaldeath<-function(rowids,dataP,asinhp,drugs,sigids,colid,antibody,patient) {
    
    nrcolors = 100
    half = 1 + nrcolors/2
    colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
    colorpalette = colorRampPalette(colpal)(nrcolors)
    Tot=length(unlist(rowids))
    sids=dataP[,dim(dataP)[2]][unlist(rowids)]
    tt=table(sids)
    total=NULL
    for (ii in 1:length(tt)) total[ii]=tt[ii]
    if(is.null(asinhp)) U=dataP[,colid]
    else U=asinh(dataP[,colid]/asinhp)
    
    minu=min(U,na.rm=TRUE)
    maxu=max(U,na.rm=TRUE)
    minx=max(abs(minu),abs(maxu))
    #key_color1=seq(-1*minx,minx,length.out=length(colid))
    key_color1=seq(minu,maxu,length.out=length(colid))
    
    if(is.null(asinhp)) Z=dataP[,sigids]
    else Z=asinh(dataP[,sigids]/asinhp)
    
    minz=min(Z,na.rm=TRUE)
    maxz=max(Z,na.rm=TRUE)
    miny=max(abs(minz),abs(maxz))
    #key_color2=seq(-1*miny,miny,length.out=length(sigids))
    key_color2=seq(minz,maxz,length.out=length(sigids))
    
    Sdata=list()
    Sdata2=list()
    cellstate <- paste("Pop",1:length(rowids))
    pProteins=antibody[sigids]
    dimn=list(pProteins,drugs,cellstate)
    Out=array(NA,c(length(pProteins),length(drugs),length(cellstate)),dimnames=dimn)
    psize2=list()
    percentmat=matrix(0,nrow=length(rowids),ncol=length(drugs))
    colnames(percentmat)=drugs
    rownames(percentmat)=cellstate
    for( j in 1:length(rowids)) {
        
        rowids0=dataP[,dim(dataP)[2]][rowids[[j]]]
        pdf(file=paste(patient,"Heatmaps for matching Treatments in Death Pop", j,".pdf",sep=""),width=12, height=12)
        
        Sdata[[j]]=list()
        pdc=NULL
        psize2[[j]]=list()
        for ( k in (as.numeric(names(table(rowids0))))) {
            rw1=which(rowids0==k)
            rowids1=rowids[[j]][rw1]
            
            BB1=as.matrix(rbind(U[rowids1,],key_color1))
            BB3=as.matrix(rbind(Z[rowids1,],key_color2))
            
            rownames(BB1)=paste(c(rowids1,"key"))
            colnames(BB1)=antibody[colid]
            rownames(BB3)=paste(c(rowids1,"key"))
            colnames(BB3)=antibody[sigids]
            
            Sdata[[j]][[k]]=BB3[-dim(BB3)[1],,drop = FALSE]
            
            t1=apply(Sdata[[j]][[k]],2,mean)
            t2=t(Sdata[[j]][[k]])
            colnames(t2)=rep(drugs[k],dim(t2)[2])
            pdc=cbind(pdc,t2)
            Out[,k,j]=t1
            
            par(mfrow=c(2,1))
            
            psize2[[j]][[k]]=round(((dim(BB1)[1]-1)/total[k])*100,4)
            percentmat[j,k]=psize2[[j]][[k]]
            image(x=1:dim(BB1)[2],y=1:dim(BB1)[1],z=as.matrix(t(BB1)),col = colorpalette,xlab = paste("Surface Markers"),xaxt="n",yaxt="n",ylab = "Cells",family="sans",cex.lab=0.8,main=paste("Death Pop",j,drugs[k],"with Relative % =",psize2[[j]][[k]]))
            mtext(rownames(BB1),at=c(1:dim(BB1)[1]),side=4,las=2,line=1, family= "sans", cex=0.4)
            mtext(colnames(BB1),at=c(1:dim(BB1)[2]),side=1,las=1,line=1, family= "sans", cex=0.4)
            
            
            image(x=1:dim(BB3)[2],y=1:dim(BB3)[1],z=as.matrix(t(BB3)),col = colorpalette,xlab = paste("Intracellular Markers"),xaxt="n",yaxt="n",ylab = "Cells",family="sans",cex.lab=0.8,main=paste("Death Pop",j,drugs[k],"with Relative % =",psize2[[j]][[k]]))
            mtext(rownames(BB3),at=c(1:dim(BB3)[1]),side=4,las=2,line=1, family= "sans", cex=0.4)
            mtext(colnames(BB3),at=c(1:dim(BB3)[2]),side=1,las=1,line=1, family= "sans", cex=0.4)
            
        }
        Sdata2[[j]]=pdc
        
        dev.off();
    }
    
    #save(Sdata, file="Sdata.rdata")
    #save(Out,file="Out.rdata")
    #save(Sdata2, file="Sdata2.rdata")
    save(percentmat, file="percentmat.rdata")
    #### Plot population mean expression profiles ####
    #pdf(file=paste(patient,"Heatmap of Mean FC Protein Expressions across Death cell states.pdf",sep=""),width=15, height=15)
    #for (t in 1:length(cellstate)) {
     #   if(sum(is.na(Out[,1,t]))>0) {
          #  BB=NULL
      #  } else {
            
       #     BB=Out[,,t]-Out[,1,t]
        #    minzz=min(BB,na.rm=TRUE)
        #    maxzz=max(BB,na.rm=TRUE)
         #   minyy=max(abs(minzz),abs(maxzz))
         #   key_color=seq(-1*minyy,minyy,length.out=dim(BB)[2])
        #    BB=as.matrix(rbind(BB,key_color))
        #    par(mar=c(3,3,3,8))
        #    image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = #"Proteins",family="sans",cex.lab=1.4,main=paste(cellstate[t],"Protein Expressions"))
    #        mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=1.5)
      #      mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=1,line=1, family= "sans", cex=1.5)
       # }
  #  }
   # dev.off();
    
    ##### Global normalization ################
   # pdf(file=paste(patient,"Heatmap of Global Mean FC Protein Expressions across Death cell states.pdf",sep=""),width=15, height=15)
    Out2=Out
    for (t in 1:length(cellstate)) {
        if(sum(is.na(Out2[,1,t]))>0) {
            BB=NULL
        } else {
            Out2[,,t]=Out2[,,t]-Out2[,1,t]
        }
    }
    mino=min(Out2,na.rm=TRUE)
    maxo=max(Out2,na.rm=TRUE)
    minoo=max(abs(mino),abs(maxo))
    key_color1=seq(-1*minoo,minoo,length.out=dim(Out2)[2])
    ##key_color1=seq(mino,maxo,length.out=dim(Out2)[2])
   # for (t in 1:length(cellstate)) {
    #    BB=Out[,,t]-Out[,1,t]
      #  BB=as.matrix(rbind(BB,key_color1))
      #  par(mar=c(3,3,3,8))
       # image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = #"Proteins",family="sans",cex.lab=1.4,main=paste(cellstate[t],"Protein Expressions"))
     #   mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=1.5)
      #  mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=1,line=1, family= "sans", cex=1.5)
    #}
   # dev.off();
    ## save(Out2,file="Out2.rdata")
    ##### RUN MNEM ###########
    ##detach("package:party", unload=TRUE)
    ###fitmnemdown(Out3=Sdata2,Outp=Out2,drugs,p=NULL,patient=patient,infer=infer,type="CONTmLL")
    
    return(list(singlecelldata=Sdata,Popmeanexprsdata=Out,Diffpopmeanexprsdata=Out2, deathdist=percentmat))
    
}
######
