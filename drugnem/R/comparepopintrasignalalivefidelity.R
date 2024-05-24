#' A function to compare both lineage and intracellular signaling changes across pre defined survival cell populations.
#'
#' This function takes as input treatment response single-cell expression data labeled by cell types and returns  a list of dataframes and heatmaps of fold changes with relative percentage of cells compared to basal for all drugs.
#'
#' @param rowids List of row ids(cell ids) for each sub population of cells stored as R object.
#' @param dataP Marker expression data matrix of pooled cells from original FCS file stored as R object.
#' @param asinhp asinh parameter p for scaling the data if necessary. Default is 1.
#' @param drugs Parameter vector corresponding to the names of all single drugs. 
#' @param sigids Vector of marker ids corresponding to the functional proteins used for building the Nested effects models.
#' @param colid column ids of lineage markers to use of equal length as colid. Must be provided by user.
#' @param antibody Vector of marker names of length equal to length of the columns in each FCS file for each channel in the FCS file.
#' @param patient Patient or subject name. Name must be included exactly as in the file names. 
#' @return 	List of dataframes with mean and fold change expression values for each drug relative to the baseline and heatmaps of fold changes with relative percentage of cells compared to basal for all drugs. The output of this functionis a list 5 items:
#' \item{singlecelldata}{List of dataframes for each subpopulation with drug response expression values for each cell across intracellular proteins }
#' \item{Popmeanexprsdata}{3D array for each subpopulation(3rd dimension) with mean expression for proteins(rows) and drugs (columns).}
#' \item{Diffpopmeanexprsdata}{3D array for each subpopulation(3rd dimension) with mean FC expression for lineage proteins(rows) and drugs (columns) relative to the basal treatment.}
#' \item{Difflineagepopmeanexprsdata}{3D array for each subpopulation(3rd dimension) with mean Fold change expression for lineage proteins(rows) and drugs (columns) relative to the basal treatment.}
#' \item{unstablepop}{Numeric vector corresponding to the lineage subpopulations likely affected by drugs.}
#' @seealso \code{\link{drugmnemmain}} which this function wraps
#' @export
#' @examples
#' res2=comparepopintrasignalalivefidelity(rowids2,dataP,asinhp=asinhp,drugs=drugs,sigids=sigids,colid=colid,antibody=antibody,patient=patient)
#'

##### Comparison of intracellular signaling across survival populations ########
comparepopintrasignalalivefidelity<-function(rowids,dataP,asinhp,drugs,sigids,colid,antibody,patient) {
    
    nrcolors = 100
    half = 1 + nrcolors/2
    colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
    colorpalette = colorRampPalette(colpal)(nrcolors)
    Tot=length(unlist(rowids))
    sids=dataP[,dim(dataP)[2]][unlist(rowids)]
    tt=table(sids)
    total=NULL
    for (ii in 1:length(tt)) total[ii]=tt[ii]
    
    if(is.null(asinhp)) { U=dataP[,colid]
    }else{
      U=asinh(dataP[,colid]/asinhp)
    }
    minu=min(U,na.rm=TRUE)
    maxu=max(U,na.rm=TRUE)
    minx=max(abs(minu),abs(maxu))
    #key_color1=seq(-1*minx,minx,length.out=length(colid))
    key_color1=seq(minu,maxu,length.out=length(colid))
    
    if(is.null(asinhp)) {
        Z=dataP[,sigids]
    }else{
        Z=asinh(dataP[,sigids]/asinhp)
    }
    
    minz=min(Z,na.rm=TRUE)
    maxz=max(Z,na.rm=TRUE)
    miny=max(abs(minz),abs(maxz))
    #key_color2=seq(-1*miny,miny,length.out=length(sigids))
    key_color2=seq(minz,maxz,length.out=length(sigids))
    
    Sdata=list()
    Sdata1=list()
    Sdata2=list()
    Sdata22=list()

    cellstate <- paste("Pop",1:length(rowids))
    pProteins=antibody[sigids]
    dimn=list(pProteins,drugs,cellstate)
    Out=array(NA,c(length(pProteins),length(drugs),length(cellstate)),dimnames=dimn)
    clustmarkers=antibody[colid]
    dimn1=list(clustmarkers,drugs,cellstate)
    Out1=array(NA,c(length(clustmarkers),length(drugs),length(cellstate)),dimnames=dimn1)
    psize1=list()

    for( j in 1:length(rowids)) {
    
        rowids0=dataP[,dim(dataP)[2]][rowids[[j]]]
        print(length(rowids0))
        print(length(rowids[[j]]))

      
        pdf(file=paste(patient,"Heatmaps for matching Treatments in Pop", j,".pdf",sep=""),width=12, height=12)
        
        Sdata[[j]]=list()
        Sdata1[[j]]=list()
        pdc=NULL
        pdc1=NULL
        psize1[[j]]=list()
        print(table(rowids0))
  
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
            Sdata1[[j]][[k]]=BB1[-dim(BB1)[1],,drop = FALSE]
            
            t1=apply(Sdata[[j]][[k]],2,mean)
            t11=apply(Sdata1[[j]][[k]],2,mean)
            t2=t(Sdata[[j]][[k]])
            t22=t(Sdata1[[j]][[k]])
            colnames(t2)=rep(drugs[k],dim(t2)[2])
            colnames(t22)=rep(drugs[k],dim(t22)[2])
            pdc=cbind(pdc,t2)
            pdc1=cbind(pdc1,t22)
            Out[,k,j]=t1
            Out1[,k,j]=t11
            par(mfrow=c(2,1))
            
            psize1[[j]][[k]]=round(((dim(BB1)[1]-1)/total[k])*100,4)
            
            image(x=1:dim(BB1)[2],y=1:dim(BB1)[1],z=as.matrix(t(BB1)),col = colorpalette,xlab = paste("Surface Markers"),xaxt="n",yaxt="n",ylab = "Cells",family="sans",cex.lab=0.8,main=paste("Pop",j,drugs[k],"with Relative % =",psize1[[j]][[k]]))
            mtext(rownames(BB1),at=c(1:dim(BB1)[1]),side=4,las=2,line=1, family= "sans", cex=0.4)
            mtext(colnames(BB1),at=c(1:dim(BB1)[2]),side=1,las=1,line=1, family= "sans", cex=0.4)
            
            
            image(x=1:dim(BB3)[2],y=1:dim(BB3)[1],z=as.matrix(t(BB3)),col = colorpalette,xlab = paste("Intracellular Markers"),xaxt="n",yaxt="n",ylab = "Cells",family="sans",cex.lab=0.8,main=paste("Pop",j,drugs[k],"with Relative % =",psize1[[j]][[k]]))
            mtext(rownames(BB3),at=c(1:dim(BB3)[1]),side=4,las=2,line=1, family= "sans", cex=0.4)
            mtext(colnames(BB3),at=c(1:dim(BB3)[2]),side=1,las=1,line=1, family= "sans", cex=0.4)
            
        }
        Sdata2[[j]]=pdc
        Sdata22[[j]]=pdc1
        dev.off();
        
    }
    
    save(Sdata, file="Sdata.rdata")
    save(Sdata1, file="Sdata1.rdata")
    save(Out,file="Out.rdata")
    save(Out1,file="Out1.rdata")
    save(Sdata2, file="Sdata2.rdata")
    save(Sdata22, file="Sdata22.rdata")
    save(psize1, file="psize1.rdata")
    #### Plot subpopulation FC mean expression profiles ####
    pdf(file=paste(patient,"Heatmap of Mean FC Protein Expressions across cell states.pdf",sep=""),width=15, height=15)
    for (t in 1:length(cellstate)) {
        if(sum(is.na(Out[,1,t]))>0) {
            BB=NULL
        } else {
            
            BB=Out[,,t]-Out[,1,t]
            minzz=min(BB,na.rm=TRUE)
            maxzz=max(BB,na.rm=TRUE)
            minyy=max(abs(minzz),abs(maxzz))
            key_color=seq(-1*minyy,minyy,length.out=dim(BB)[2])
            BB=as.matrix(rbind(BB,key_color))
            par(mar=c(3,3,3,8))
            image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "Proteins",family="sans",cex.lab=1.4,main=paste(cellstate[t],"Protein Expressions"))
            mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=1.5)
            mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=1,line=1, family= "sans", cex=1.5)
        }
    }
    dev.off();
    
    pdf(file=paste(patient,"Heatmap of Mean FC Lineage Protein Expressions across cell states.pdf",sep=""),width=15, height=15)
    for (t in 1:length(cellstate)) {
        if(sum(is.na(Out1[,1,t]))>0) {
            BB=NULL
        } else {
            
            BB=Out1[,,t]-Out1[,1,t]
            minzz=min(BB,na.rm=TRUE)
            maxzz=max(BB,na.rm=TRUE)
            minyy=max(abs(minzz),abs(maxzz))
            key_color=seq(-1*minyy,minyy,length.out=dim(BB)[2])
            BB=as.matrix(rbind(BB,key_color))
            par(mar=c(3,3,3,8))
            image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "Lineage Proteins",family="sans",cex.lab=1.4,main=paste(cellstate[t],"lineage Protein Expressions"))
            mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=1.5)
            mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=1,line=1, family= "sans", cex=1.5)
        }
    }
    dev.off();
    ##### Global protein normalization ################
    #pdf(file=paste(patient,"Heatmap of Global Mean FC Protein Expressions across cell states.pdf",sep=""),width=15, height=15)
    Out2=Out
    for (t in 1:length(cellstate)) {
        if(sum(is.na(Out2[,,t]))>0) {
            Out2[,,t]=NA
        } else {
            Out2[,,t]=Out2[,,t]-Out2[,1,t]
        }
    }
    #mino=min(Out2,na.rm=TRUE)
    # maxo=max(Out2,na.rm=TRUE)
    # minoo=max(abs(mino),abs(maxo))
    # key_color1=seq(-1*minoo,minoo,length.out=dim(Out2)[2])
    ###### lineage global #####
    Out11=Out1
    for (t in 1:length(cellstate)) {
        if(sum(is.na(Out11[,,t]))>0) {
            Out11[,,t]=NA
        } else {
            Out11[,,t]=Out11[,,t]-Out11[,1,t]
        }
    }
    #mino=min(Out11,na.rm=TRUE)
    #maxo=max(Out11,na.rm=TRUE)
    # minoo=max(abs(mino),abs(maxo))
    # key_color1=seq(-1*minoo,minoo,length.out=dim(Out11)[2])

    ##key_color1=seq(mino,maxo,length.out=dim(Out2)[2])
  #  for (t in 1:length(cellstate)) {
   #     BB=Out[,,t]-Out[,1,t]
   #     BB=as.matrix(rbind(BB,key_color1))
   #     par(mar=c(3,3,3,8))
    #    image(x=1:dim(BB)[2],y=1:dim(BB)[1],z=as.matrix(t(BB)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = #"Proteins",family="sans",cex.lab=1.4,main=paste(cellstate[t],"Protein Expressions"))
   #     mtext(rownames(BB),at=c(1:dim(BB)[1]),side=4,las=2,line=1, family= "sans", cex=1.5)
   #     mtext(colnames(BB),at=c(1:dim(BB)[2]),side=1,las=1,line=1, family= "sans", cex=1.5)
   # }
    #dev.off();
    save(Out2,file="Out2.rdata")
    save(Out11,file="Out11.rdata")
    ##### Identify unstable lineage subpopulations ###########
    load("Sdata22.rdata")
    poorid = fidelity(Sdata22,drugs=drugs,method="simes",alpha1=0.05)
    if(!is.null(poorid)) {
        print(paste("POP",poorid," cells were likely affected by all drugs",sep=""))
    }else{
        print("No subpopulation of cells were affected by all drugs")
    }

return(list(singlecelldata=Sdata,Popmeanexprsdata=Out,Diffpopmeanexprsdata=Out2,Difflineagepopmeanexprsdata=Out11, unstablepop=poorid))
    
}

