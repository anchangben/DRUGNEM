#' A function to identify lineage subpopulations that are affected by drugs.
#'
#' This function takes as input a list of expression values as data matrices of length equal to the number of lineage subpopulations, then applies a simes global test of no change to each subpopulations using all lineage markers and returns the numeric ids for the subpopulations that show statistically significant differences.
#'
#' @param x List of expression values as data matrices of length equal to the number of lineage subpopulations stored as R object. Each matrix has rows as lineage markers and columns as cells labeled as treatment conditions including Basal.
#' @param method character string of either "Simes" (simes) or bonferoni (bonfer).Default is Simes.
#' @param alpha1 cuttoff p-value parameter for benchmarking adjusted p-values . Default is 0.05.

#' @return 	A numeric vector of subpopulation ids associated with statistically significant differences across all drugs:
#' @seealso \code{\link{comparepopintrasignalalivefidelity}} which uses this function.
#' @export
#' @examples
#' poorid= fidelity(Sdata22,method="simes")
#'

##### Comparison of intracellular signaling across survival populations ########
######## Generating sime's matrix
fidelity<-function(x,drugs,method="simes",alpha1=0.05) {
    E3=list()
    for ( i in 1:length(x)){
        
        dat2=x[[i]]
        drugname=c(unique(colnames(dat2)))
        if(all(drugs%in% drugname)&(dim(dat2)[2]>2*length(drugname))) {
        E3[[i]]=drugnem:::make.R.matrix1(dat2, wt=drugs[1], pi1 = 0.01)[[3]] ###pvalue
        } else {
         E3[[i]]=NULL
        }
    }
    if(any(unlist(lapply(E3,function(x) is.null(x))))) {
        print("Atleast one subpopulation has missing data")
        poorid=NULL
    } else  {
    E4=matrix(NA, nrow=length(E3), ncol=dim(E3[[1]])[2])
    rownames(E4)=paste("POP",1:length(E3))
    colnames(E4)=drugs1[-1]
    E5=matrix(NA, nrow=length(E3), ncol=dim(E3[[1]])[2])
    rownames(E5)=paste("POP",1:length(E3))
    colnames(E5)=drugs1[-1]
    if (method=="simes") {
        for ( j in 1:length(E3)){
            
            dat3=E3[[j]]
            E4[j,]=round(apply(dat3,2,function(x) simes(x,alpha=alpha1)),4)
            E5[j,]=round(apply(dat3,2,function(x) simes(x,alpha=alpha1,returnstat = TRUE)[2]),4)
        }
    } else {
        for ( j in 1:length(E3)){
            
            dat3=E3[[j]]
            E4[j,]=round(apply(dat3,2,function(x) bonfer(x,alpha=alpha1)),4)
            E5[j,]=round(apply(dat3,2,function(x) bonfer(x,alpha=alpha1,returnstat = TRUE)[2]),4)
        }
        
    }
    
    E6=ifelse(E4<E5,0,1)
    
    nrcolors = 100
    half = 1 + nrcolors/2
    colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
    colorpalette = colorRampPalette(colpal)(nrcolors)
    ##### Plot heatmaps ####
    BB1=as.matrix(E6)
    print(E4)
    print(E5)
    print(E6)
    pdf(file=paste(patient,"Heatmap for p-values for fidelity analysis.pdf",sep=""),width=12, height=12)
    par(mar=c(3,3,3,8))
    image(x=1:dim(BB1)[2],y=1:dim(BB1)[1],z=as.matrix(t(BB1)),col = colorpalette,xlab = paste("Treatments"),xaxt="n",yaxt="n",ylab = "Subpopulations",family="sans",cex.lab=0.8,main="P-values for fidelity analysis")
    mtext(rownames(BB1),at=c(1:dim(BB1)[1]),side=4,las=2,line=1, family= "sans", cex=0.8)
    mtext(colnames(BB1),at=c(1:dim(BB1)[2]),side=1,las=1,line=1, family= "sans", cex=0.8)
    dev.off()
    popid=which(apply(E6,1,function(x) sum(x))==0)
    if(length(popid)==0) poorid=NULL
    else poorid=popid
    }
    
    return(poorid)
}
