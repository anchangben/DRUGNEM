#' A function for importing and density downsampling of several FCS files from an input directory.
#'
#' @param file A file path or directory path name with all fcs files of interest.
#' @param transformlogic An character string that defines the type of transformation. See function read.FCS in R packgae flowCore for more information. Valid values are linearize (default) or scale.The linearize transformation applies the appropriate power transform to the data while the scale transformation scales all columns to $[0,10^decades]$. defaulting to decades=0 as in the FCS4 specification. A logical can also be used: TRUE is equal to linearize and FALSE(which is the default) corresponds to no transformation.
#' @param asinhp asinh cofactor parameter p for scaling the data if necessary. Default is 1.
#' @param colid index of markers of interest to be used by ccast algorithm.
#' @param coln column names of markers to use of equal length as colid. If not provided the default fcs column names will be used.
#' @param rown row names of data to use of equal length as number of cells. If not provided the default fcs file names will be used.
#' @param patient Patient or subject name. Name must be included in the file names.
#' @param subanalysis A logical to downsample original FCS file(s). Default is FALSE. Recommended if total number of cells across all FCS files > 100000.
#' @param subsamplesize Number of cells to subsample using density downsampling approach by Qiu et al. 2011
#' @return 	A FCS file with all cells pooled from all input FCS files and the marker expressions and list of 3 items:
#' \item{biomarkerdata }{Dataframe of pooled cell expressions for all selected lineage markers}
#' \item{allexprsdata }{Dataframe of marker expression of  all pooled cells from original FCS files}
#' \item{fcsfile}{Flowframe object of all marker expression for all combined cells from original FCS files} 
#' @seealso \code{\link{drugmnemmain}} which this function wraps
#' @export
#' @examples
#' Dall1=fcsimport3(file=file,transformlogic=transformlogic,asinhp=asinhp,colid=colid,coln=coln,rown=rown,patient,subanalysis=subanalysis,subsamplesize=subsamplesize)

fcsimport3 <-
function(file,transformlogic=FALSE,asinhp=asinhp, colid,coln=NULL,rown=NULL,patient,subanalysis=FALSE,subsamplesize=subsamplesize) {

        ff=list.files(file)
     patientff<-grep(patient,ff)
    filelist <- paste(file,ff[patientff],sep="")
    
    d1=list()
    d2=list()
    Dall=NULL
    dat1=NULL
    ctype=NULL
    
    for ( i in 1:length(filelist)) {
        if(subanalysis==TRUE) {
            if(is.null(asinhp)) stop("asinhp must not be NULL when downsampling.")
            d0= SPADE.addDensity.downsample(filelist[i],cols = colid, arcsinh_cofactor = asinhp,desired_samples=subsamplesize)
            dat=d0
            d2=dat[,colid]
            if(!is.null(coln)) colnames(d2)=coln
            if(!is.null(rown))  {rownames(d2) = rown
            } else {
                rownames(d2) =rep(patient,dim(d2)[1])
            }
            D=d2
            Dall=rbind(Dall,D)
            dat1=rbind(dat1,dat)
            ctype=c(ctype,rep(i,dim(d2)[1]))
            
        } else {
        d1=read.FCS(filelist[i],transformation=transformlogic)
        dat=exprs(d1)
        d2=dat[,colid]
        if(!is.null(coln)) colnames(d2)=coln
        if(!is.null(rown))  {rownames(d2) = rown
        } else { rownames(d2) =rep(file,dim(d2)[1]) }
        D=d2
        Dall=rbind(Dall,D)
        dat1=rbind(dat1,dat)
        ctype=c(ctype,rep(i,dim(d2)[1]))
        
        }
    } 
    
    
    d1=read.FCS(filelist[i],transformation=transformlogic)
    in_data <- dat1
    params <- flowCore::parameters(d1)
    desc <- description(d1)
    
    pd <- pData(params)
    typevec <- ctype
        
    firstout <- "pooled_treatment.fcs"
    
    channel_number <- ncol(in_data) + 1
    channel_id <- paste("$P", channel_number, sep = "")
    channel_name <- "treatment"
    channel_range <- length(filelist)
    plist <- matrix(c(channel_name, channel_name, channel_range,
    0, channel_range - 1))
    rownames(plist) <- c("name", "desc", "range", "minRange",
    "maxRange")
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    pData(params) <- pd
    out_data <- cbind(in_data, treatment = typevec)
    out_frame <- flowFrame(out_data, params, description = desc)
    keyval <- list()
    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
    keyword(out_frame) <- keyval
    
    write.FCS(out_frame,firstout)
   
    return(list(biomarkerdata=Dall,allexprsdata=out_data,fcsfile=out_frame))
    
}

