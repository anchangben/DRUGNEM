#'  A function for importing and density downsampling of several manually gated FCS files from an input directory.
#'
#' @param file A file path or directory path name with all fcs files of interest.
#' @param drugs Parameter vector corresponding to the names of all single drugs. Paremeter must maintain the order and names exactly as in the input folder with the FCS files with the Basal or control file name coming first.
#' @param celltype Parameter vector corresponding to the manually gated cell types. Celltype names and ordering must be exact as in the manuall gated file names. Default is NULL if there are no manually gated data.
#' @param clustids index of lineage markers of interest to be used by ccast algorithm.
#' @param transformlogic An character string that defines the type of transformation. See function read.FCS in R packgae flowCore for more information. Valid values are linearize (default) or scale.The linearize transformation applies the appropriate power transform to the data while the scale transformation scales all columns to $[0,10^decades]$. defaulting to decades=0 as in the FCS4 specification. A logical can also be used: TRUE is equal to linearize and FALSE(which is the default) corresponds to no transformation.
#' @param coln column names of markers to use of equal length as clustids. If not provided the default fcs column names will be used.
#' @param rown row names of data to use of equal length as number of cells. If not provided the default fcs file names will be used.
#' @param patient Patient or subject name. Name must be included in the file names.
#' @param subanalysis A logical to downsample original FCS file(s). Default is FALSE. Recommended if total number of cells across all FCS files > 100000.
#' @param subsamplesize Number of cells to subsample using density downsampling approach by Qiu et al. 2011
#' @return 	A FCS file with all cells pooled from all input FCS files and the marker expressions and list of 4 items:
#' \item{subsetdata }{Dataframe of pooled cell expressions for all selected lineage markers}
#' \item{allexprsdata }{Dataframe of marker expression of  all pooled cells from original FCS files}
#' \item{fcsfile}{Flowframe object of all marker expression for all combined cells from original FCS files}
#' \item{groupids}{List of row ids(cell ids) corresponding to subpopulation of cells}
#' @seealso \code{\link{drugmnemmain}} which this function wraps
#' @export
#' @examples
#' out=fcsimportmanual(file,drugs,celltype,clustids=colid,transformlogic=FALSE,coln=NULL,rown=NULL,patient,subanalysis=FALSE,subsamplesize=subsamplesize)

fcsimportmanual <- function(file,drugs,celltype,clustids=NULL,transformlogic=FALSE,coln=NULL,rown=NULL,patient,subanalysis=FALSE,subsamplesize=NULL) {

     ff=list.files(file)
     patientff<-grep(patient,ff)
    d1=list()
    d2=list()
    Dall=NULL
    dat1=NULL
    dtype=NULL
    ctype=NULL
    filelist=list()
     if(subanalysis==TRUE) {
    for ( j in 1:length(drugs)) {

        treatmentff<-grep(drugs[j],ff[patientff])
        filelist[[j]]=list()
            for ( k in 1:length(celltype)) {
            drugcellff<-grep(celltype[k],ff[treatmentff])
            filelist[[j]][[k]] <- paste(file,ff[treatmentff][drugcellff],sep="")

        d1=read.FCS(filelist[[j]][[k]],transformation=transformlogic)
        dat=exprs(d1)
        subs=sample(1:dim(dat)[1],subsamplesize,replace = FALSE)
        d2=dat[subs,]
        if(!is.null(coln)) colnames(d2)=coln
        if(!is.null(rown))  {rownames(d2) =rown
        } else {
                rownames(d2) =rep(patient,dim(d2)[1])
            }

        if (is.null(clustids)) {
        D=d2
        } else {
        D=d2[,clustids]
        }
        Dall=rbind(Dall,D)
        dat1=rbind(dat1,dat)
        dtype=c(dtype,rep(j,dim(d2)[1]))
        ctype=c(ctype,rep(k,dim(d2)[1]))

              }

    }

     } else {
     for ( j in 1:length(drugs)) {
        treatmentff<-grep(drugs[j],ff[patientff])

            filelist[[j]]=list()
            for ( k in 1:length(celltype)) {
            drugcellff<-grep(celltype[k],ff[treatmentff])

            filelist[[j]][[k]] <- paste(file,ff[treatmentff][drugcellff],sep="")

        d1=read.FCS(filelist[[j]][[k]],transformation=transformlogic)
        dat=exprs(d1)

        d2=dat
        if(!is.null(coln)) colnames(d2)=coln
        if(!is.null(rown))  {rownames(d2) =rown
        } else {
                rownames(d2) =rep(patient,dim(d2)[1])
            }

        if (is.null(clustids)) {
        D=d2
        } else {
        D=d2[,clustids]
        }
        Dall=rbind(Dall,D)
        dat1=rbind(dat1,dat)
        dtype=c(dtype,rep(j,dim(d2)[1]))
        ctype=c(ctype,rep(k,dim(d2)[1]))

              }
        }
        }


    firstout <- paste(patient,"pooled_treatment.fcs",sep="")
    in_data <- dat1
    params <- flowCore::parameters(d1)
    desc <- description(d1)
    pd <- pData(params)

    channel_number <- ncol(in_data) + 1
    channel_id <- paste("$P", channel_number, sep = "")
    channel_name <- "treatment"
    channel_range <- length(unlist(filelist))
    plist <- matrix(c(channel_name, channel_name, channel_range,
    0, channel_range - 1))
    rownames(plist) <- c("name", "desc", "range", "minRange",
    "maxRange")
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    pData(params) <- pd
    out_data <- cbind(in_data, treatment=dtype)
    out_frame <- flowFrame(out_data, params, description = desc)

    keyval <- list()
    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
    keyword(out_frame) <- keyval


    in_data2 <- out_data
    params <- flowCore::parameters(out_frame)
    desc <- description(out_frame)
    pd <- pData(params)

   channel_number <- ncol(in_data2) + 1
    channel_id <- paste("$P", channel_number, sep = "")
    channel_name <- "groups"
    channel_range <- length(unlist(filelist))
    plist <- matrix(c(channel_name, channel_name, channel_range,
    0, channel_range - 1))
    rownames(plist) <- c("name", "desc", "range", "minRange",
    "maxRange")
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    pData(params) <- pd
    out_data <- cbind(in_data2, groups= ctype)
    out_frame <- flowFrame(out_data, params, description = desc)
    keyval <- list()
    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
    keyword(out_frame) <- keyval


    write.FCS(out_frame,firstout)

   rowids=list()
   last1=dim(out_data)[2]
   alldata=out_data[,-last1]
   index1=unique(out_data[,last1])
   for (i in 1:max(index1)) {
     id=which(out_data[,last1]==index1[i])
     rowids[[i]]=id
         }
    save(rowids,file="rowids.rdata")
    save(alldata,file="Allexprsdata.rda")

    return(list(subsetdata=Dall,allexprsdata=out_data,fcsfile=out_frame,groupids=rowids))

}

