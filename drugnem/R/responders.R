get#' A function to filter out death or dying cells and return responder cell ids.
#'
#' @param rowids List of row ids(cell ids) for each sub population of cells stored as R object.
#' @param dataP Marker expression data matrix of pooled cells from original FCS file stored as R object.
#' @param drugs Parameter vector corresponding to the names of all single drugs.
#' @param asinhp asinh parameter p for scaling the data if necessary. Default is 1.
#' @param colid vector apoptotic marker ids of interest to be used to define responders and nonresponders
#' @param antibody Vector of marker names of length equal to length of the columns in each FCS file.
#' @param patient Patient or subject name. Name must be included exactly as in the file names.
#' @param outfile Name of working directory. Drugmenem will use this directory to create other subdirectories for filtering analysis results.
#' @param respondermarker Vector string of apoptotic markers to specified by user. Must contain atleast one element. Default is NULL

#' @return 	List of dataframes with mean and fold change expression values for each drug relative to the baseline and heatmaps of fold changes with relative percentage of cells compared to basal for all drugs. The output of this functionis a list 3 items:
#' \item{SensitiveCCASTanalysis}{ccast output object containing a list of 8 items with the last corresponding to the responder and non-responder cells derived from the input FCS fileid}
#' \item{pooledresponders}{List of cell ids of length equal to the number of cell types corresponding to responders}
#' \item{poolednonresponders}{List of cell ids of length equal to the number of cell types corresponding to responders.}
#' @seealso \code{\link{drugmnemmain}} which this function wraps
#' @export
#' @examples
#' res11=responders(rowids,dataP,drugs=drugs,asinhp=asinhp,colid=colid,antibody=antibody,patient=patient,outfile=dirc,respondermarker=respondermarker)

##### Classification of responders ##############
responders<-function(rowids,dataP,drugs=drugs,asinhp=asinhp,colid=colid,antibody=antibody,patient=patient,outfile=dirc,respondermarker=NULL){

    respondermarker1=NULL
    for ( l in 1:length(rowids)) respondermarker1[[l]]=respondermarker

    res=list()
    resids=list()
    respondids=list()
    nonrespondids=list()
    resmarker=list()

    for( j in 1:length(rowids)) {
        file1=dataP[rowids[[j]],]
        sids2=file1[,dim(dataP)[2]]
        tt2=table(sids2)
        if((length(drugs) == length(tt2))&(all(tt2>=5))) {
            ##### GET responders for all subpopulations ###########
            resmarker[[j]]=respondermarker1[[j]]
            resids[[j]]=match(resmarker[[j]],antibody)

            dirc2=paste(outfile,j,sep="/")
            dir.create(dirc2)
            setwd(dirc2)

            res[[j]]=ccast_filter6(file=dataP[rowids[[j]],],asinhp=asinhp,colid=resids[[j]],coln=resmarker[[j]],rown=NULL,k=2,patient=patient)
            respondids[[j]]=res[[j]][[8]][[1]]
            nonrespondids[[j]]=res[[j]][[8]][[2]]

        } else {

            res[[j]]=NULL
            resmarker[[j]]= NULL
            respondids[[j]]= NULL
            nonrespondids[[j]]= NULL
        }

    }
    setwd(outfile)
    ##detach("package:party", unload=TRUE)


    return(list(SensitiveCCASTanalysis=res,pooledresponders=respondids,poolednonresponders=nonrespondids))

}


