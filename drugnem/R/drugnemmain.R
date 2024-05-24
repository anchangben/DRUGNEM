#' The main function for analyzing heterogeneous drug response data across a population of single cells to optimize combination therapy based on desired intracellular effects.
#'
#' This algorithm combines Population identification, Nested effects modeling and Rank analysis to identify and rank multi-target drug combinations from single drug effects measured at the level of single cells. Drugnemmain takes as input either (1) a folder containing all drug response FCS files or (2) A data matrix of pooled drug responses for all proteins(columns) across all cells(rows) or (3) a folder with manually gated cells for drug responses from mass cytometry comprising of n drugs and m target (marker) responses on t cells taken from single or  several patient samples and returns several pdf files including both the drug target nested network and the sorted list of all possible drug combinations and their scores with the best drug regimen at the top derived from pre-specified desired effects such as down regulation effects or effects associated with a desired response of a desired marker e.g death marker or from a model with no effects. The algorithm also varies the desired effects, computes the DRUG-NEM score and chooses the model with the desired effects that produce the maximal score. 
#'
#' @param patient Patient or subject name for data. Name must be included exactly as in the FCS file names.
#' @param file A file path or directory path name with all FCS files of interest or a cxm flow cytometry expression matrix of c number of cells and m number of markers. The last 2 columns for the data matrix must be labeled  "treatment" and "groups" respectively corresponding to the index combination of drug label and cell cluster for each cell.
#' @param outputfile Name of output directory. drugnemain will set this as the working directory for storing all output files.
#' @param manualgates A logical indicating whether input folder contains manual gated data or not. Default is FALSE.
#' @param celltype Parameter vector corresponding to the manually gated cell types. Celltype names and ordering must be exact as in the manuall gated file names. Default is NULL if there are no manually gated data.
#' @param drugs Parameter vector corresponding to the names of all single drugs. Paremeter must maintain the order and names exactly as in the input folder with the FCS files with the Basal or control file name first.
#' @param drugids Integer vector of subset index of drug ids to be used for DRUGNEM analysis derived from the previous drugs parameter to optimize a subset of drugs.
#' @param antibody Vector of marker names of length equal to length of the columns in each FCS file for each channel in the FCS file.
#' @param pProteins Vector of markers corresponding to the functional proteins used for building the Nested effects models.
#' @param transformlogic An character string that defines the type of transformation. See function read.FCS in R packgae flowCore for more information. Valid values are linearize (default) or scale.The linearize transformation applies the appropriate power transform to the data while the scale transformation scales all columns to $[0,10^decades]$. defaulting to decades=0 as in the FCS4 specification. A logical can also be used: TRUE is equal to linearize and FALSE(which is the default) corresponds to no transformation.
#' @param asinhp asinh cofactor parameter p for scaling the data if necessary or NULL. Default is 1.
#' @param colid index of lineage markers of interest to be used by ccast algorithm to define cell clusters. Must be provided by user.
#' @param coln column names of lineage markers to use of equal length as colid. Must be provided by user.
#' @param rown row names of data to use of equal length as number of cells. If not provided the default FCS file names will be used.
#' @param k minimum number of expected clusters at baseline. Must be specified by user.
#' @param param string of parameter names from an FCS file. Default is NULL.
#' @param fcsdes Description of parameter names from an FCS file. Default is NULL.
#' @param respondermarker Set of apoptotic or desired endpoint markers speficied by the user. Names must correspond to names in the input FCS files.
#' @param infer Network search approach to score the best drug network. Default is "search" for exhaustive search. For large networks(>5 nodes), we recommend "nem.greedy"
#' @param ylabel y axis marker variable used to generate biaxial plots for all inner node marker variables from the optimized CCAST decision tree. Must be one of the names from the lineage marker list.
#' @param subanalysis A logical to downsample original FCS file(s). Default is FALSE. Recommended if total number of cells across all FCS files > 100000 to save time.
#' @param subsamplesize Number of cells to subsample using density downsampling approach by Qiu et al. 2011
#' @param type Parameter defining the type of data and likelihood for scoring the Nested effect models: Default "CONTmLL" called internally by the nem function in nem r-package correspond to the effect probabilities while "CONTmLLBayes" correspond to logodds that each protein is differentially expressed.
#' @param filter A logical to filter out death or dying cells or not. Default is TRUE.
#' @param weighttype parameters defining the type of desired phenotype information to use for optimizing treatment effects. Values include: (i) 'T-stat': Regulation using T-stat of intracellular signaling effects  (ii) 'FC'(default): Regulation based on Foldchange of intracellular signaling effects, and (iii) 'deathmarker': effects associated with regulation of a death effector protein or any other desired marker.
#' @param targetprior A logical to filter out outlier or noisy or alien drug targets. Default is TRUE.
#' @param Fidelity A logical to test for fidelity of pre vs post treatment characterisation of surface or lineage subpopulations using Simes statistics and identify potential unstable baseline subpopulations. Default is FALSE.
#' @param optimizenoeffect A logical to optmize drug combinations for No effect DRUG-NEM model. Useful for normal data. Default is FALSE.
#' @return 	The output of this function is a list of the following 33 items.  Item 2 corresponds to the DRUGNEM ranking using "Any effect or no prior", item 10 corresponds to the DRUGNEM ranking using prior "Up regulation effects", item 12 correspond to the DRUGNEM ranking using prior "Down regulation effects", item 14 corresponds to the DRUGNEM ranking using prior "No effects" and item 33 corresponds the best ranking after comparing priors "Up effects" and "Down effects". It also generates several pdf(Tiff) files corresponding to various diagnostic plots including plots of initial and final optimized CCAST decision  trees for all cell  lineage sub groups, Heatmaps of fold changes for all functional or intracellular  markers across all sub populations and DRUGNEM network and ranking list of all drug combinations with corresponding score.
#' \item{DRUGNEMeffectmodel}{Graphical model showing the nested relationships between drugs and their target proteins using no prior weighting of desired effect data}
#' \item{Rankingeffectmodel}{Table showing the rankings and scores for all drug combinations using no prior desired drug effects.}
#' \item{DRUGNEMupmodel}{Graphical model showing the nested relationships between drugs and their target proteins using upregulation desired effect data.}
#' \item{Rankingupmodel}{Table showing the rankings and scores for all drug combination using upregulation desired effect data.}
#' \item{DRUGNEMdownmodel}{Graphical model showing the nested relationships between drugs and their target proteins using down regulation effect data.}
#' \item{Rankingdownmodel}{Table showing the rankings and scores for all drug combination using down regulation desired effect data.}
#' \item{DRUGNEMnoeffectmodel}{Graphical model showing the nested relationships between drugs and their target proteins using no effect desired effect data.}
#' \item{Rankingnoeffectmodel}{Table showing the rankings and scores for all drug combination using no effect desired effect data.}
#' \item{DRUGNEMupmarkermodel}{Graphical model showing the nested relationships between drugs and their target proteins using only upregulation desired effect data.}
#' \item{Rankingupmarkermodel}{Table showing the rankings and scores for all drug combination using only upregulation desired effect data.}
#' \item{DRUGNEMdownmarkermodel}{Graphical model showing the nested relationships between drugs and their target proteins using only down regulation desired effect data.}
#' \item{Rankingdownmarkermodel}{Table showing the rankings and scores for all drug combination using only down regulation desired effect data.}
#' \item{DRUGNEMnoeffectmarkermodel}{Graphical model showing the nested relationships between drugs and their target proteins using only no effect desired effect data.}
#' \item{Rankingnoeffectmarkermodel}{Table showing the rankings and scores for all drug combination using only no effect desired effect data.}
#' \item{DRUGNEMbesteffectmodel}{Graphical model showing the nested relationships between drugs and their target proteins using best desired effect data.}
#' \item{Rankingbesteffectmodel}{Table showing the rankings and scores for all drug combination using best desired effect data.}
#' \item{Up_object}{DRUGNEM R object output using upregulation desired priors.}
#' \item{Down_object}{DRUGNEM R object output using down regulation desired priors.}
#' \item{Noeffect_object}{DRUGNEM R object output using no effect desired priors.}
#' \item{Upmarker_object}{DRUGNEM R object output using only upregulation desired priors.}
#' \item{Downmarker_object}{DRUGNEM R object output using only down regulation desired priors.}
#' \item{Noeffectmarker_object}{DRUGNEM R object output using only no effect desired priors.}
#' \item{Besteffect_object}{DRUGNEM R object output using the most probable desired effects.}
#' \item{Rankingindependencemodel}{Table showing the rankings and scores for all drug combinations using no prior desired effects based on independent drug effects.}
#' \item{Rankingupindependencemodel}{Table showing the rankings and scores for all drug combination using down regulation desired effects data for independence model.}
#' \item{Rankingdownindependencemodel}{Table showing the rankings and scores for all drug combination using up regulateion desired effects data for independence model.}
#' \item{Rankingnoeffectindependencemodel}{Table showing the rankings and scores for all drug combination using no effect desired effects data for independence model.}
#' \item{Rankingupmarkerindependencemodel}{Table showing the rankings and scores for all drug combination using only down desired effects data for independence model.}
#' \item{Rankingdownmarkerindependencemodel}{Table showing the rankings and scores for all drug combination using only up desired effects data for independence model.}
#' \item{Rankingnoeffectmarkerindependencemodel}{Table showing the rankings and scores for all drug combination using only no effect desired effects data for independence model.}
#' \item{Rankingbestindependencemodel}{Table showing the rankings and scores for all drug combination using best desired effects data for independence model.}
#' \item{Bestdesired_independencemodel}{list of (i) Optimal desired type of effect and (ii) table showing the rankings and scores for all drug combination using the most likely desired effects data for the independence model.}
#' \item{Bestdesired_nestedmodel}{list of (i) Optimal desired type of effect and (ii) table showing the rankings and scores for all drug combination using the most likely desired effects data for the nested model.}
#' @export
#' @examples
#' drugnemmain(patient,file=input,outputfile=getwd(), manualgates=FALSE,celltype=NULL,drugs=drugs1,drugids=c(1,2,5,7),antibody=antibody1,pProteins=pProteins1,transformlogic=FALSE,asinhp=1,colid=clustids,coln=clustmarkers,rown=NULL,k=clustersize,param=NULL,respondermarker="cPARP",infer="search",fcsdes=NULL,ylabel=ylabel,subanalysis=TRUE,subsamplesize=10000,type="CONTmLL",filter=FALSE,weighttype="downFC",targetprior=TRUE,Fidelity=FALSE,optimizenoeffect=FALSE)

drugnemmain<-function(patient,file,outputfile, manualgates=FALSE,celltype=NULL,drugs,drugids,antibody,pProteins,transformlogic=FALSE,asinhp=1,colid=NULL,coln=NULL,rown=NULL,k=NULL,param=NULL,respondermarker=NULL,infer="search",fcsdes=NULL,ylabel=ylabel,subanalysis=FALSE,subsamplesize=20000,type="CONTmLL", filter=TRUE,weighttype="FC",targetprior=TRUE,Fidelity=FALSE,optimizenoeffect=FALSE ) {

    dirc=paste(outputfile,patient,sep="/")
    dir.create(dirc)
    setwd(dirc)

    if ((!is.character(file))&(!is.matrix(file)))
        stop("Input 'filename' must be a directory path or matrix")

    if((is.matrix(file)) & (manualgates==TRUE))
       stop("Input 'filename' cannot be both a directory path and a matrix")

    if((is.matrix(file)) & (manualgates==FALSE)) {
    last1=dim(file)[2]
    last2=last1-1
    if(colnames(file)[last1]!="groups")
        stop("Last column of input file must be a vector of integers with column name 'groups'.")

    if(colnames(file)[last2]!="treatment")
        stop("Second to the last column of input file must be a vector of integers with column name 'treatment'.")
   rowids=list()
   index1=unique(file[,last1])
   for (i in 1:max(index1)) {
     id=which(file[,last1]==index1[i])
     rowids[[i]]=id
         }
   save(rowids,file="rowids.rdata")
    alldata=file[,-last1]
    rownames(alldata)=1:dim(file)[1]
    sigids=match(pProteins,antibody)
    dataP=alldata
    save(alldata,file="Allexprsdata.rda")
    }
    if(!is.matrix(file)) {
       if (manualgates==FALSE) {
    require(party)
    res1=ccast_filter5(file,transformlogic=FALSE,asinhp=asinhp,colid=colid,coln=coln,rown=NULL,k=k,param=NULL,fcsdes=NULL,ylabel=ylabel,subanalysis=subanalysis,patient=patient,subsamplesize=subsamplesize)
    load("rowids.rdata")
    load("Allexprsdata.rda")
    dataP=alldata
    sigids=match(pProteins,antibody)

    } else {
   out=fcsimportmanual(file,drugs,celltype,clustids=NULL,transformlogic=FALSE,coln=NULL,rown=NULL,patient,subanalysis=FALSE,subsamplesize=subsamplesize)
    load("rowids.rdata")
    load("Allexprsdata.rda")
    dataP=alldata
    sigids=match(pProteins,antibody)

       }
    }

if (filter == TRUE) {
 res11=responders(rowids,dataP,drugs=drugs,asinhp=asinhp,antibody=antibody,patient=patient,outfile=dirc,respondermarker=respondermarker)

    rowids2=res11[[2]]
    rowids3=res11[[3]]
    
    dataP=alldata
     require(nem)
     ##require(RBGL)
     ##require(PerformanceAnalytics)
if(Fidelity==FALSE ) {
    res2=comparepopintrasignalalive(rowids2,dataP,asinhp=asinhp,drugs=drugs,sigids=sigids,colid=colid,antibody=antibody,patient=patient)

    res22=comparepopintrasignaldeath(rowids3,dataP,asinhp=asinhp,drugs=drugs,sigids=sigids,colid=colid,antibody=antibody,patient=patient)
    
     } else {
         
    res2=comparepopintrasignalalivefidelity(rowids2,dataP,asinhp=asinhp,drugs=drugs,sigids=sigids,colid=colid,antibody=antibody,patient=patient)
         
    res22=comparepopintrasignaldeath(rowids3,dataP,asinhp=asinhp,drugs=drugs,sigids=sigids,colid=colid,antibody=antibody,patient=patient)
     }

    load("Sdata2.rdata")
    load("Out2.rdata")
    dids=drugs[drugids]

    if (!weighttype %in% c("T-stat","FC","deathmarker"))
        stop(" weighttype must be either 'T-stat','FC', 'deathmarker'.")
        
        if (weighttype=="T-stat")  {
            
    res3=fitnemupdownnoeffecttstat(Out3=Sdata2,Outp=Out2,drugs=dids,p=respondermarker[1],patient=patient,infer=infer,type=type,targetprior=targetprior)
    y1=list()
    y2=list()
    y1[[1]]=res3[[10]]
    y1[[2]]=res3[[12]]
    y1[[3]]=res3[[14]]
    y2[[1]]=res3[[28]]
    y2[[2]]=res3[[29]]
    y2[[3]]=res3[[30]]
    best_scores1=cbind(as.numeric(y1[[1]][1,1]),cbind(as.numeric(y1[[2]][1,1]),as.numeric(y1[[3]][1,1])))
    colnames(best_scores1)=c("Up","Down","No Effect")
    best_scores2=cbind(as.numeric(y2[[1]][1,1]),cbind(as.numeric(y2[[2]][1,1]),as.numeric(y2[[3]][1,1])))
    colnames(best_scores2)=c("Up","Down","No Effect")
    if (optimizenoeffect==FALSE) {
        
        if((all(best_scores1[,1:2]>0))&(all(best_scores2[,1:2]>0))) {
        ORindex1=which.max(c((best_scores1[,1]/best_scores1[,2]),(best_scores1[,2]/best_scores1[,1])))
        OR_bestmodel1=list()
        OR_bestmodel1[[1]]=best_scores1[,ORindex1]
        OR_bestmodel1[[2]]=y1[[ORindex1]]
        ORindex2=which.max(c((best_scores2[,1]/best_scores2[,2]),(best_scores2[,2]/best_scores2[,1])))
        OR_bestmodel2=list()
        OR_bestmodel2[[1]]=best_scores2[,ORindex2]
        OR_bestmodel2[[2]]=y2[[ORindex2]]
        } else {
            print("Calculation of odds ratio not possible due to NaN values")
            OR_bestmodel1=list()
            OR_bestmodel2=list()
            OR_bestmodel1[[1]]=NULL
            OR_bestmodel1[[2]]=NULL
            OR_bestmodel2[[1]]=NULL
            OR_bestmodel2[[2]]=NULL
        }
    } else {
        OR_bestmodel1=list()
        OR_bestmodel1[[1]]=best_scores1[,3]
        OR_bestmodel1[[2]]=y1[[3]]
        OR_bestmodel2=list()
        OR_bestmodel2[[1]]=best_scores2[,3]
        OR_bestmodel2[[2]]=y2[[3]]
    }
    
    }
        if (weighttype=="FC")  {

    res3=fitnemupdownnoeffectFC(Out3=Sdata2,Outp=Out2,drugs=dids,p=respondermarker[1],patient=patient,infer=infer,type=type,targetprior=targetprior)
    y1=list()
    y2=list()
    y1[[1]]=res3[[10]]
    y1[[2]]=res3[[12]]
    y1[[3]]=res3[[14]]
    y2[[1]]=res3[[28]]
    y2[[2]]=res3[[29]]
    y2[[3]]=res3[[30]]
    best_scores1=cbind(as.numeric(y1[[1]][1,1]),cbind(as.numeric(y1[[2]][1,1]),as.numeric(y1[[3]][1,1])))
    colnames(best_scores1)=c("Up","Down","No Effect")
    best_scores2=cbind(as.numeric(y2[[1]][1,1]),cbind(as.numeric(y2[[2]][1,1]),as.numeric(y2[[3]][1,1])))
    colnames(best_scores2)=c("Up","Down","No Effect")
    if (optimizenoeffect==FALSE) {
        
        if((all(best_scores1[,1:2]>0))&(all(best_scores2[,1:2]>0))) {
        ORindex1=which.max(c((best_scores1[,1]/best_scores1[,2]),(best_scores1[,2]/best_scores1[,1])))
        OR_bestmodel1=list()
        OR_bestmodel1[[1]]=best_scores1[,ORindex1]
        OR_bestmodel1[[2]]=y1[[ORindex1]]
        ORindex2=which.max(c((best_scores2[,1]/best_scores2[,2]),(best_scores2[,2]/best_scores2[,1])))
        OR_bestmodel2=list()
        OR_bestmodel2[[1]]=best_scores2[,ORindex2]
        OR_bestmodel2[[2]]=y2[[ORindex2]]
        } else {
            print("Calculation of odds ratio not possible due to NaN values")
            OR_bestmodel1=list()
            OR_bestmodel2=list()
            OR_bestmodel1[[1]]=NULL
            OR_bestmodel1[[2]]=NULL
            OR_bestmodel2[[1]]=NULL
            OR_bestmodel2[[2]]=NULL
        }
    } else {
        OR_bestmodel1=list()
        OR_bestmodel1[[1]]=best_scores1[,3]
        OR_bestmodel1[[2]]=y1[[3]]
        OR_bestmodel2=list()
        OR_bestmodel2[[1]]=best_scores2[,3]
        OR_bestmodel2[[2]]=y2[[3]]
    }

    }
        if (weighttype=="deathmarker")  {
            
    res3=fitnemupdownnoeffectdeathmarker(Out3=Sdata2,Outp=Out2,drugs=dids,p=respondermarker[1],patient=patient,infer=infer,type=type,targetprior=targetprior)
    y1=list()
    y2=list()
    y1[[1]]=res3[[10]]
    y1[[2]]=res3[[12]]
    y1[[3]]=res3[[14]]
    y2[[1]]=res3[[28]]
    y2[[2]]=res3[[29]]
    y2[[3]]=res3[[30]]
    best_scores1=cbind(as.numeric(y1[[1]][1,1]),cbind(as.numeric(y1[[2]][1,1]),as.numeric(y1[[3]][1,1])))
    colnames(best_scores1)=c("Up","Down","No Effect")
    best_scores2=cbind(as.numeric(y2[[1]][1,1]),cbind(as.numeric(y2[[2]][1,1]),as.numeric(y2[[3]][1,1])))
    colnames(best_scores2)=c("Up","Down","No Effect")
    if (optimizenoeffect==FALSE) {
        
        if((all(best_scores1[,1:2]>0))&(all(best_scores2[,1:2]>0))) {
        ORindex1=which.max(c((best_scores1[,1]/best_scores1[,2]),(best_scores1[,2]/best_scores1[,1])))
        OR_bestmodel1=list()
        OR_bestmodel1[[1]]=best_scores1[,ORindex1]
        OR_bestmodel1[[2]]=y1[[ORindex1]]
        ORindex2=which.max(c((best_scores2[,1]/best_scores2[,2]),(best_scores2[,2]/best_scores2[,1])))
        OR_bestmodel2=list()
        OR_bestmodel2[[1]]=best_scores2[,ORindex2]
        OR_bestmodel2[[2]]=y2[[ORindex2]]
        } else {
            print("Calculation of odds ratio not possible due to NaN values")
            OR_bestmodel1=list()
            OR_bestmodel2=list()
            OR_bestmodel1[[1]]=NULL
            OR_bestmodel1[[2]]=NULL
            OR_bestmodel2[[1]]=NULL
            OR_bestmodel2[[2]]=NULL
        }
    } else {
        OR_bestmodel1=list()
        OR_bestmodel1[[1]]=best_scores1[,3]
        OR_bestmodel1[[2]]=y1[[3]]
        OR_bestmodel2=list()
        OR_bestmodel2[[1]]=best_scores2[,3]
        OR_bestmodel2[[2]]=y2[[3]]
    }

    }


 } else {
    rowids2=rowids
    rowids3=NULL
   ## if (!is.matrix(file)) detach("package:party", unload=TRUE)

      dataP=alldata
     require(nem)
     ##require(RBGL)
     ##require(PerformanceAnalytics)
     if(Fidelity==FALSE ) {
         res2=comparepopintrasignalalive(rowids2,dataP,asinhp=asinhp,drugs=drugs,sigids=sigids,colid=colid,antibody=antibody,patient=patient)
         
     } else {
         
         res2=comparepopintrasignalalivefidelity(rowids2,dataP,asinhp=asinhp,drugs=drugs,sigids=sigids,colid=colid,antibody=antibody,patient=patient)
         
     }

    load("Sdata2.rdata")
    load("Out2.rdata")
    dids=drugs[drugids]

    if (!weighttype %in% c("T-stat","FC","deathmarker"))
        stop(" weighttype must be either 'T-stat','FC','deathmarker'.")
        
        if (weighttype=="T-stat")  {
            
            res3=fitnemupdownnoeffecttstat(Out3=Sdata2,Outp=Out2,drugs=dids,p=respondermarker[1],patient=patient,infer=infer,type=type,targetprior=targetprior)
            y1=list()
            y2=list()
            y1[[1]]=res3[[10]]
            y1[[2]]=res3[[12]]
            y1[[3]]=res3[[14]]
            y2[[1]]=res3[[28]]
            y2[[2]]=res3[[29]]
            y2[[3]]=res3[[30]]
            best_scores1=cbind(as.numeric(y1[[1]][1,1]),cbind(as.numeric(y1[[2]][1,1]),as.numeric(y1[[3]][1,1])))
            colnames(best_scores1)=c("Up","Down","No Effect")
            best_scores2=cbind(as.numeric(y2[[1]][1,1]),cbind(as.numeric(y2[[2]][1,1]),as.numeric(y2[[3]][1,1])))
            colnames(best_scores2)=c("Up","Down","No Effect")
            if (optimizenoeffect==FALSE) {
        
                if((all(best_scores1[,1:2]>0))&(all(best_scores2[,1:2]>0))) {
                ORindex1=which.max(c((best_scores1[,1]/best_scores1[,2]),(best_scores1[,2]/best_scores1[,1])))
                OR_bestmodel1=list()
                OR_bestmodel1[[1]]=best_scores1[,ORindex1]
                OR_bestmodel1[[2]]=y1[[ORindex1]]
                ORindex2=which.max(c((best_scores2[,1]/best_scores2[,2]),(best_scores2[,2]/best_scores2[,1])))
                OR_bestmodel2=list()
                OR_bestmodel2[[1]]=best_scores2[,ORindex2]
                OR_bestmodel2[[2]]=y2[[ORindex2]]
                } else {
                    print("Calculation of odds ratio not possible due to NaN values")
                    OR_bestmodel1=list()
                    OR_bestmodel2=list()
                    OR_bestmodel1[[1]]=NULL
                    OR_bestmodel1[[2]]=NULL
                    OR_bestmodel2[[1]]=NULL
                    OR_bestmodel2[[2]]=NULL
                }
            } else {
                OR_bestmodel1=list()
                OR_bestmodel1[[1]]=best_scores1[,3]
                OR_bestmodel1[[2]]=y1[[3]]
                OR_bestmodel2=list()
                OR_bestmodel2[[1]]=best_scores2[,3]
                OR_bestmodel2[[2]]=y2[[3]]
            }

        }
        if (weighttype=="FC")  {
            
            res3=fitnemupdownnoeffectFC(Out3=Sdata2,Outp=Out2,drugs=dids,p=respondermarker[1],patient=patient,infer=infer,type=type,targetprior=targetprior)
            y1=list()
            y2=list()
            y1[[1]]=res3[[10]]
            y1[[2]]=res3[[12]]
            y1[[3]]=res3[[14]]
            y2[[1]]=res3[[28]]
            y2[[2]]=res3[[29]]
            y2[[3]]=res3[[30]]
            best_scores1=cbind(as.numeric(y1[[1]][1,1]),cbind(as.numeric(y1[[2]][1,1]),as.numeric(y1[[3]][1,1])))
            colnames(best_scores1)=c("Up","Down","No Effect")
            best_scores2=cbind(as.numeric(y2[[1]][1,1]),cbind(as.numeric(y2[[2]][1,1]),as.numeric(y2[[3]][1,1])))
            colnames(best_scores2)=c("Up","Down","No Effect")
            if (optimizenoeffect==FALSE) {
                
                if((all(best_scores1[,1:2]>0))&(all(best_scores2[,1:2]>0))) {
                ORindex1=which.max(c((best_scores1[,1]/best_scores1[,2]),(best_scores1[,2]/best_scores1[,1])))
                OR_bestmodel1=list()
                OR_bestmodel1[[1]]=best_scores1[,ORindex1]
                OR_bestmodel1[[2]]=y1[[ORindex1]]
                ORindex2=which.max(c((best_scores2[,1]/best_scores2[,2]),(best_scores2[,2]/best_scores2[,1])))
                OR_bestmodel2=list()
                OR_bestmodel2[[1]]=best_scores2[,ORindex2]
                OR_bestmodel2[[2]]=y2[[ORindex2]]
                } else {
                    print("Calculation of odds ratio not possible due to NaN values")
                    OR_bestmodel1=list()
                    OR_bestmodel2=list()
                    OR_bestmodel1[[1]]=NULL
                    OR_bestmodel1[[2]]=NULL
                    OR_bestmodel2[[1]]=NULL
                    OR_bestmodel2[[2]]=NULL
                }
            } else {
                OR_bestmodel1=list()
                OR_bestmodel1[[1]]=best_scores1[,3]
                OR_bestmodel1[[2]]=y1[[3]]
                OR_bestmodel2=list()
                OR_bestmodel2[[1]]=best_scores2[,3]
                OR_bestmodel2[[2]]=y2[[3]]
            }

        }
        if (weighttype=="deathmarker")  {
            
            res3=fitnemupdownnoeffectdeathmarker(Out3=Sdata2,Outp=Out2,drugs=dids,p=respondermarker[1],patient=patient,infer=infer,type=type,targetprior=targetprior)
            y1=list()
            y2=list()
            y1[[1]]=res3[[10]]
            y1[[2]]=res3[[12]]
            y1[[3]]=res3[[14]]
            y2[[1]]=res3[[28]]
            y2[[2]]=res3[[29]]
            y2[[3]]=res3[[30]]
            best_scores1=cbind(as.numeric(y1[[1]][1,1]),cbind(as.numeric(y1[[2]][1,1]),as.numeric(y1[[3]][1,1])))
            colnames(best_scores1)=c("Up","Down","No Effect")
            best_scores2=cbind(as.numeric(y2[[1]][1,1]),cbind(as.numeric(y2[[2]][1,1]),as.numeric(y2[[3]][1,1])))
            colnames(best_scores2)=c("Up","Down","No Effect")
            if (optimizenoeffect==FALSE) {
                
                if((all(best_scores1[,1:2]>0))&(all(best_scores2[,1:2]>0))) {
                ORindex1=which.max(c((best_scores1[,1]/best_scores1[,2]),(best_scores1[,2]/best_scores1[,1])))
                OR_bestmodel1=list()
                OR_bestmodel1[[1]]=best_scores1[,ORindex1]
                OR_bestmodel1[[2]]=y1[[ORindex1]]
                ORindex2=which.max(c((best_scores2[,1]/best_scores2[,2]),(best_scores2[,2]/best_scores2[,1])))
                OR_bestmodel2=list()
                OR_bestmodel2[[1]]=best_scores2[,ORindex2]
                OR_bestmodel2[[2]]=y2[[ORindex2]]
                } else {
                    print("Calculation of odds ratio not possible due to NaN values")
                    OR_bestmodel1=list()
                    OR_bestmodel2=list()
                    OR_bestmodel1[[1]]=NULL
                    OR_bestmodel1[[2]]=NULL
                    OR_bestmodel2[[1]]=NULL
                    OR_bestmodel2[[2]]=NULL
                }
            } else {
                OR_bestmodel1=list()
                OR_bestmodel1[[1]]=best_scores1[,3]
                OR_bestmodel1[[2]]=y1[[3]]
                OR_bestmodel2=list()
                OR_bestmodel2[[1]]=best_scores2[,3]
                OR_bestmodel2[[2]]=y2[[3]]
            }

        }
        

    }

return(list(DRUGNEMeffectmodel=res3[[1]],Rankingeffectmodel=res3[[2]],DRUGNEMupmodel=res3[[3]],Rankingupmodel=res3[[4]],DRUGNEMdownmodel=res3[[5]],Rankingdownmodel=res3[[6]],DRUGNEMnoeffectmodel=res3[[7]],Rankingnoeffectmodel=res3[[8]],DRUGNEMupmarkermodel=res3[[9]],Rankingupmarkermodel=res3[[10]],DRUGNEMdownmarkermodel=res3[[11]],Rankingdownmarkermodel=res3[[12]],DRUGNEMnoeffectmarkermodel=res3[[13]],Rankingnoeffectmarkermodel=res3[[14]],DRUGNEMbesteffectmodel=res3[[15]],Rankingbesteffectmodel=res3[[16]],Up_object=res3[[17]],Down_object=res3[[18]], Noeffect_object=res3[[19]],Upmarker_object=res3[[20]],Downmarker_object=res3[[21]], Noeffectmarker_object=res3[[22]],Besteffect_object=res3[[23]],Rankingindependencemodel=res3[[24]],Rankingupindependencemodel=res3[[25]],Rankingdownindependencemodel=res3[[26]],Rankingnoeffectindependencemodel=res3[[27]],Rankingupmarkerindependencemodel=res3[[28]],Rankingdownmarkerindependencemodel=res3[[29]],Rankingnoeffectmarkerindependencemodel=res3[[30]],Rankingbestindependencemodel=res3[[31]],Bestdesired_independencemodel=OR_bestmodel2,Bestdesired_nestedmodel=OR_bestmodel1))

}


