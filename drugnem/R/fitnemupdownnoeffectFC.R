#' A function to predict the DRUGNEM network, scores drug combinations and plot diagnostic heatmaps for combination therapy using weighttype=="FC" option.
#'
#' @param Out3 R object corresponding to a list of dataframes for each subpopulation with protein expressions for each cell across functional proteins
#' @param Outp R object corresponding to a list of dataframes for each subpopulation with Fold change expressions for proteins(rows) and drugs (columns) relative to the basal treatment.
#' @param drugs Integer vector of subset index of drug ids to optimize drug combinations from all available drugs.
#' @param p Apoptotic marker from the respondermarker parameter. Default is NULL.
#' @param patient Patient or subject name. Name must be included exactly as in the file names.
#' @param infer Network search approach to score the best drug network. Default is "search" for exhaustive search. For large networks use "nem.greedy"
#' @param type Parameter defining the type of data and likelihood for scoring the Nested effect models: Default "CONTmLL" called internally by the nem function in nem package correspond to the effect probabilities while "CONTmLLBayes" correspond to logodds that each protein is differentially expressed.
#' @param targetprior A logical to filter out outlier or noisy drug targets. Default is FALSE.
#' @return 	The output of this function is a list 31 items:
#' \item{DRUGNEMeffectmodel}{Graphical model showing the nested relationships between drugs and their target proteins using an effect as a desired effect}
#' \item{Rankingeffectmodel}{Table showing the rankings and scores for all drug combinations using using an effect as a desired effect.}
#' \item{DRUGNEMupmodel}{Graphical model showing the nested relationships between drugs and their target proteins using effects associated with positive t-stat as desired effects.}
#' \item{Rankingupmodel}{Table showing the rankings and scores for all drug combination using   effects associated with positive t-stat as desired effects.}
#' \item{DRUGNEMdownmodel}{Graphical model showing the nested relationships between drugs and their target proteins using  effects associated with negative t-stat as desired effects.}
#' \item{Rankingdownmodel}{Table showing the rankings and scores for all drug combination using  effects associated with negative t-stat as as desired effects.}
#  \item{DRUGNEMnoeffectmodel}{Graphical model showing the nested relationships between drugs and their target proteins using  markers associated with no effects as desired effects.}
#' \item{Rankingnoeffectmodel}{Table showing the rankings and scores for all drug combination using  markers associated with no effects as desired effects.}
#' \item{DRUGNEMupmarkermodel}{Graphical model showing the nested relationships between drugs and the markers that will produce an upregulation desired effect.}
#' \item{Rankingupmarkermodel}{Table showing the rankings and scores for all drug combination using the markers that will produce an upregulation desired effect.}
#' \item{DRUGNEMdownmarkermodel}{Graphical model showing the nested relationships between drugs and the markers that will produce a downregulation desired effect.}
#' \item{Rankingdownmarkermodel}{Table showing the rankings and scores for all drug combination using the markers that will produce a downregulation desired effect.}
#'  \item{DRUGNEMnoeffectmarkermodel}{Graphical model showing the nested relationships between drugs and the markers that will produce a no effect as desired effects.}
#' \item{Rankingnoeffectmarkermodel}{Table showing the rankings and scores for all drug combination using the markers that will produce a no effect as desired effects.}
#' \item{DRUGNEMbesteffectmodel}{Graphical model showing the nested relationships between drugs and their target proteins using most likely observed regulation as desired effects.}
#' \item{Rankingbesteffectmodel}{Table showing the rankings and scores for all drug combination using most likely observed regulation as desired effects.}
#' \item{Up_object}{DRUGNEM R object output using positive t-stat as desired priors.}
#' \item{Down_object}{DRUGNEM R object output using negative t-stat desired priors.}
#' \item{Noeffect_object}{DRUGNEM R object output using no effect desired priors.}
#' \item{Upmarker_object}{DRUGNEM R object output using only the markers that will produce an upregulation desired effect.}
#' \item{Downmarker_object}{DRUGNEM R object output using only the markers that will produce a downregulation desired effect.}
#' \item{Noeffectmarker_object}{DRUGNEM R object output using only the markers that will produce a no effect as desired effects.}
#' \item{Besteffect_object}{DRUGNEM R object output using the most likely observed regulation as desired effects.}
#' \item{Rankingindependencemodel}{Table showing the rankings and scores for all drug combinations using effects associated with positive t-stat as desired effects based on independent drug effects.}
#' \item{Rankingupindependencemodel}{Table showing the rankings and scores for all drug combination using  effects associated with positive t-stat as desired effects for independence model.}
#' \item{Rankingdownindependencemodel}{Table showing the rankings and scores for all drug combination using  effects associated with negative t-stat as desired effects for independence model.}
#' \item{Rankingnoeffectindependencemodel}{Table showing the rankings and scores for all drug combination using no effects as effect desired effects data for independence model.}
#' \item{Rankingupmarkerindependencemodel}{Table showing the rankings and scores for all drug combination using only the markers that will produce an upregulation desired effect for independence model.}
#' \item{Rankingdownmarkerindependencemodel}{Table showing the rankings and scores for all drug combination using only the markers that will produce an downregulation desired effect for independence model.}
#' \item{Rankingnoeffectmarkerindependencemodel}{Table showing the rankings and scores for all drug combination using the markers that will produce an no effect as desired effect for independence model.}
#' \item{Rankingbestindependencemodel}{Table showing the rankings and scores for all drug combination using the most likely observed regulation as desired effects for independence model.}
#' @seealso \code{\link{drugnemmain}} which this function wraps
#' @export
#' @examples
#' res3=fitmnemdown(Out3=Sdata2,Outp=Out2,drugs=dids,p=NULL,patient=patient,infer=infer,type=type)

###### Fit,score and Plot DRUGNEM with heatmap on single pdf##########

fitnemupdownnoeffectFC<-function(Out3,Outp,drugs,p=NULL,patient,infer,type="CONTmLL",targetprior){
    Deffect=list()
    Deffect2=list()
    Deffect3=list()
    Deffect4=list()
    Deffect4b=list()
    Deffect5=list()
    Deffect5b=list()
    Deffect6=list()
    Deffect6b=list()
    Deffect7=list()
    Deffect_label=list()
    Deffect_labelp=list()
    pop.a.c=list()
    
    celltype <- paste("Pop",1:length(Out3))
    nrcolors = 100
    half = 1 + nrcolors/2
    colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9,"OrRd")[1:9])
    colorpalette = colorRampPalette(colpal)(nrcolors)
    
    
    #####Fit subpopulation NEM ########
    cdata1=NULL
    cdata2=NULL
    cdata2b=NULL
    cdata3=NULL
    cdata3b=NULL
    cdata6=NULL
    cdata6b=NULL
    cdata7=NULL
   
    Esize=matrix(0,nrow=length(celltype),ncol=(length(drugs)-1))
     colnames(Esize)=drugs[-1]
    rownames(Esize)=celltype
    for (t in 1:length(celltype)) {
        
        fcp19=Out3[[t]]
        deathid=which(rownames(fcp19)==p)
         drugname=c(unique(colnames(fcp19)))
         did=match(drugs[-1],drugname[-1])
         if(is.null(fcp19)) {
             Deffect[[t]]=NULL
             Deffect2[[t]]=NULL
             Deffect3[[t]]=NULL
             Deffect4[[t]]=NULL
             Deffect4b[[t]]=NULL
             Deffect5[[t]]=NULL
             Deffect5b[[t]]=NULL
             Deffect6[[t]]=NULL
             Deffect6[[t]]=NULL
             Deffect6b[[t]]=NULL
             Deffect7[[t]]=NULL
             Deffect_label[[t]]=NULL
             Deffect_labelp[[t]]=NULL
             rr3=NULL
             rr4=NULL
             rr4b=NULL
             rr5=NULL
             rr5b=NULL
             rr6=NULL
             rr6b=NULL
             rr7=NULL
            
         } else {
        if(all(drugs%in% drugname)&(dim(fcp19)[2]>2*length(drugname))) {
            dat1=fcp19[deathid,,drop=FALSE]
            Deffect[[t]]=make.R.matrix2(dat1, wt=drugs[1], pi1 = 0.01)
            Deffect3[[t]]=Deffect[[t]][[2]][,did,drop=FALSE]
            Esize[t,]=Deffect3[[t]]
            dat2=fcp19[-deathid,]
            hdat=make.R.matrix1(dat2, wt=drugs[1], pi1 = 0.01)[[1]]
            hdat2=ifelse(Outp[-deathid,-1,t]>0,hdat,rnorm(1,-10,0.05))
            hdat2b=ifelse(Outp[-deathid,-1,t]<0,hdat,rnorm(1,-10,0.05))
            hdat2c=ifelse(Outp[-deathid,-1,t]==0,-hdat,rnorm(1,-10,0.05))
            hdat2p=exp(hdat2)/(1+exp(hdat2))
            hdat2bp=exp(hdat2b)/(1+exp(hdat2b))
            hdat2cp=1/(1+exp(hdat))
            hdat3=pmax(hdat2c,pmax(hdat2,hdat2b))
            hdat3p=pmax(hdat2cp,pmax(hdat2p,hdat2bp))
            hdat3_labelp=hdat3p
            for (a in 1:dim(hdat3p)[1]) {
                for ( b in 1:dim(hdat3p)[2]) {
                    hdat3_labelp[a,b]=match(hdat3p[a,b],c(hdat2p[a,b],hdat2bp[a,b],hdat2cp[a,b]))
                }
            }
            hdat3_label=hdat3
            for (a in 1:dim(hdat3)[1]) {
                for ( b in 1:dim(hdat3)[2]) {
                    hdat3_label[a,b]=match(hdat3[a,b],c(hdat2[a,b],hdat2b[a,b],hdat2c[a,b]))
                }
            }
            Deffect_label[[t]]=hdat3_label[,did]
            Deffect_labelp[[t]]=hdat3_labelp[,did]
            if(!(type %in% c("mLL", "CONTmLL"))) {
                dat33=hdat2[,drugs[-1]]
                dat33b=hdat2b[,drugs[-1]]
                dat33c=hdat2c[,drugs[-1]]
                dat33d=hdat3[,drugs[-1]]
                Deffect2[[t]]=hdat[,drugs[-1]]
                Deffect4[[t]]=dat33
                Deffect5[[t]]=dat33b
                Deffect6[[t]]=dat33c
                Deffect7[[t]]=dat33d
                if(targetprior==TRUE) {
                    rr3=nem(Deffect2[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect2[[t]])),type="CONTmLLBayes"), verbose=FALSE)
                    rr4=nem(Deffect4[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect4[[t]])),type="CONTmLLBayes"), verbose=FALSE)
                    #Deffect4b[[t]]=ifelse(Deffect_label[[t]]==1,Deffect2[[t]],rnorm(1,-10,0.05))
                    Deffect4b[[t]]=ifelse(Deffect_label[[t]]==1,Deffect7[[t]],rnorm(1,-10,0.05))
                    rr4b=nem(Deffect4b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect4b[[t]])),type="CONTmLLBayes"), verbose=FALSE)
                    rr5=nem(Deffect5[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect5[[t]])),type="CONTmLLBayes"), verbose=FALSE)
                    #Deffect5b[[t]]=ifelse(Deffect_label[[t]]==2,Deffect2[[t]],rnorm(1,-10,0.05))
                    Deffect5b[[t]]=ifelse(Deffect_label[[t]]==2,Deffect7[[t]],rnorm(1,-10,0.05))
                    rr5b=nem(Deffect5b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect5b[[t]])),type="CONTmLLBayes"), verbose=FALSE)
                    rr6=nem(Deffect6[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect6[[t]])),type="CONTmLLBayes"), verbose=FALSE)
                    #Deffect6b[[t]]=ifelse(Deffect_label[[t]]==3,Deffect2[[t]],rnorm(1,-10,0.05))
                    Deffect6b[[t]]=ifelse(Deffect_label[[t]]==3,Deffect7[[t]],rnorm(1,-10,0.05))
                    rr6b=nem(Deffect6b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect6b[[t]])),type="CONTmLLBayes"), verbose=FALSE)
                    rr7=nem(Deffect7[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect7[[t]])),type="CONTmLLBayes"), verbose=FALSE)
                    
                } else {
                    rr3=nem(Deffect2[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect2[[t]])),type="CONTmLLBayes",delta=0), verbose=FALSE)
                    rr4=nem(Deffect4[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect4[[t]])),type="CONTmLLBayes",delta=0), verbose=FALSE)
                    #Deffect4b[[t]]=ifelse(Deffect_label[[t]]==1,Deffect2[[t]],rnorm(1,-10,0.05))
                    Deffect4b[[t]]=ifelse(Deffect_label[[t]]==1,Deffect7[[t]],rnorm(1,-10,0.05))
                    rr4b=nem(Deffect4b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect4b[[t]])),type="CONTmLLBayes",delta=0), verbose=FALSE)
                    rr5=nem(Deffect5[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect5[[t]])),type="CONTmLLBayes",delta=0), verbose=FALSE)
                    #Deffect5b[[t]]=ifelse(Deffect_label[[t]]==2,Deffect2[[t]],rnorm(1,-10,0.05))
                    Deffect5b[[t]]=ifelse(Deffect_label[[t]]==2,Deffect7[[t]],rnorm(1,-10,0.05))
                    rr5b=nem(Deffect5b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect5b[[t]])),type="CONTmLLBayes",delta=0), verbose=FALSE)
                    rr6=nem(Deffect6[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect6[[t]])),type="CONTmLLBayes",delta=0), verbose=FALSE)
                    #Deffect6b[[t]]=ifelse(Deffect_label[[t]]==3,Deffect2[[t]],rnorm(1,-10,0.05))
                    Deffect6b[[t]]=ifelse(Deffect_label[[t]]==3,Deffect7[[t]],rnorm(1,-10,0.05))
                    rr6b=nem(Deffect6b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect6b[[t]])),type="CONTmLLBayes",delta=0), verbose=FALSE)
                    rr7=nem(Deffect7[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect7[[t]])),type="CONTmLLBayes",delta=0), verbose=FALSE)
                }
            } else {
                dat33=hdat2p[,drugs[-1]]
                dat33b=hdat2bp[,drugs[-1]]
                dat33c=hdat2cp[,drugs[-1]]
                dat33d=hdat3p[,drugs[-1]]
                dat31=round(dat33,6)
                dat31b=round(dat33b,6)
                dat31c=round(dat33c,6)
                dat31d=round(dat33d,6)
                numdata=length(dat31)
                numdatab=length(dat31b)
                numdatac=length(dat31c)
                numdatad=length(dat31d)
                
                ss=abs(rnorm(numdata,0,0.005))
                ssb=abs(rnorm(numdatab,0,0.005))
                ssc=abs(rnorm(numdatac,0,0.005))
                ssd=abs(rnorm(numdatad,0,0.005))
                dat32=ifelse(is.na(dat31),1-ss,dat31)
                dat32b=ifelse(is.na(dat31b),1-ss,dat31b)
                dat32c=ifelse(is.na(dat31c),1-ss,dat31c)
                dat32d=ifelse(is.na(dat31d),1-ss,dat31d)
                dat34=ifelse(dat32>=0.99,dat32-ss,dat32)[,drugs[-1]]
                dat34b=ifelse(dat32b>=0.99,dat32b-ssb,dat32b)[,drugs[-1]]
                dat34c=ifelse(dat32c>=0.99,dat32c-ssb,dat32c)[,drugs[-1]]
                dat34d=ifelse(dat32d>=0.99,dat32d-ssb,dat32d)[,drugs[-1]]
                Deffect4[[t]]=dat34
                Deffect5[[t]]=dat34b
                Deffect6[[t]]=dat34c
                Deffect7[[t]]=dat34d
                
                ddat31=round(exp(hdat[,drugs[-1]])/(1+exp(hdat[,drugs[-1]])),6)
                numdata=length(ddat31)
                sss=abs(rnorm(numdata,0,0.005))
                
                ddat32=ifelse(is.na(ddat31),1-sss,ddat31)
                ddat33=ifelse(ddat32>=0.99,ddat32-ss,ddat32)[,drugs[-1]]
                Deffect2[[t]]=ddat33
                
                Deffect2[[t]]=ifelse(Deffect2[[t]]==0,abs(rnorm(1,0,0.005)),Deffect2[[t]])
                Deffect4[[t]]=ifelse(Deffect4[[t]]==0,abs(rnorm(1,0,0.005)),Deffect4[[t]])
                Deffect5[[t]]=ifelse(Deffect5[[t]]==0,abs(rnorm(1,0,0.005)),Deffect5[[t]])
                Deffect6[[t]]=ifelse(Deffect6[[t]]==0,abs(rnorm(1,0,0.005)),Deffect6[[t]])
                Deffect7[[t]]=ifelse(Deffect7[[t]]==0,abs(rnorm(1,0,0.005)),Deffect7[[t]])
                
                if(targetprior==TRUE) {
                    rr3=nem(Deffect2[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect2[[t]])),type=type), verbose=FALSE)
                    rr4=nem(Deffect4[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect4[[t]])),type=type), verbose=FALSE)
                    #Deffect4b[[t]]=ifelse(Deffect_label[[t]]==1,Deffect2[[t]],abs(rnorm(1,0,0.05)))
                    Deffect4b[[t]]=ifelse(Deffect_labelp[[t]]==1,Deffect7[[t]],abs(rnorm(1,0,0.005)))
                    rr4b=nem(Deffect4b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect4b[[t]])),type=type), verbose=FALSE)
                    rr5=nem(Deffect5[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect5[[t]])),type=type), verbose=FALSE)
                    #Deffect5b[[t]]=ifelse(Deffect_labelp[[t]]==2,Deffect2[[t]],abs(rnorm(1,0,0.05)))
                    Deffect5b[[t]]=ifelse(Deffect_labelp[[t]]==2,Deffect7[[t]],abs(rnorm(1,0,0.005)))
                    rr5b=nem(Deffect5b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect5b[[t]])),type=type), verbose=FALSE)
                    rr6=nem(Deffect6[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect6[[t]])),type=type), verbose=FALSE)
                    #Deffect6b[[t]]=ifelse(Deffect_labelp[[t]]==3,Deffect2[[t]],abs(rnorm(1,0,0.05)))
                    Deffect6b[[t]]=ifelse(Deffect_labelp[[t]]==3,Deffect7[[t]],abs(rnorm(1,0,0.005)))
                    rr6b=nem(Deffect6b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect6b[[t]])),type=type), verbose=FALSE)
                    rr7=nem(Deffect7[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect7[[t]])),type=type), verbose=FALSE)
                } else {
                    rr3=nem(Deffect2[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect2[[t]])),type=type,delta=0), verbose=FALSE)
                    rr4=nem(Deffect4[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect4[[t]])),type=type,delta=0), verbose=FALSE)
                    #Deffect4b[[t]]=ifelse(Deffect_labelp[[t]]==1,Deffect2[[t]],abs(rnorm(1,0,0.05)))
                    Deffect4b[[t]]=ifelse(Deffect_labelp[[t]]==1,Deffect7[[t]],abs(rnorm(1,0,0.005)))
                    rr4b=nem(Deffect4b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect4b[[t]])),type=type,delta=0), verbose=FALSE)
                    rr5=nem(Deffect5[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect5[[t]])),type=type,delta=0), verbose=FALSE)
                    #Deffect5b[[t]]=ifelse(Deffect_labelp[[t]]==2,Deffect2[[t]],abs(rnorm(1,0,0.05)))
                    Deffect5b[[t]]=ifelse(Deffect_labelp[[t]]==2,Deffect7[[t]],abs(rnorm(1,0,0.005)))
                    rr5b=nem(Deffect5b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect5b[[t]])),type=type,delta=0), verbose=FALSE)
                    rr6=nem(Deffect6[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect6[[t]])),type=type,delta=0), verbose=FALSE)
                    #Deffect6b[[t]]=ifelse(Deffect_labelp[[t]]==3,Deffect2[[t]],abs(rnorm(1,0,0.05)))
                    Deffect6b[[t]]=ifelse(Deffect_labelp[[t]]==3,Deffect7[[t]],abs(rnorm(1,0,0.005)))
                    rr6b=nem(Deffect6b[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect6b[[t]])),type=type,delta=0), verbose=FALSE)
                    rr7=nem(Deffect7[[t]], inference = infer , control = set.default.parameters(unique(colnames(Deffect7[[t]])),type=type,delta=0), verbose=FALSE)
                }
                
            }
            
            ####nemselected = rr3$selected
            if (infer=="search") {
                nemselected1= unique(unlist(rr3$mappos[[which.max(rr3$mLL)]]))
                nemselected2= unique(unlist(rr4$mappos[[which.max(rr4$mLL)]]))
                nemselected2b= unique(unlist(rr4b$mappos[[which.max(rr4b$mLL)]]))
                nemselected3= unique(unlist(rr5$mappos[[which.max(rr5$mLL)]]))
                nemselected3b= unique(unlist(rr5b$mappos[[which.max(rr5b$mLL)]]))
                nemselected4= unique(unlist(rr6$mappos[[which.max(rr6$mLL)]]))
                nemselected4b= unique(unlist(rr6b$mappos[[which.max(rr6b$mLL)]]))
                nemselected5= unique(unlist(rr7$mappos[[which.max(rr7$mLL)]]))
            } else {
               nemselected1= unique(unlist(rr3$mappos))
               nemselected2= unique(unlist(rr4$mappos))
               nemselected2b= unique(unlist(rr4b$mappos))
               nemselected3= unique(unlist(rr5$mappos))
               nemselected3b= unique(unlist(rr5b$mappos))
               nemselected4= unique(unlist(rr6$mappos))
               nemselected4b= unique(unlist(rr6b$mappos))
               nemselected5= unique(unlist(rr7$mappos))
            }
            
            if(length(nemselected1)!=0) {
                TM0=Deffect2[[t]][nemselected1,,drop = FALSE]
                pic2=TM0
            } else {
                TM0=Deffect2[[t]]
                pic2=TM0
            }
            if(length(nemselected2)!=0) {
                TM2=Deffect4[[t]][nemselected2,,drop = FALSE]
                pic3=TM2
            } else {
                TM2=Deffect4[[t]]
                pic3=TM2
            }
            if(length(nemselected2b)!=0) {
                TM2b=Deffect4b[[t]][nemselected2b,,drop = FALSE]
                pic3b=TM2b
            } else {
                TM2b=Deffect4b[[t]]
                pic3b=TM2b
            }
            if(length(nemselected3)!=0) {
                TM3=Deffect5[[t]][nemselected3,,drop = FALSE]
                pic4=TM3
            } else {
                TM3=Deffect5[[t]]
                pic4=TM3
            }
            if(length(nemselected3b)!=0) {
                TM3b=Deffect5b[[t]][nemselected3b,,drop = FALSE]
                pic4b=TM3b
            } else {
                TM3b=Deffect5b[[t]]
                pic4b=TM3b
            }
            if(length(nemselected4)!=0) {
                TM6=Deffect6[[t]][nemselected4,,drop = FALSE]
                pic6=TM6
            } else {
                TM6=Deffect6[[t]]
                pic6=TM6
            }
            if(length(nemselected4b)!=0) {
                TM6b=Deffect6b[[t]][nemselected4b,,drop = FALSE]
                pic6b=TM6b
            } else {
                TM6b=Deffect6b[[t]]
                pic6b=TM6b
            }
            if(length(nemselected5)!=0) {
                TM7=Deffect7[[t]][nemselected5,,drop = FALSE]
                pic7=TM7
            } else {
                TM7=Deffect7[[t]]
                pic7=TM7
            }
            
            rownames(TM0)=paste(rownames(TM0),celltype[t])
            cdata1=rbind(cdata1,TM0)
            rownames(TM2)=paste(rownames(TM2),celltype[t])
            cdata2=rbind(cdata2,TM2)
            rownames(TM2b)=paste(rownames(TM2b),celltype[t])
            cdata2b=rbind(cdata2b,TM2b)
            rownames(TM3)=paste(rownames(TM3),celltype[t])
            cdata3=rbind(cdata3,TM3)
            rownames(TM3b)=paste(rownames(TM3b),celltype[t])
            cdata3b=rbind(cdata3b,TM3b)
            rownames(TM6)=paste(rownames(TM6),celltype[t])
            cdata6=rbind(cdata6,TM6)
            rownames(TM6b)=paste(rownames(TM6b),celltype[t])
            cdata6b=rbind(cdata6b,TM6b)
            rownames(TM7)=paste(rownames(TM7),celltype[t])
            cdata7=rbind(cdata7,TM7)
            
            TM=graph2adj(rr3$graph)
            TMs=graph2adj(rr4$graph)
            TMs2=graph2adj(rr5$graph)
            TMs3=graph2adj(rr6$graph)
            TMsb=graph2adj(rr4b$graph)
            TMs2b=graph2adj(rr5b$graph)
            TMs3b=graph2adj(rr6b$graph)
            TMs4=graph2adj(rr7$graph)
            
        } else {
            Deffect[[t]]=NULL
            Deffect2[[t]]=NULL
            Deffect3[[t]]=NULL
            Deffect4[[t]]=NULL
            Deffect5[[t]]=NULL
            Deffect6[[t]]=NULL
            Deffect4b[[t]]=NULL
            Deffect5b[[t]]=NULL
            Deffect6b[[t]]=NULL
            Deffect7[[t]]=NULL
            Deffect_label[[t]]=NULL
            rr3=NULL
            rr4=NULL
            rr5=NULL
            rr6=NULL
            rr4b=NULL
            rr5b=NULL
            rr6b=NULL
            rr7=NULL
        }
        
    }
    }
    save(Deffect,file="Deffect.rdata") ### Limma output for each subpopulation
    save(Deffect2,file="Deffect2.rdata") ### Protein effect for each drug for each subpopulation
    save(Deffect3,file="Deffect3.rdata")  ### One-sided p-value for death marker(s) for each subpopulation
    save(Deffect4,file="Deffect4.rdata") ### up effect matrix1 for each subpopulation
    save(Deffect5,file="Deffect5.rdata")### down effect matrix2 for each subpopulation
    save(Deffect6,file="Deffect6.rdata")### noeffect effect matrix2 for each subpopulation
    save(Deffect4b,file="Deffect4b.rdata") ### up effect matrix1b for each subpopulation
    save(Deffect5b,file="Deffect5b.rdata")### down effect matrix2b for each subpopulation
    save(Deffect6b,file="Deffect6b.rdata")### noeffect effect matrix2b for each subpopulation
    save(Deffect_label,file="Deffect_label.rdata")### Regulatory effect matrix for each subpopulation
     save(Deffect_labelp,file="Deffect_labelp.rdata")### Regulatory prob effect matrix for each subpopulation
    save(Deffect7,file="Deffect7.rdata")### max effect matrix2 for each subpopulation
    print("Saving Esize matrix")
    save(Esize,file="Esize.rdata") ### subpopulation death associated effect matrix
    ##### Fit population NEM ########
    if(!(type %in% c("mLL", "CONTmLL"))) {
        rownames(cdata1)=paste(rownames(cdata1),1:dim(cdata1)[1])
        rownames(cdata2)=paste(rownames(cdata2),1:dim(cdata2)[1])
        rownames(cdata3)=paste(rownames(cdata3),1:dim(cdata3)[1])
        rownames(cdata6)=paste(rownames(cdata6),1:dim(cdata6)[1])
        rownames(cdata2b)=paste(rownames(cdata2b),1:dim(cdata2b)[1])
        rownames(cdata3b)=paste(rownames(cdata3b),1:dim(cdata3b)[1])
        rownames(cdata6b)=paste(rownames(cdata6b),1:dim(cdata6b)[1])
        rownames(cdata7)=paste(rownames(cdata7),1:dim(cdata7)[1])
        #rownames(cdata6)=paste(rownames(cdata6),1:dim(cdata6)[1])
        
        if(targetprior==TRUE) {
            rr=nem(cdata1, inference = infer , control = set.default.parameters(unique(colnames(cdata1)),type="CONTmLLBayes"), verbose=FALSE)
            cdata11=round(exp(cdata1)/(1+exp(cdata1)),6)
            
            rr1=nem(cdata2, inference = infer , control = set.default.parameters(unique(colnames(cdata2)),type="CONTmLLBayes"), verbose=FALSE)
            cdata22=round(exp(cdata2)/(1+exp(cdata2)),6)
            
            rr2=nem(cdata3, inference = infer , control = set.default.parameters(unique(colnames(cdata3)),type="CONTmLLBayes"), verbose=FALSE)
            cdata33=round(exp(cdata3)/(1+exp(cdata3)),6)
            
            rr66=nem(cdata6, inference = infer , control = set.default.parameters(unique(colnames(cdata6)),type="CONTmLLBayes"), verbose=FALSE)
            cdata66a=round(exp(cdata6)/(1+exp(cdata6)),6)
            
            rr1b=nem(cdata2b, inference = infer , control = set.default.parameters(unique(colnames(cdata2b)),type="CONTmLLBayes"), verbose=FALSE)
            cdata22b=round(exp(cdata2b)/(1+exp(cdata2b)),6)
            
            rr2b=nem(cdata3b, inference = infer , control = set.default.parameters(unique(colnames(cdata3b)),type="CONTmLLBayes"), verbose=FALSE)
            cdata33b=round(exp(cdata3b)/(1+exp(cdata3b)),6)
            
            rr66b=nem(cdata6b, inference = infer , control = set.default.parameters(unique(colnames(cdata6b)),type="CONTmLLBayes"), verbose=FALSE)
            cdata66b=round(exp(cdata6b)/(1+exp(cdata6b)),6)
            
            rr77=nem(cdata7, inference = infer , control = set.default.parameters(unique(colnames(cdata7)),type="CONTmLLBayes"), verbose=FALSE)
            cdata77a=round(exp(cdata7)/(1+exp(cdata7)),6)
        } else {
            rr=nem(cdata1, inference = infer , control = set.default.parameters(unique(colnames(cdata1)),type="CONTmLLBayes",delta=0), verbose=FALSE)
            cdata11=round(exp(cdata1)/(1+exp(cdata1)),6)
            
            rr1=nem(cdata2, inference = infer , control = set.default.parameters(unique(colnames(cdata2)),type="CONTmLLBayes",delta=0), verbose=FALSE)
            cdata22=round(exp(cdata2)/(1+exp(cdata2)),6)
            
            rr2=nem(cdata3, inference = infer , control = set.default.parameters(unique(colnames(cdata3)),type="CONTmLLBayes",delta=0), verbose=FALSE)
            cdata33=round(exp(cdata3)/(1+exp(cdata3)),6)
            
            rr66=nem(cdata6, inference = infer , control = set.default.parameters(unique(colnames(cdata6)),type="CONTmLLBayes",delta=0), verbose=FALSE)
            cdata66a=round(exp(cdata6)/(1+exp(cdata6)),6)
            
            rr1b=nem(cdata2b, inference = infer , control = set.default.parameters(unique(colnames(cdata2b)),type="CONTmLLBayes",delta=0), verbose=FALSE)
            cdata22b=round(exp(cdata2b)/(1+exp(cdata2b)),6)
            
            rr2b=nem(cdata3b, inference = infer , control = set.default.parameters(unique(colnames(cdata3b)),type="CONTmLLBayes",delta=0), verbose=FALSE)
            cdata33b=round(exp(cdata3b)/(1+exp(cdata3b)),6)
            
            rr66b=nem(cdata6b, inference = infer , control = set.default.parameters(unique(colnames(cdata6b)),type="CONTmLLBayes",delta=0), verbose=FALSE)
            cdata66b=round(exp(cdata6b)/(1+exp(cdata6b)),6)
            
            rr77=nem(cdata7, inference = infer , control = set.default.parameters(unique(colnames(cdata7)),type="CONTmLLBayes",delta=0), verbose=FALSE)
            cdata77a=round(exp(cdata7)/(1+exp(cdata7)),6)
        }
        ###### Plot MNEM ##########
        TM=graph2adj(rr$graph)
        TMs=graph2adj(rr1$graph)
        TMs2=graph2adj(rr2$graph)
        TMs3=graph2adj(rr66$graph)
        TMsb=graph2adj(rr1b$graph)
        TMs2b=graph2adj(rr2b$graph)
        TMs3b=graph2adj(rr66b$graph)
        TMs4=graph2adj(rr77$graph)
        if((sum(TM!=0)>1)&(edgemode(rr$graph)=="directed")) {
            plot.nem2(rr,filename=paste(patient,"DRUGNEM effect network and Heatmap for all states",infer,type,".pdf",sep=""),main="DRUGNEM effect Network",PDF=TRUE,what = "graph",D=cdata1,draw.lines = FALSE)
           # if(any(rr$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
         #   else plot.nem(rr,filename=paste(patient,"Maxlikelihood for top networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"DRUGNEM effect network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr$graph,main="DRUGNEM effect network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr$mappos[[which.max(rr$mLL)]]))
            } else {
                nemselected= unique(unlist(rr$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic2=cdata1[nemselected,,drop = FALSE]
                image(x=1:dim(pic2)[2],y=1:dim(pic2)[1],z=as.matrix(t(pic2)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","-","pooled"," drug effects", sep=""))
                mtext(rownames(pic2),at=c(1:dim(pic2)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic2),at=c(1:dim(pic2)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        
        if((sum(TMs!=0)>1)&(edgemode(rr1$graph)=="directed")) {
            plot.nem2(rr1,filename=paste(patient,"Up DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Up DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata2,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
          #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Up DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr1$graph,main="Up DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr1$mappos[[which.max(rr1$mLL)]]))
            } else {
                nemselected= unique(unlist(rr1$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic3=cdata22[nemselected,,drop = FALSE]
                image(x=1:dim(pic3)[2],y=1:dim(pic3)[1],z=as.matrix(t(pic3)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","-","pooled"," drug effects", sep=""))
                mtext(rownames(pic3),at=c(1:dim(pic3)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic3),at=c(1:dim(pic3)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        
        if((sum(TMs2!=0)>1)&(edgemode(rr2$graph)=="directed")) {
            plot.nem2(rr2,filename=paste(patient,"Down DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Down DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata3,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Down DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr2$graph,main="Down DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr2$mappos[[which.max(rr2$mLL)]]))
            } else {
                nemselected= unique(unlist(rr2$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic4=cdata33[nemselected,,drop = FALSE]
                image(x=1:dim(pic4)[2],y=1:dim(pic4)[1],z=as.matrix(t(pic4)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","-","pooled"," drug effects", sep=""))
                mtext(rownames(pic4),at=c(1:dim(pic4)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic4),at=c(1:dim(pic4)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        if((sum(TMs3!=0)>1)&(edgemode(rr66$graph)=="directed")) {
            plot.nem2(rr66,filename=paste(patient,"No effect DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="No effect DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata66a,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"No effect DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr66$graph,main="No effect DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr66$mappos[[which.max(rr66$mLL)]]))
            } else {
                nemselected= unique(unlist(rr66$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic6=cdata66a[nemselected,,drop = FALSE]
                image(x=1:dim(pic6)[2],y=1:dim(pic6)[1],z=as.matrix(t(pic6)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","pooled"," drug No effects", sep=""))
                mtext(rownames(pic6),at=c(1:dim(pic6)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic6),at=c(1:dim(pic6)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        
        if((sum(TMsb!=0)>1)&(edgemode(rr1b$graph)=="directed")) {
            plot.nem2(rr1b,filename=paste(patient,"Up markers DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Up markers DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata2b,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Up markers DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr1b$graph,main="Up markers DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr1b$mappos[[which.max(rr1b$mLL)]]))
            } else {
                nemselected= unique(unlist(rr1b$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic3b=cdata22b[nemselected,,drop = FALSE]
                image(x=1:dim(pic3b)[2],y=1:dim(pic3b)[1],z=as.matrix(t(pic3b)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","-","pooled"," drug effects", sep=""))
                mtext(rownames(pic3b),at=c(1:dim(pic3b)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic3b),at=c(1:dim(pic3b)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        
        if((sum(TMs2b!=0)>1)&(edgemode(rr2b$graph)=="directed")) {
            plot.nem2(rr2b,filename=paste(patient,"Down markers DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Down DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata3b,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Down markers DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr2b$graph,main="Down markers DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr2b$mappos[[which.max(rr2b$mLL)]]))
            } else {
                nemselected= unique(unlist(rr2b$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic4b=cdata33b[nemselected,,drop = FALSE]
                image(x=1:dim(pic4b)[2],y=1:dim(pic4b)[1],z=as.matrix(t(pic4b)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","-","pooled"," drug effects", sep=""))
                mtext(rownames(pic4b),at=c(1:dim(pic4b)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic4b),at=c(1:dim(pic4b)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        if((sum(TMs3b!=0)>1)&(edgemode(rr66b$graph)=="directed")) {
            plot.nem2(rr66b,filename=paste(patient,"No effect markers DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="No effect DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata6b,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"No effect markers DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr66b$graph,main="No effect markers DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr66b$mappos[[which.max(rr66b$mLL)]]))
            } else {
                nemselected= unique(unlist(rr66b$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic6b=cdata66b[nemselected,,drop = FALSE]
                image(x=1:dim(pic6b)[2],y=1:dim(pic6b)[1],z=as.matrix(t(pic6b)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","pooled"," drug No effects", sep=""))
                mtext(rownames(pic6b),at=c(1:dim(pic6b)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic6b),at=c(1:dim(pic6b)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }

        if((sum(TMs4!=0)>1)&(edgemode(rr77$graph)=="directed")) {
            plot.nem2(rr77,filename=paste(patient,"Best effect DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Best effect DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata77a,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Best effect DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr77$graph,main="Best effect DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr77$mappos[[which.max(rr77$mLL)]]))
            } else {
                nemselected= unique(unlist(rr77$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic7=cdata77a[nemselected,,drop = FALSE]
                image(x=1:dim(pic7)[2],y=1:dim(pic7)[1],z=as.matrix(t(pic7)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","pooled"," drug best effects", sep=""))
                mtext(rownames(pic7),at=c(1:dim(pic7)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic7),at=c(1:dim(pic7)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }


        if (infer=="search") {
            mappos1=rr$mappos[[which.max(rr$mLL)]]
            mappos2=rr1$mappos[[which.max(rr1$mLL)]]
            mappos3=rr2$mappos[[which.max(rr2$mLL)]]
            mappos4=rr66$mappos[[which.max(rr66$mLL)]]
            mappos2b=rr1b$mappos[[which.max(rr1b$mLL)]]
            mappos3b=rr2b$mappos[[which.max(rr2b$mLL)]]
            mappos4b=rr66b$mappos[[which.max(rr66b$mLL)]]
            mappos5=rr77$mappos[[which.max(rr77$mLL)]]
            pos1=rr$pos[[which.max(rr$mLL)]]
            pos2=rr1$pos[[which.max(rr1$mLL)]]
            pos3=rr2$pos[[which.max(rr2$mLL)]]
            pos4=rr66$pos[[which.max(rr66$mLL)]]
            pos2b=rr1b$pos[[which.max(rr1b$mLL)]]
            pos3b=rr2b$pos[[which.max(rr2b$mLL)]]
            pos4b=rr66b$pos[[which.max(rr66b$mLL)]]
            pos5=rr77$pos[[which.max(rr77$mLL)]]
        } else {
            mappos1=rr$mappos
            mappos2=rr1$mappos
            mappos3=rr2$mappos
            mappos4=rr66$mappos
            mappos2b=rr1b$mappos
            mappos3b=rr2b$mappos
            mappos4b=rr66b$mappos
            mappos5=rr77$mappos
            pos1=rr$pos
            pos2=rr1$pos
            pos3=rr2$pos
            pos4=rr66$pos
            pos2b=rr1b$pos
            pos3b=rr2b$pos
            pos4b=rr66b$pos
            pos5=rr77$pos
        }
        
        if(length(unlist(strsplit(colnames(cdata11),"+",fixed=TRUE)))>length(drugs[-1])){
            
            Mat01=NULL
            Mat01s=NULL
            Mat01t=NULL
            Mat01u=NULL
            Mat01sb=NULL
            Mat01tb=NULL
            Mat01ub=NULL
            Mat01v=NULL
            
            Mat1=NULL
            Mat1s=NULL
            Mat1t=NULL
            Mat1u=NULL
            Mat1sb=NULL
            Mat1tb=NULL
            Mat1ub=NULL
            Mat1v=NULL
            #Mat2=NULL
            #Mat2s=NULL
            #Mat2t=NULL
            
        } else {
            out22=bestdrugtarget(G=rr$graph,mappos=mappos1,data=pos1)
            out22s=bestdrugtarget(G=rr1$graph,mappos=mappos2,data=pos2)
            out22t=bestdrugtarget(G=rr2$graph,mappos=mappos3,data=pos3)
            out22u=bestdrugtarget(G=rr66$graph,mappos=mappos4,data=pos4)
            out22sb=bestdrugtarget(G=rr1b$graph,mappos=mappos2b,data=pos2b)
            out22tb=bestdrugtarget(G=rr2b$graph,mappos=mappos3b,data=pos3b)
            out22ub=bestdrugtarget(G=rr66b$graph,mappos=mappos4b,data=pos4b)
            out22v=bestdrugtarget(G=rr77$graph,mappos=mappos5,data=pos5)
            out3=scoredrugcombbliss(cdata11)
            out3s=scoredrugcombbliss(cdata22)
            out3t=scoredrugcombbliss(cdata33)
            out3u=scoredrugcombbliss(cdata66a)
            out3sb=scoredrugcombbliss(cdata22b)
            out3tb=scoredrugcombbliss(cdata33b)
            out3ub=scoredrugcombbliss(cdata66b)
            out3v=scoredrugcombbliss(cdata77a)
            #out4=scoredrugcombbliss(cdata11)
            #out4s=scoredrugcombbliss(cdata22)
            
            
            pdf(file=paste(drugs,"DRUG-COMBINATION scores all states.pdf", sep=""), width=13, height=13)


            Mat01=matrix(out22[[2]],nrow=length(out22[[2]]))
            rownames(Mat01)=names(out22[[2]])
            textplot(round(Mat01,4),cex.main=3,cex=2.5,family="sans")
            title("DRUGNEM effect scores")
            
            Mat01s=matrix(out22s[[2]],nrow=length(out22s[[2]]))
            rownames(Mat01s)=names(out22s[[2]])
            textplot(round(Mat01s,4),cex.main=3,cex=2.5,family="sans")
            title("Up DRUGNEM scores")
            
            Mat01t=matrix(out22t[[2]],nrow=length(out22t[[2]]))
            rownames(Mat01t)=names(out22t[[2]])
            textplot(round(Mat01t,4),cex.main=3,cex=2.5,family="sans")
            title("Down DRUGNEM scores")
            
            Mat01u=matrix(out22u[[2]],nrow=length(out22u[[2]]))
            rownames(Mat01u)=names(out22u[[2]])
            textplot(round(Mat01u,4),cex.main=3,cex=2.5,family="sans")
            title("No effect DRUGNEM scores")
            
            Mat01sb=matrix(out22sb[[2]],nrow=length(out22sb[[2]]))
            rownames(Mat01sb)=names(out22sb[[2]])
            textplot(round(Mat01sb,4),cex.main=3,cex=2.5,family="sans")
            title("Up markers DRUGNEM scores")
            
            Mat01tb=matrix(out22tb[[2]],nrow=length(out22tb[[2]]))
            rownames(Mat01tb)=names(out22tb[[2]])
            textplot(round(Mat01tb,4),cex.main=3,cex=2.5,family="sans")
            title("Down markers DRUGNEM scores")
            
            Mat01ub=matrix(out22ub[[2]],nrow=length(out22ub[[2]]))
            rownames(Mat01ub)=names(out22ub[[2]])
            textplot(round(Mat01ub,4),cex.main=3,cex=2.5,family="sans")
            title("No effect markers DRUGNEM scores")
            
            Mat01v=matrix(out22v[[2]],nrow=length(out22v[[2]]))
            rownames(Mat01v)=names(out22v[[2]])
            textplot(round(Mat01v,4),cex.main=3,cex=2.5,family="sans")
            title("Best effect DRUGNEM scores")
            
            Mat1=matrix(out3[[2]],nrow=length(out3[[2]]))
            rownames(Mat1)=names(out3[[2]])
            textplot(round(Mat1,4),cex.main=3,cex=2.5,family="sans")
            title("Drug Independence effect Scores")
            
            Mat1s=matrix(out3s[[2]],nrow=length(out3s[[2]]))
            rownames(Mat1s)=names(out3s[[2]])
            textplot(round(Mat1s,4),cex.main=3,cex=2.5,family="sans")
            title("Up Drug Independence Scores")
            
            Mat1t=matrix(out3t[[2]],nrow=length(out3t[[2]]))
            rownames(Mat1t)=names(out3t[[2]])
            textplot(round(Mat1t,4),cex.main=3,cex=2.5,family="sans")
            title("Down Drug Independence Scores")
            
            Mat1u=matrix(out3u[[2]],nrow=length(out3u[[2]]))
            rownames(Mat1u)=names(out3u[[2]])
            textplot(round(Mat1u,4),cex.main=3,cex=2.5,family="sans")
            title("No effect Drug Independence Scores")
            
            Mat1sb=matrix(out3sb[[2]],nrow=length(out3sb[[2]]))
            rownames(Mat1sb)=names(out3sb[[2]])
            textplot(round(Mat1sb,4),cex.main=3,cex=2.5,family="sans")
            title("Up markers Drug Independence Scores")
            
            Mat1tb=matrix(out3tb[[2]],nrow=length(out3tb[[2]]))
            rownames(Mat1tb)=names(out3tb[[2]])
            textplot(round(Mat1tb,4),cex.main=3,cex=2.5,family="sans")
            title("Down markers Drug Independence Scores")
            
            Mat1ub=matrix(out3ub[[2]],nrow=length(out3ub[[2]]))
            rownames(Mat1ub)=names(out3ub[[2]])
            textplot(round(Mat1ub,4),cex.main=3,cex=2.5,family="sans")
            title("No effect markers Drug Independence Scores")
            
            Mat1v=matrix(out3v[[2]],nrow=length(out3v[[2]]))
            rownames(Mat1v)=names(out3v[[2]])
            textplot(round(Mat1v,4),cex.main=3,cex=2.5,family="sans")
            title("Best effect Drug Independence Scores")
            
            #Mat2=matrix(out4[[2]],nrow=length(out4[[2]]))
            #rownames(Mat2)=names(out4[[2]])
            #textplot(round(Mat2,4),cex.main=3,cex=2.5,family="sans")
            #title("Additivescore")
            
            #Mat2s=matrix(out4s[[2]],nrow=length(out4s[[2]]))
            #rownames(Mat2s)=names(out4s[[2]])
            #textplot(round(Mat2s,4),cex.main=3,cex=2.5,family="sans")
            #title("Weighted Additivescore")
            
            dev.off();
            
        }
        
    } else {
        
        rownames(cdata1)=paste(rownames(cdata1),1:dim(cdata1)[1])
        rownames(cdata2)=paste(rownames(cdata2),1:dim(cdata2)[1])
        rownames(cdata3)=paste(rownames(cdata3),1:dim(cdata3)[1])
        rownames(cdata6)=paste(rownames(cdata6),1:dim(cdata6)[1])
        rownames(cdata2b)=paste(rownames(cdata2b),1:dim(cdata2b)[1])
        rownames(cdata3b)=paste(rownames(cdata3b),1:dim(cdata3b)[1])
        rownames(cdata6b)=paste(rownames(cdata6b),1:dim(cdata6b)[1])
        rownames(cdata7)=paste(rownames(cdata7),1:dim(cdata7)[1])
        
        if(targetprior==TRUE) {
            rr=nem(cdata1, inference = infer , control = set.default.parameters(unique(colnames(cdata1)),type=type), verbose=FALSE)
            
            rr1=nem(cdata2, inference = infer , control = set.default.parameters(unique(colnames(cdata2)),type=type), verbose=FALSE)
            
            rr2=nem(cdata3, inference = infer , control = set.default.parameters(unique(colnames(cdata3)),type=type), verbose=FALSE)
            
            rr66=nem(cdata6, inference = infer , control = set.default.parameters(unique(colnames(cdata6)),type=type), verbose=FALSE)
            
            rr1b=nem(cdata2b, inference = infer , control = set.default.parameters(unique(colnames(cdata2b)),type=type), verbose=FALSE)
            
            rr2b=nem(cdata3b, inference = infer , control = set.default.parameters(unique(colnames(cdata3b)),type=type), verbose=FALSE)
            
            rr66b=nem(cdata6b, inference = infer , control = set.default.parameters(unique(colnames(cdata6b)),type=type), verbose=FALSE)
            
            rr77=nem(cdata7, inference = infer , control = set.default.parameters(unique(colnames(cdata7)),type=type), verbose=FALSE)
            
        } else {
            rr=nem(cdata1, inference = infer , control = set.default.parameters(unique(colnames(cdata1)),type=type,delta=0), verbose=FALSE)
            
            rr1=nem(cdata2, inference = infer , control = set.default.parameters(unique(colnames(cdata2)),type=type,delta=0), verbose=FALSE)
          
            rr2=nem(cdata3, inference = infer , control = set.default.parameters(unique(colnames(cdata3)),type=type,delta=0), verbose=FALSE)
            
            rr66=nem(cdata6, inference = infer , control = set.default.parameters(unique(colnames(cdata6)),type=type,delta=0), verbose=FALSE)
            
            rr1b=nem(cdata2b, inference = infer , control = set.default.parameters(unique(colnames(cdata2b)),type=type,delta=0), verbose=FALSE)
            
            rr2b=nem(cdata3b, inference = infer , control = set.default.parameters(unique(colnames(cdata3b)),type=type,delta=0), verbose=FALSE)
            
            rr66b=nem(cdata6b, inference = infer , control = set.default.parameters(unique(colnames(cdata6b)),type=type,delta=0), verbose=FALSE)

            
            rr77=nem(cdata7, inference = infer , control = set.default.parameters(unique(colnames(cdata7)),type=type,delta=0), verbose=FALSE)
            
        }
        
        TM=graph2adj(rr$graph)
        TMs=graph2adj(rr1$graph)
        TMs2=graph2adj(rr2$graph)
        TMs3=graph2adj(rr66$graph)
        TMsb=graph2adj(rr1b$graph)
        TMs2b=graph2adj(rr2b$graph)
        TMs3b=graph2adj(rr66b$graph)
        TMs4=graph2adj(rr77$graph)


        if((sum(TM!=0)>1)&(edgemode(rr$graph)=="directed")) {
            ###### Plot MNEM ##########
            plot.nem2(rr,filename=paste(patient,"DRUGNEM effect network and Heatmap for all states",infer,type,".pdf",sep=""),main="DRUGNEM effect Network",PDF=TRUE,what = "graph",D=cdata1,draw.lines = FALSE)
           # if(any(rr$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
           # else plot.nem(rr,filename=paste(patient,"Maxlikelihood for top networks for all states.pdf",sep=""),main="DRUGNEM # Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
            
        } else {
            pdf(file=paste(patient,"DRUGNEM effect network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr$graph,main="DRUGNEM effect network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr$mappos[[which.max(rr$mLL)]]))
            } else {
                nemselected= unique(unlist(rr$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic2=cdata1[nemselected,,drop = FALSE]
                image(x=1:dim(pic2)[2],y=1:dim(pic2)[1],z=as.matrix(t(pic2)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap of drug effects for all states",infer,type, sep=""))
                mtext(rownames(pic2),at=c(1:dim(pic2)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic2),at=c(1:dim(pic2)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
            
        }
        
        if((sum(TMs!=0)>1)&(edgemode(rr1$graph)=="directed")) {
            ###### Plot MNEM ##########
            plot.nem2(rr1,filename=paste(patient,"Up DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Up DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata2,draw.lines = FALSE)
           # if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
          #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood-Survival for top networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
            
        } else {
            pdf(file=paste(patient,"Up DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr1$graph,main="Up DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr1$mappos[[which.max(rr1$mLL)]]))
            } else {
                nemselected= unique(unlist(rr1$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic3=cdata2[nemselected,,drop = FALSE]
                image(x=1:dim(pic3)[2],y=1:dim(pic3)[1],z=as.matrix(t(pic3)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap of drug effects for all states",infer,type, sep=""))
                mtext(rownames(pic3),at=c(1:dim(pic3)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic3),at=c(1:dim(pic3)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
            
        }
        
        if((sum(TMs2!=0)>1)&(edgemode(rr2$graph)=="directed")) {
            plot.nem2(rr2,filename=paste(patient,"Down DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Down DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata3,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Down DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr2$graph,main="Down DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr2$mappos[[which.max(rr2$mLL)]]))
            } else {
                nemselected= unique(unlist(rr2$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic4=cdata3[nemselected,,drop = FALSE]
                image(x=1:dim(pic4)[2],y=1:dim(pic4)[1],z=as.matrix(t(pic4)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","-","pooled"," drug effects", sep=""))
                mtext(rownames(pic4),at=c(1:dim(pic4)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic4),at=c(1:dim(pic4)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        if((sum(TMs3!=0)>1)&(edgemode(rr66$graph)=="directed")) {
            plot.nem2(rr66,filename=paste(patient,"No effect DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="No effect DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata6,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"No effect DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr66$graph,main="No effect DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr66$mappos[[which.max(rr66$mLL)]]))
            } else {
                nemselected= unique(unlist(rr66$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic6=cdata6[nemselected,,drop = FALSE]
                image(x=1:dim(pic6)[2],y=1:dim(pic6)[1],z=as.matrix(t(pic6)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","pooled"," drug No effects", sep=""))
                mtext(rownames(pic6),at=c(1:dim(pic6)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic6),at=c(1:dim(pic6)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        
        if((sum(TMsb!=0)>1)&(edgemode(rr1b$graph)=="directed")) {
            plot.nem2(rr1b,filename=paste(patient,"Up markers DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Up markers DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata2b,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Up markers DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr1b$graph,main="Up markers DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr1b$mappos[[which.max(rr1b$mLL)]]))
            } else {
                nemselected= unique(unlist(rr1b$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic3b=cdata2b[nemselected,,drop = FALSE]
                image(x=1:dim(pic3b)[2],y=1:dim(pic3b)[1],z=as.matrix(t(pic3b)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","-","pooled"," drug effects", sep=""))
                mtext(rownames(pic3b),at=c(1:dim(pic3b)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic3b),at=c(1:dim(pic3b)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        
        if((sum(TMs2b!=0)>1)&(edgemode(rr2b$graph)=="directed")) {
            plot.nem2(rr2b,filename=paste(patient,"Down markers DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Down markers DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata3b,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Down markers DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr2b$graph,main="Down markers DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr2b$mappos[[which.max(rr2b$mLL)]]))
            } else {
                nemselected= unique(unlist(rr2b$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic4b=cdata3b[nemselected,,drop = FALSE]
                image(x=1:dim(pic4b)[2],y=1:dim(pic4b)[1],z=as.matrix(t(pic4b)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","-","pooled"," drug effects", sep=""))
                mtext(rownames(pic4b),at=c(1:dim(pic4b)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic4b),at=c(1:dim(pic4b)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }
        if((sum(TMs3b!=0)>1)&(edgemode(rr66b$graph)=="directed")) {
            plot.nem2(rr66b,filename=paste(patient,"No effect markers DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="No effect DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata6b,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"No effect markers DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr66b$graph,main="No effect markers DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr66b$mappos[[which.max(rr66b$mLL)]]))
            } else {
                nemselected= unique(unlist(rr66b$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic6b=cdata6b[nemselected,,drop = FALSE]
                image(x=1:dim(pic6b)[2],y=1:dim(pic6b)[1],z=as.matrix(t(pic6b)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","pooled"," drug No effects", sep=""))
                mtext(rownames(pic6b),at=c(1:dim(pic6b)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic6b),at=c(1:dim(pic6b)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }

        if((sum(TMs4!=0)>1)&(edgemode(rr77$graph)=="directed")) {
            plot.nem2(rr77,filename=paste(patient,"Best effect DRUGNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="Best effect DRUGNEM Network",PDF=TRUE,what = "graph",D=cdata7,draw.lines = FALSE)
            #if(any(rr1$mLL==Inf)) print("No Maximum likelihood profile plotted due to Inf values")
            #  else plot.nem(rr1,filename=paste(patient,"Maxlikelihood for top survival networks for all states.pdf",sep=""),main="DRUGNEM Network",PDF=TRUE,what = "mLL",D=NULL,draw.lines = TRUE)
            
        } else {
            pdf(file=paste(patient,"Best effect DRUGNEM network and Heatmap of pooled cells  with",infer,type,patient, ".pdf", sep=""),width=15, height=15)
            
            plot(rr77$graph,main="Best effect DRUGNEM network")
            par(mar=c(5,2,5,5))
            if (infer=="search") {
                nemselected= unique(unlist(rr77$mappos[[which.max(rr77$mLL)]]))
            } else {
                nemselected= unique(unlist(rr77$mappos))
            }
            
            if(length(nemselected)!=0) {
                pic7=cdata7[nemselected,,drop = FALSE]
                image(x=1:dim(pic7)[2],y=1:dim(pic7)[1],z=as.matrix(t(pic7)),col = colorpalette,xlab = "Drugs",xaxt="n",yaxt="n",ylab = "target-proteins",family="sans",cex.lab=1.4,main=paste("Heatmap","pooled"," drug best effects", sep=""))
                mtext(rownames(pic7),at=c(1:dim(pic7)[1]),side=4,las=2,line=1, family= "sans", cex=0.5)
                mtext(colnames(pic7),at=c(1:dim(pic7)[2]),side=1,las=1,line=1, family= "sans", cex=1)
                dev.off();
            }
        }



        if (infer=="search") {
            mappos1=rr$mappos[[which.max(rr$mLL)]]
            mappos2=rr1$mappos[[which.max(rr1$mLL)]]
            mappos3=rr2$mappos[[which.max(rr2$mLL)]]
            mappos4=rr66$mappos[[which.max(rr66$mLL)]]
            mappos2b=rr1b$mappos[[which.max(rr1b$mLL)]]
            mappos3b=rr2b$mappos[[which.max(rr2b$mLL)]]
            mappos4b=rr66b$mappos[[which.max(rr66b$mLL)]]
            mappos5=rr77$mappos[[which.max(rr77$mLL)]]
            pos1=rr$pos[[which.max(rr$mLL)]]
            pos2=rr1$pos[[which.max(rr1$mLL)]]
            pos3=rr2$pos[[which.max(rr2$mLL)]]
            pos4=rr66$pos[[which.max(rr66$mLL)]]
            pos2b=rr1b$pos[[which.max(rr1b$mLL)]]
            pos3b=rr2b$pos[[which.max(rr2b$mLL)]]
            pos4b=rr66b$pos[[which.max(rr66b$mLL)]]
            pos5=rr77$pos[[which.max(rr77$mLL)]]
        } else {
            mappos1=rr$mappos
            mappos2=rr1$mappos
            mappos3=rr2$mappos
            mappos4=rr66$mappos
            mappos2b=rr1b$mappos
            mappos3b=rr2b$mappos
            mappos4b=rr66b$mappos
            mappos5=rr77$mappos
            pos1=rr$pos
            pos2=rr1$pos
            pos3=rr2$pos
            pos4=rr66$pos
            pos2b=rr1b$pos
            pos3b=rr2b$pos
            pos4b=rr66b$pos
            pos5=rr77$pos
        }
      
      
        
        if(length(unlist(strsplit(colnames(cdata1),"+",fixed=TRUE)))>length(drugs[-1])){
            
            Mat01=NULL
            Mat01s=NULL
            Mat01t=NULL
            Mat01u=NULL
            Mat01sb=NULL
            Mat01tb=NULL
            Mat01ub=NULL
            Mat01v=NULL
            Mat1=NULL
            Mat1s=NULL
            Mat1t=NULL
            Mat1u=NULL
            Mat1sb=NULL
            Mat1tb=NULL
            Mat1ub=NULL
            Mat1v=NULL
            #Mat2=NULL
            #Mat2s=NULL
            #Mat2t=NULL
            
        } else {
            
            out22=bestdrugtarget(G=rr$graph,mappos=mappos1,data=pos1)
            out22s=bestdrugtarget(G=rr1$graph,mappos=mappos2,data=pos2)
            out22t=bestdrugtarget(G=rr2$graph,mappos=mappos3,data=pos3)
            out22u=bestdrugtarget(G=rr66$graph,mappos=mappos4,data=pos4)
            out22sb=bestdrugtarget(G=rr1b$graph,mappos=mappos2b,data=pos2b)
            out22tb=bestdrugtarget(G=rr2b$graph,mappos=mappos3b,data=pos3b)
            out22ub=bestdrugtarget(G=rr66b$graph,mappos=mappos4b,data=pos4b)
            out22v=bestdrugtarget(G=rr77$graph,mappos=mappos5,data=pos5)

            out3=scoredrugcombbliss(cdata1)
            out3s=scoredrugcombbliss(cdata2)
            out3t=scoredrugcombbliss(cdata3)
            out3u=scoredrugcombbliss(cdata6)
            out3sb=scoredrugcombbliss(cdata2b)
            out3tb=scoredrugcombbliss(cdata3b)
            out3ub=scoredrugcombbliss(cdata6b)
            out3v=scoredrugcombbliss(cdata7)
            #out4=scoredrugcombbliss(cdata1)
            #out4s=scoredrugcombbliss(cdata2)
            
            pdf(file=paste(drugs,"DRUG-COMBINATION scores all states.pdf", sep=""), width=13, height=13)
            
            Mat01=matrix(out22[[2]],nrow=length(out22[[2]]))
            rownames(Mat01)=names(out22[[2]])
            textplot(round(Mat01,4),cex.main=3,cex=2.5,family="sans")
            title("DRUGNEM effect scores")
            
            Mat01s=matrix(out22s[[2]],nrow=length(out22s[[2]]))
            rownames(Mat01s)=names(out22s[[2]])
            textplot(round(Mat01s,4),cex.main=3,cex=2.5,family="sans")
            title("Up DRUGNEM scores")
            
            Mat01t=matrix(out22t[[2]],nrow=length(out22t[[2]]))
            rownames(Mat01t)=names(out22t[[2]])
            textplot(round(Mat01t,4),cex.main=3,cex=2.5,family="sans")
            title("Down DRUGNEM scores")
            
            Mat01u=matrix(out22u[[2]],nrow=length(out22u[[2]]))
            rownames(Mat01u)=names(out22u[[2]])
            textplot(round(Mat01u,4),cex.main=3,cex=2.5,family="sans")
            title("No effect DRUGNEM scores")
            
            Mat01sb=matrix(out22sb[[2]],nrow=length(out22sb[[2]]))
            rownames(Mat01sb)=names(out22sb[[2]])
            textplot(round(Mat01sb,4),cex.main=3,cex=2.5,family="sans")
            title("Up markers DRUGNEM scores")
            
            Mat01tb=matrix(out22tb[[2]],nrow=length(out22tb[[2]]))
            rownames(Mat01tb)=names(out22tb[[2]])
            textplot(round(Mat01tb,4),cex.main=3,cex=2.5,family="sans")
            title("Down markers DRUGNEM scores")
            
            Mat01ub=matrix(out22ub[[2]],nrow=length(out22ub[[2]]))
            rownames(Mat01ub)=names(out22ub[[2]])
            textplot(round(Mat01ub,4),cex.main=3,cex=2.5,family="sans")
            title("No effect markers DRUGNEM scores")
            
            Mat01v=matrix(out22v[[2]],nrow=length(out22v[[2]]))
            rownames(Mat01v)=names(out22v[[2]])
            textplot(round(Mat01v,4),cex.main=3,cex=2.5,family="sans")
            title("Best effect DRUGNEM scores")
            
            Mat1=matrix(out3[[2]],nrow=length(out3[[2]]))
            rownames(Mat1)=names(out3[[2]])
            textplot(round(Mat1,4),cex.main=3,cex=2.5,family="sans")
            title("Drug Independence effect Scores")
            
            Mat1s=matrix(out3s[[2]],nrow=length(out3s[[2]]))
            rownames(Mat1s)=names(out3s[[2]])
            textplot(round(Mat1s,4),cex.main=3,cex=2.5,family="sans")
            title("Up Drug Independence Scores")
            
            Mat1t=matrix(out3t[[2]],nrow=length(out3t[[2]]))
            rownames(Mat1t)=names(out3t[[2]])
            textplot(round(Mat1t,4),cex.main=3,cex=2.5,family="sans")
            title("Down Drug Independence Scores")
            
            Mat1u=matrix(out3u[[2]],nrow=length(out3u[[2]]))
            rownames(Mat1u)=names(out3u[[2]])
            textplot(round(Mat1u,4),cex.main=3,cex=2.5,family="sans")
            title("No effect Drug Independence Scores")
            
            Mat1sb=matrix(out3sb[[2]],nrow=length(out3sb[[2]]))
            rownames(Mat1sb)=names(out3sb[[2]])
            textplot(round(Mat1sb,4),cex.main=3,cex=2.5,family="sans")
            title("Up markers Drug Independence Scores")
            
            Mat1tb=matrix(out3tb[[2]],nrow=length(out3tb[[2]]))
            rownames(Mat1tb)=names(out3tb[[2]])
            textplot(round(Mat1tb,4),cex.main=3,cex=2.5,family="sans")
            title("Down markers Drug Independence Scores")
            
            Mat1ub=matrix(out3ub[[2]],nrow=length(out3ub[[2]]))
            rownames(Mat1ub)=names(out3ub[[2]])
            textplot(round(Mat1ub,4),cex.main=3,cex=2.5,family="sans")
            title("No effect markers Drug Independence Scores")
            
            Mat1v=matrix(out3v[[2]],nrow=length(out3v[[2]]))
            rownames(Mat1v)=names(out3v[[2]])
            textplot(round(Mat1v,4),cex.main=3,cex=2.5,family="sans")
            title("Best effect Drug Independence Scores")
            #Mat2=matrix(out4[[2]],nrow=length(out4[[2]]))
            #rownames(Mat2)=names(out4[[2]])
            #textplot(round(Mat2,4),cex.main=3,cex=2.5,family="sans")
            #title("Additivescore")
            
            #Mat2s=matrix(out4s[[2]],nrow=length(out4s[[2]]))
            #rownames(Mat2s)=names(out4s[[2]])
            #textplot(round(Mat2s,4),cex.main=3,cex=2.5,family="sans")
            #title("Weighted Additivescore")
            
            dev.off();

        }
        
        
    }
    return(list(DRUGNEMeffectmodel=rr$graph,Rankingeffectmodel=Mat01,DRUGNEMupmodel=rr1$graph,Rankingupmodel=Mat01s,DRUGNEMdownmodel=rr2$graph,Rankingdownmodel=Mat01t,DRUGNEMnoeffectmodel=rr66$graph,Rankingnoeffectmodel=Mat01u,DRUGNEMupmarkermodel=rr1b$graph,Rankingupmarkermodel=Mat01sb,DRUGNEMdownmarkermodel=rr2b$graph,Rankingdownmarkermodel=Mat01tb,DRUGNEMnoeffectmarkermodel=rr66b$graph,Rankingnoeffectmarkermodel=Mat01ub,DRUGNEMbesteffectmodel=rr77$graph,Rankingbesteffectmodel=Mat01v,Up_object=rr1,Down_object=rr2, Noeffect_object=rr66,Upmarker_object=rr1b,Downmarker_object=rr2b, Noeffectmarker_object=rr66b,Besteffect_object=rr77,Rankingindependencemodel=Mat1,Rankingupindependencemodel=Mat1s,Rankingdownindependencemodel=Mat1t,Rankingnoeffectindependencemodel=Mat1u,Rankingupmarkerindependencemodel=Mat1sb,Rankingdownmarkerindependencemodel=Mat1tb,Rankingnoeffectmarkerindependencemodel=Mat1ub,Rankingbestindependencemodel=Mat1v))
}
