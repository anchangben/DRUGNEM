#' Scoring functions for all drug combinations using the nested effects models
#'
#' This function takes a nem graph object on drugs, a list of Drug tagets for each node and a data matrix with posterior weights for attaching a target protein to each node in the drug network. It produes a list of drug regimen(s) with the best score and a ranking of all drug combinations and their corresponding scores.
#'
#' @param G a graphnel object derived from nem.
#' @param mappos List of names for each drug node in the network corresponding to the proteins connected to each node also derived from a nem output.
#' @param data 2D matrix of posterior probabilities for each drug-protein pair derived from a nem output object.
#' @return List of the drug regimen(s) with the highest score and a ranking of all drug combinations and their corresponding scores:
#  \item{MAXcomb }{Drug regimen(s) with the highest score}
#' \item{Ordercomb}{Ranking of all drug combinations and their corresponding scores} 
#' @seealso \code{\link{drugmnemmain}} which this function wraps
#' @export
#' @examples
#' out22=bestdrugtarget(G=rr$graph,mappos=mappos,data=pos)
#'
################# nested combination scoring function drugmnmemscore 2 
bestdrugtarget<-function(G,mappos,data) {
    n=numNodes(G)
    COMB=list()
    for ( k in 1:(n-1)) {
        COMB[[k]]=combn(graph:::nodes(G),k)
        #names(COMB[[k]])<-apply(COMB[[k]], 2, function(x) paste(x,collapse=" "))
    }
    COMB[[n]]= matrix(graph:::nodes(G),ncol=1)
    ###names(COMB[[n]])<-apply(COMB[[n]], 2, function(x) paste(x,collapse=" "))
    res=list()
    for ( s in 1: (length(COMB)-1)) {
        ###for ( s in 1: length(COMB)) {
        res[[s]]<-apply(COMB[[s]],2,function(x) bestcombscore(G,mappos,data,x)[[1]])
        names(res[[s]])<-apply(COMB[[s]], 2, function(x) paste(x,collapse=" "))
    }
    scores=unlist(res)
    scores2=sort(scores,decreasing=TRUE)
    maxScore=scores[which.max(unlist(res))]
    
    return(list("MAXcomb"=maxScore,"Ordercomb"=scores2))
}

############
#### Calculating the expected additive effect under bliss independent ####
addbliss<-function(D) {
    if ((dim(D)[2])==1) {
        ee=(sum(D))
    } else {
        fmla<-as.formula(paste(" ~ 0+", paste("(",paste(colnames(D), collapse= "+"),")","^",dim(D)[2],sep="")))
        A= model.matrix(fmla,D)
        B=A[,1:dim(D)[2],drop=FALSE]
        C=A[,(dim(D)[2]+1):dim(A)[2],drop=FALSE]
        E1=apply(B,1,sum)
        E2=apply(C,1,function(x) sum(-x))
        ee=sum(E1+E2)
    }
    return(ee)
}

#### Using names of drugs for combination ###
scoredrugcombbliss<-function(D){
    COMB=list()
    for ( k in 1:(dim(D)[2]-1)) {
        COMB[[k]]=combn(colnames(D),k)
        
    }
    COMB[[dim(D)[2]]]= matrix(colnames(D),ncol=1)
    res=list()
    for ( s in 1: (length(COMB)-1)) {
        ###for ( s in 1: length(COMB)) {
        res[[s]]<-apply(COMB[[s]],2,function(x) addbliss(as.data.frame(D[,x])))
        names(res[[s]])<-apply(COMB[[s]], 2, function(x) paste(x,collapse=" "))
    }
    scores=unlist(res)
    scores2=sort(scores,decreasing=TRUE)
    maxScore=scores[which.max(unlist(res))]
    
    return(list("MAXcomb"=maxScore,"Ordercomb"=scores2))
}

##### Internal functions used to score and rank drug combinations #######
#### Minimum Independence score ########
Robf<- function(D) {
    ee1=apply(D,1,function(x) (1-(prod(1-x))))
    ee=sum(ee1)
    return(ee)
}
################

################
#### Minimum Idependence score: Using names of drugs for combination ###
scoredrugcomb2<-function(D){
    COMB=list()
    for ( k in 1:(dim(D)[2]-1)) {
        COMB[[k]]=combn(colnames(D),k)
        
    }
    COMB[[dim(D)[2]]]= matrix(colnames(D),ncol=1)
    res=list()
    for ( s in 1: (length(COMB)-1)) {
        ###for ( s in 1: length(COMB)) {
        res[[s]]<-apply(COMB[[s]],2,function(x) Robf(as.matrix(D[,x])))
        names(res[[s]])<-apply(COMB[[s]], 2, function(x) paste(x,collapse=" "))
    }
    scores=unlist(res)
    scores2=sort(scores,decreasing=TRUE)
    maxScore=scores[which.max(unlist(res))]
    
    return(list("MAXcomb"=maxScore,"Ordercomb"=scores2))
}
############################


################# Mutually exclusisve drugmnemscore 1####
drugtargetcount3<-function(G,mappos,pos) {
    n=numNodes(G)
    COMB=list()
    for ( k in 1:(n-1)) {
        COMB[[k]]=combn(graph:::nodes(G),k)
        #names(COMB[[k]])<-apply(COMB[[k]], 2, function(x) paste(x,collapse=" "))
    }
    COMB[[n]]= matrix(graph:::nodes(G),ncol=1)
    ###names(COMB[[n]])<-apply(COMB[[n]], 2, function(x) paste(x,collapse=" "))
    res=list()
    for ( s in 1: (length(COMB)-1)) {
        ###for ( s in 1: length(COMB)) {
        res[[s]]<-apply(COMB[[s]],2,function(x) combscore2(G,mappos,pos,x)[[1]])
        names(res[[s]])<-apply(COMB[[s]], 2, function(x) paste(x,collapse=" "))
    }
    scores=unlist(res)
    scores2=sort(scores,decreasing=TRUE)
    maxScore=scores[which.max(unlist(res))]
    
    return(list("MAXcomb"=maxScore,"Ordercomb"=scores2))
}
#################
############

########mutual exclusive sum posterior of total number of targets for drug combinations#####
combscore2<-function(G,mappos,pos,comb) {
  if (length(comb)==1) {
      tset= unique(scoreonedrug(comb,mappos,G))
      s.matrix2=tset
      
      scoretargetset=sum(apply(pos[tset,-dim(pos)[2],drop=FALSE], 1, max))
      } else {
    s.matrix2<-list()
    scoretargetset=NULL
for (ii in comb) {
      s.matrix2[[ii]]<-unique(scoreonedrug(ii,mappos,G))
      scoretargetset=sum(scoretargetset,sum(apply(pos[unlist(s.matrix2[[ii]]),-dim(pos)[2],drop=FALSE], 1, max)))
                    }
  tset=unique(unlist(s.matrix2))
     }
    
return(list(score=scoretargetset,combset=s.matrix2))
}

#### score using max of combining posterior effects under dependence and mutually exclusive assumption #####
bestcombscore<-function(G,mappos,data,comb) {
    data=data[,-dim(data)[2]]
    if (length(comb)==1) {
        tset= unique(scoreonedrug(comb,mappos,G))
        s.matrix2=tset
        
        scoretargetset=sum(apply(data[tset,,drop=FALSE], 1, max))
    } else {
        s.matrix2<-list()
        ##scoretargetset=NULL
        set1=NULL
        set2=NULL
        for (ii in comb) {
            
            s.matrix2[[ii]]<-unique(scoreonedrug(ii,mappos,G))
            
            ###if (length(intersect(s.matrix2[[ii]],set2))>0)  {
            ind1=s.matrix2[[ii]]
            ind2=unique(c(ind1,set2))
            score1=(sum(apply(data[ind1,,drop=FALSE],1,max)))
            score2=(sum(apply(data[ind2,,drop=FALSE],1,max)))
            set1=max(c(score1,score2))
            ##scoretargetset=sum(scoretargetset,length(unlist(s.matrix2[[ii]])))
            
            ###} else {
            ### set1=sum(c(set1,sum(apply(data[unlist(s.matrix2[[ii]]),,drop=FALSE],1,max))))
            ###set1=sum(c(set1,sum(data[unlist(s.matrix2[[ii]]),ii,drop=FALSE])))
            
            set2=c(set2,s.matrix2[[ii]])
        }
        scoretargetset=set1
        tset=unique(unlist(s.matrix2))
    }

    
    return(list(score=scoretargetset,combset=s.matrix2))
}
########
########

##### Calculate the set of targets for each node #####
scoreonedrug<-function(id,mappos,G) {
    m<-graph2adj(G)
    n=numNodes(G)
    m1 <- diag(n)
    m2<-matrix(0,nrow=n,ncol=2*n)
    rownames(m1)<-c(letters[1:n])
    rownames(m2)<-c(letters[1:n])
    colnames(m1)<-c(letters[1:n])
    m3<-cbind(m,m1)
    m4<-rbind(m3,m2)
            s.matrix<-list()
            s.matrix2<-list()
            i= match(id,graph:::nodes(G))
            s.matrix[[i]]<-list()
            s.matrix2[[i]]<-list()
            targetset=NULL
            for (j in 1:n) {
                #print("j=",j)
                x1=paths(m4,i,n+j)
                x2=which(x1==i)
                if((length(x2)==1)|(length(x2)==0)) {
                    
                    if (length(x2)==1)  {
                        R2<-numtarget2(x1[-length(x1)],mappos)
                        s.matrix[[i]][[j]]<-x1
                        s.matrix2[[i]][[j]]<-R2
                        DrugTarget<-R2
                        targetset=c(targetset,DrugTarget)
                        #print(targetset)
                    }
                    if (length(x2)==0)  {
                        R2<-NULL
                        s.matrix[[i]][[j]]<-list()
                        s.matrix2[[i]][[j]]<-NULL
                        DrugTarget<-R2
                        targetset=c(targetset,DrugTarget)
                        # print(targetset)
                    }
                } else {
                    
                    R=list()
                    R2=NULL
                    for ( k in 1: (length(x2))){
                        if(k<length(x2))  {
                            R[[k]]<-x1[x2[k]:(x2[k+1]-1)]
                            R2<-c(R2,numtarget2(R[[k]][-length(R[[k]])],mappos))
                        } else {
                            R[[k]]<-x1[x2[k]:length(x1)]
                            R2<-c(R2,numtarget2(R[[k]][-length(R[[k]])],mappos))
                        }
                    }
                    s.matrix[[i]][[j]]<-R
                    s.matrix2[[i]][[j]]<-R2
                    DrugTarget<-R2
                    targetset=c(targetset,DrugTarget)
                    #print(targetset)
                }
            }
            
            return(targetset)
        }
############
######## Find all paths between 2 nodes using igraph ########
paths <- function (amat, st, en, path = c()){
    indices = 1:nrow(amat)
    if(st == en)     return(c(path, st)) ### st is 'node' in recursive calls
    if(sum(amat[st,]) == 0 )
    return(NULL)
    paths = c()
    ne = indices[amat[st,]==1] #### Boundary of x. Assumes that amat is symmetric
    
    for(node in ne){
        
        if(!is.element(node, c(path, st))){
            newpaths = paths(amat, node, en, c(path, st))
            for(newpath in newpaths){
                paths = c(paths,newpath)
            }
        }
        
        
    }
    return(paths)
}


numtarget2<-function(vec,mappos) {
            rr=NULL
            for (ii in 1:length(vec)) {
                rr=unique(c(rr,mappos[[vec[ii]]]))
            }
            return(rr)
        }

#### A function to turn an adjacency matrix into a graph ######
adj2graph <- function(adj.matrix) {
    V   <- rownames(adj.matrix)
    edL <- vector("list", length=nrow(adj.matrix))
    names(edL) <- V
    for (i in 1:nrow(adj.matrix)) {
        edL[[i]] <- list(edges=which(!adj.matrix[i,]==0),
        weights=adj.matrix[i,!adj.matrix[i,]==0])
    }
    gR <- new("graphNEL",nodes=V,edgeL=edL,edgemode="directed")
    return(gR)
}

##
#### A function to turn a graph into an adjacency matrix.
#### Written by Benedict Anchang.
##

graph2adj <- function(gR) {
	require(graph)
    adj.matrix <- matrix(0,
    length(graph:::nodes(gR)),
    length(graph:::nodes(gR))
    )
    rownames(adj.matrix) <- graph:::nodes(gR)
    colnames(adj.matrix) <- graph:::nodes(gR)
    for (i in 1:length(graph:::nodes(gR))) {
        adj.matrix[graph:::nodes(gR)[i],adj(gR,graph:::nodes(gR)[i])[[1]]] <- 1
    }
    
    return(adj.matrix)
}

#### LIMMA ######
make.R.matrix1<-function (dat, wt, pi1 = 0.01)
{
    ###require(limma)
    cn <- colnames(dat)
    cn <- setdiff(cn, wt)
    colnames(dat) <- make.names(colnames(dat))
    exps <- factor(colnames(dat))
    design <- model.matrix(~0 + factor(colnames(dat)))
    fit1 <- lmFit(dat, design)
    contrast.matrix <- makeContrasts(contrasts = paste(setdiff(exps,
    wt), "-", wt, sep = ""), levels = exps)
    fit2 <- contrasts.fit(fit1, contrast.matrix)
    fit3 <- eBayes(fit2,proportion=pi1)
    lods <- fit3$lods
    tscore<-fit3$t
    pscore<-fit3$p.value
    pdiff <- exp(lods)/(1+exp(lods))
    out=topTable(fit3,number=100)
    colnames(lods) <- cn
    res<-list(logodds=lods,tstat=tscore,pvalue=pscore,top=out,prob=pdiff)
    return(res)
}

##### LIMMA one-sided #######
######### One sided limma ######
limmaonesidedtest<-function (fit, lower = TRUE)
{
    se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled
    df.total <- fit$df.prior + fit$df.residual
    pt(fit$t, df = df.total, lower.tail = lower)
}


make.R.matrix2<-function (dat, wt, pi1 = 0.01)
{
    ###require(limma)
    cn <- colnames(dat)
    cn <- setdiff(cn, wt)
    colnames(dat) <- make.names(colnames(dat))
    exps <- factor(colnames(dat))
    design <- model.matrix(~0 + factor(colnames(dat)))
    fit1 <- lmFit(dat, design)
    contrast.matrix <- makeContrasts(contrasts = paste(setdiff(exps,
    wt), "-", wt, sep = ""), levels = exps)
    fit2 <- contrasts.fit(fit1, contrast.matrix)
    fit3 <- eBayes(fit2,proportion=pi1)
    #### test upregulation for each treatment #####
    ft<-limmaonesidedtest(fit3, lower = FALSE)
    lods <- fit3$lods
    tscore<-fit3$t
    pscore<-fit3$p.value
    pdiff <- exp(lods)/(1+exp(lods))
    out=topTable(fit3,number=100)
    colnames(lods) <- cn
    res<-list(logodds=lods,tstat=tscore,pvalue=pscore,top=out,prob=pdiff,downregpvalue=ft)
    return(res)
}

make.R.matrix3<-function (dat, wt, pi1 = 0.01)
{
    ###require(limma)
    cn <- colnames(dat)
    cn <- setdiff(cn, wt)
    colnames(dat) <- make.names(colnames(dat))
    exps <- factor(colnames(dat))
    design <- model.matrix(~0 + factor(colnames(dat)))
    fit1 <- lmFit(dat, design)
    contrast.matrix <- makeContrasts(contrasts = paste(setdiff(exps,
    wt), "-", wt, sep = ""), levels = exps)
    fit2 <- contrasts.fit(fit1, contrast.matrix)
    fit3 <- eBayes(fit2,proportion=pi1)
    #### test upregulation for each treatment #####
    ft<-limmaonesidedtest(fit3, lower = TRUE)
    lods <- fit3$lods
    tscore<-fit3$t
    pscore<-fit3$p.value
    pdiff <- exp(lods)/(1+exp(lods))
    out=topTable(fit3,number=100)
    colnames(lods) <- cn
    res<-list(logodds=lods,tstat=tscore,pvalue=pscore,top=out,prob=pdiff,downregpvalue=ft)
    return(res)
}
###### end one-sided #####

#####
### distribution of bivariate logistic distribution with logistic marginals under independence ###
logoddsf<-function(x,y) {
    a=exp(x+y)
    b=1+exp(x)+exp(y)+a
    z=a/b
    z1=log(z/(1-z))
    return(z1)
}

logoddsf2<-function(x,y) {
    a=exp(x+y)
    b=1+exp(x)+exp(y)+a
    c=exp(y)/(1+exp(y))
    z=c-(a/b)
    z1=log(z/(1-z))
    return(z1)
}

### NEM CI with no prior on NULL S-gene ####
nemcombindex<-function(nemr) {
    pos=nemr$pos[[which.max(nemr$mLL)]]
    mappos=nemr$mappos[[which.max(nemr$mLL)]]
    A=graph2adj(nemr$graph)
    diag(A)=1
    B=t(A)
    pos1=pos[,-dim(pos)[2]]
    nn=combn(colnames(A),2)
    C=matrix(0, nrow=dim(A)[1],ncol=dim(A)[2])
    rownames(C)=c(paste("E", colnames(A)))
    colnames(C)<-apply(nn, 2, function(x) paste(x,collapse="+"))
    for (j in 1:dim(nn)[2]) {
        C[,j]=ifelse(((B[,nn[1,j],drop=FALSE]+B[,nn[2,j],drop=FALSE])>=1),1,0)
    }
    print(B)
    print(C)
    TM0=NULL
    for ( j in 1:(length(mappos)-1)) {
        if (length(mappos[[j]])>0) {
            for ( r in 1:length(mappos[[j]])){
                TM0=rbind(TM0,C[j,,drop=FALSE])
            }
            
        } else{
            TM1=NULL
            TM0=rbind(TM0,TM1)
        }
    }
    for ( j in length(mappos)) {
        if (length(mappos[[j]])>0) {
            for ( r in 1:length(mappos[[j]])){
                TM0=rbind(TM0,pos1[mappos[[j]][r],,drop=FALSE])
                ##TM0=rbind(TM0,pos1[mappos[[j]][r],,drop=FALSE]/sum(pos1[mappos[[j]][r],,drop=FALSE]))
            }
            ###rownames(TM0)=mappos[[j]]
        } else{
            TM1=NULL
            TM0=rbind(TM0,TM1)
        }
    }
    rownames(TM0)=unlist(mappos)
    numdata=length(TM0)
    ss=abs(rnorm(numdata,0,0.05))
    
    ddat32=ifelse(TM0==0,TM0+ss,TM0)
    ddat33=ifelse(ddat32==1,ddat32-ss,ddat32)
    combindexdata=ddat33
    ###combindexdata=TM0*dat[rownames(TM0),]
    return (combindexdata)
}

##### NEM CI with prior on NULL S-gene####
nemcombindex2<-function(nemr) {
    pos=nemr$pos[[which.max(nemr$mLL)]]
    mappos=nemr$mappos[[which.max(nemr$mLL)]]
    A=graph2adj(nemr$graph)
    diag(A)=1
    B=t(A)
    pos1=pos[,-dim(pos)[2]]
    nn=combn(colnames(A),2)
    C=matrix(0, nrow=dim(A)[1],ncol=dim(A)[2])
    rownames(C)=c(paste("E", colnames(A)))
    colnames(C)<-apply(nn, 2, function(x) paste(x,collapse="+"))
    for (j in 1:dim(nn)[2]) {
        C[,j]=ifelse(((B[,nn[1,j],drop=FALSE]+B[,nn[2,j],drop=FALSE])>=1),1,0)
    }
    print(B)
    print(C)
    TM0=NULL
    for ( j in 1:(length(mappos)-1)) {
        if (length(mappos[[j]])>0) {
            for ( r in 1:length(mappos[[j]])){
                TM0=rbind(TM0,C[j,,drop=FALSE])
            }
            
        } else{
            TM1=NULL
            TM0=rbind(TM0,TM1)
        }
    }
    for ( j in length(mappos)) {
        if (length(mappos[[j]])>0) {
            ss=abs(rnorm((length(mappos[[j]])*dim(C)[2]),0,0.05))
            CC=matrix(ss,nrow=length(mappos[[j]]),ncol=dim(C)[2])
            TM0=rbind(TM0,CC)
            ###TM0=rbind(TM0,pos1[mappos[[j]][r],,drop=FALSE]/sum(pos1[mappos[[j]][r],,drop=FALSE]))
            ###rownames(TM0)=mappos[[j]]
        } else{
            TM1=NULL
            TM0=rbind(TM0,TM1)
        }
    }
    rownames(TM0)=unlist(mappos)
    combindexdata=TM0
    return (combindexdata)
}

#### Simes test function
simes<-function (x, returnstat = FALSE,alpha=0.05)
{
    r = rank(x)
    T = min(length(x) * x/r)
    ###cutoff=which(c(length(x) * x/r)==T)
    if (returnstat)
    c(T, alpha/length(x))
    else T
}

###Bonferonni test function
bonfer<-function(x,returnstat = FALSE,alpha=0.05){
    require(graphics)
    p1=p.adjust(x, method = "bonferroni")
    T = length(x)* min(p1)
    if (returnstat)
    c(T, alpha/length(x))
    else T
    
}


