#' Optimize the initial decision tree with height L for homogenous clusters.
#'
#' @section description: This function optimizes the decision tree by maximizing the size of the homogenous clusters.
#'
#' @param t1 A tree object from the function ccast_L corresponding to homogeneous subgroups in the leaf nodes of the initial decision tree.
#' @param dD1 a dataframe of biomarker expression derived  from fcs files containing an extra column for cluster assignments labeled as "groups".
#' @param s maximum level L for prunning decision tree
#' @return 	A pdf file labeled "finaltree.pdf" corresponding to the final ccast decision tree. In addition a list of 3 items:
#' \item{tree }{final ccast tree output with distinct homogeneous bins in the leaf nodes.}
#' \item{fprediction }{list representing the distributions of events in each leaf node in the optimized ccast decision tree.}
#' \item{puredata}{list of all updated filtered data corresponding to each iteration in the optimization step. } 
#' @seealso \code{\link{ccast_filter5}} which this function wraps
#' @export
#' @examples
#' optccastree<-ccast_optimize(t1,dD1,s)
#'
ccast_optimize <-
function(t1,dD1,s){

last=dim(dD1)[2]
if(colnames(dD1)[last]!="groups") colnames(dD1)[last]<-"groups"
b=unique(unlist(lapply(t1, function(x) max(x$prediction))))
iter=1
DD=list()
DD[[iter]]=dD1
minleaf=1:max(dD1[,last])
    if (length(b)==1) {
         tree2 <- ctree(formula=as.factor(groups) ~ ., data = dD1,controls = ctree_control(maxdepth = s))
        save(tree2,file="finaltree.rdata")
        pdf(file="finaltree.pdf",width=15, height=10)
        plot(tree2,xlab="Celltypes",main="Final CCAST tree")
        dev.off();
    } else {

while(((sum(b==1)!=length(b))&(iter<=25))) {
print(paste("iteration: ",iter))
print(paste("Distribution of cells at leaf nodes"))
w3=NULL
minl=NULL
  for( j in 1:length(t1)) {
  mx=which.max(t1[[j]]$prediction)
  w1=which(t1[[j]]$weights==1)
  w2=w1[which(DD[[iter]][[last]][w1]==mx)]
  minl2=names(table(DD[[iter]][[last]][w2]))
  minl=c(minl,minl2)
  print(table(DD[[iter]][[last]][w2]))
 w3=c(w3,w2)
 }
 if(length(unique(minl))==1) {
 	tree2 <- ctree(formula=as.factor(groups) ~ ., data = dD1,controls = ctree_control(maxdepth = s))
 	
 break
 } else { 
if(all(minleaf %in% minl)) { 
 iter=iter+1
DD[[iter]]=as.data.frame(as.matrix(DD[[iter-1]])[w3,,drop=FALSE])
tree2 <- party::ctree(as.factor(groups) ~ ., data = DD[[iter]],controls = ctree_control(maxdepth = s)) 
t1=party::nodes(tree2, unique(where(tree2)))
b=unique(unlist(lapply(t1, function(x) max(x$prediction))))
} else {
DD[[iter]]=DD[[iter-1]]
tree2 <- party::ctree(as.factor(groups) ~ ., data = DD[[iter]],controls = ctree_control(maxdepth = s)) 
    t1=party::nodes(tree2, unique(where(tree2)))
    b=unique(unlist(lapply(t1, function(x) max(x$prediction))))

break
}
}
}
save(tree2,file="finaltree.rdata")
pdf(file="finaltree.pdf",width=15, height=10)
plot(tree2,xlab="Celltypes",main="Final CCAST tree")
dev.off();
        
 }
return(list(tree=tree2,fprediction=b,puredata=DD))
}
