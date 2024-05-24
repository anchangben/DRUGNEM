#' Determining the optimal prunning level L of the initial ccast tree
#'
#'This function takes a dataframe of biomarkers and cluster assignments and produces a ccast decision tree with L levels
#'
#' @param DD a dataframe of biomarker expression derived  from fcs files containing an extra column for cluster assignments labeled as "groups".
#' @param mm maximum number of predicted clusters
#' @return A pdf file labeled "initial_tree.pdf" corresponding to the initial ccast decision tree. In addition a list of 3 items:
#' \item{treenodes }{list of data of all leaf nodes of the ccast tree corresponding to homogeneous subgroups}
#' \item{level }{maximum level L for prunning decision tree}
#' \item{initialtree}{initial ccast tree output with less homogeneous bins in the leaf nodes} 
#' @seealso \code{\link{ccast_filter5}} which this function wraps
#' @export
#' @examples
#' ccastree<-ccast_L(dD1,mm)
#'
ccast_L <-
function(DD,mm){
a=1
for ( s in 1:10) {
     if(length(a) < mm) {
         tree1 <- ctree(formula=as.factor(groups) ~ ., data = DD,controls = ctree_control(minbucket = 50,maxdepth = s))
         t1=party::nodes(tree1, unique(where(tree1)))
         L=length(t1)
         a=unique(unlist(lapply(t1, function(x) which.max(x$prediction))))
         print(a)
     }
     else {
         print(a)
         break
     }
 }
    #pdf(file="initial_tree.pdf",width=15, height=10)
    #plot(tree1,xlab="Celltypes",main="Initial CCAST tree")
    #dev.off();
return(list(treenodes=t1,level=s-1,initialtree=tree1))
}

#################
ccast_L2 <-
function(DD,mm){
    a=1
    for ( s in 1:10) {
        if(length(a) < mm) {
            tree1 <- ctree(formula=as.factor(groups) ~ ., data = DD,controls = ctree_control(minbucket = 50,maxdepth = s))
            t1=party::nodes(tree1, unique(where(tree1)))
            L=length(t1)
            a=unique(unlist(lapply(t1, function(x) which.max(x$prediction))))
            print(a)
        }
        else {
            print(a)
            break
        }
    }
    return(list(treenodes=t1,level=s-1))
}