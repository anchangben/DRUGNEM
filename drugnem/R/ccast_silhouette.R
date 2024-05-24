#' Optimizing the density of cell subpopulations.
#'
#' This function uses Silhouette coefficients to optimize both the number of cell clusters as well as cell cluster distributions.
#'
#' @param d distant measurement. Object  dist()
#' @param fileout pdf output file. Default labeled "silhouette.pdf"
#' @param nCluster.init number of initial clusters
#' @param assign cluster assignments from any clustering algorithm
#' @param iter.max max number of iteration
#' @param diagnosis logical for ploting silouette diagnostic plots. Default is TRUE
#' @return 	A pdf file labeled "silhouette.pdf" corresponding to the silhouette diagnostic plots. In addition:
#' \item{orgFit }{Output from initial silhouette run}
#' \item{optFit }{Output from final silhouette run}
#' \item{optFitCluster }{Dataframe of final cluster assignments, final optimal number of clusters, optimal number of iteration and optimal cost function statistics} 
#' @seealso \code{\link{ccast_filter5}} which this function wraps
#' @export
#' @examples
#' assign <- cutree(hc, k=4)
#' out1=ccast_silhouette(d,fileout="silhouette.pdf", nCluster.init = 4, assign, iter.max = 10, diagnosis = TRUE)

ccast_silhouette <-
function(d,fileout, nCluster.init = k,assign, iter.max = 10, diagnosis = TRUE){
        suppressMessages(require(ggplot2))
        
        groups = assign #plot(h); rect.hclust(h, k=10, border="red")
        sk0 = silhouette(groups, d)
        #rownames(sk0) = h$labels # h$labels == names(groups) == rownames(d)
        ## profiling: choose opt value
        sk=sk0
        i = 1
        nc = vector()
        sw = vector()
        ssw = vector()
        while(any(sk[,3] < 0) && i <= iter.max){
            nc[i] = length(unique(sk[, "cluster"]))
            ssw[i] = sum(sk[,"sil_width"][sk[,"sil_width"] < 0])
            assign2 = ifelse(sk[,"sil_width"]<0, sk[,"neighbor"], sk[,"cluster"])
            sk <- silhouette(assign2, d)
            i = i + 1
        }
        print(ssw)
        ssw.opt = max(ssw)
        p1 = ggplot(data.frame("numberOfClusters" = nc, "Iteration" =  seq_along(nc)), aes(x=Iteration, y = numberOfClusters)) +
        geom_point() +
        labs(x = "Iteration", y = "Number of clusters", title="Profiling")
        dmat = data.frame("ssw" = ssw, "Iteration" =  seq_along(ssw))
        p2 = ggplot(dmat, aes(x=Iteration, y = ssw)) +
        geom_point() +
        geom_point(data = subset(dmat, ssw == ssw.opt), size = 5, colour = alpha("red", 1/2)) +
        labs(x = "Iteration", y = "Sum of negative Silhouette width", title="Profiling")
        
        ## opt fit
        sk=sk0
        i = 1
        sw = -1
        ssw = -1000000
        nc = vector()
        ssw.vec = vector()
        while(ssw < ssw.opt && i <= iter.max){
            nc.opt = length(unique(sk[, "cluster"]))
            nc[i] = nc.opt
            ssw = sum(sk[,"sil_width"][sk[,"sil_width"] < 0])
            ssw.vec[i] = ssw
            assign2 = ifelse(sk[,"sil_width"]<0, sk[,"neighbor"], sk[,"cluster"])
            sk <- silhouette(assign2, d)
            i = i + 1
        }
        ##rownames(sk) = rownames(sk0)
        
        p3 = ggplot(data.frame("numberOfClusters" = nc, "Iteration" =  seq_along(nc)), aes(x=Iteration, y = numberOfClusters)) +
        geom_point() +
        labs(x = "Iteration", y = "Number of clusters", title="Iteration to optimal fit")
        dmat = data.frame("ssw" = ssw.vec, "Iteration" =  seq_along(ssw.vec))
        p4 = ggplot(dmat, aes(x=Iteration, y = ssw.opt)) +
        geom_point() +
        geom_point(data = subset(dmat, ssw == ssw.opt), size = 5, colour = alpha("red", 1/2)) +
        geom_text(data = subset(dmat, ssw == ssw.opt), aes(x = Iteration, y = ssw, label = paste("(", Iteration,",", round(ssw.opt,2), ")", sep = "")), hjust=1, vjust = -1) +
        #	geom_text(data = subset(dmat, ssw == ssw.opt), x = 0, y = ssw.opt, label = round(ssw.opt,2)) +
        labs(x = "Iteration", y = "Sum of negative silhouette width", title="Iterate to optimal fit") +
        ylim(c(min(ssw.vec), ifelse(max(ssw.vec) > 0, ssw.vec, 0)))# + xlim(c(0, length(nc)+5))
        if(diagnosis){
            pdf(fileout)
            plot(sk0, main = "Silhouette plot: Original clusters")
            print(p1)
            print(p2)
            plot(sk,  main = "Silhouette plot: Final clusters")
            print(p3)
            #print(p4)
            dev.off()
        }
        
        gc()
        
        return(list("orgFit" = sk0, "optFit" = sk,"optFitCluster" = data.frame("cluster" = sk[,"cluster"], stringsAsFactors = F),"optNumCluster" = nc.opt, "stopAtIteration" = i, "optSumOfNegSW" = ssw.opt))
    }
