#' Plotting Biaxial plots for all inner nodes of ccast decision tree output
#'
#' This function takes as input an optimal ccast decision tree output from the ccast_optimize function and any biomarker label not used in building the tree provided by the user and produces a pdf file with all pairwise scatter plots of all the inner nodes in the decision tree with density contours and split point cut off estimates used in building the tree.
#'
#' @param optccastreeoutput An optimal ccast tree object produced by the function ccast_optimize
#' @param ylabel Any biomarker label provided by the user preferably one that is not used to build the tree
#' @return A pdf file labeled "Biaxial tree node plots.pdf" with biaxial plots for all inner nodes from the optimal decision tree.
#' @seealso \code{\link{ccast_filter5}} which this function wraps
#' @export
#' @examples
#' ylabel<-"M1"
#' data(optccastree)
#' ccast_biaxialplot(optccastreeoutput=optccastree,ylabel=ylabel)
#'
ccast_biaxialplot <-
function(optccastreeoutput,ylabel) {
    ##tiff(file="CCAST Biaxial plots%d.tiff",width=360, height=360)
    ##require(MASS)
    pdf(file="Biaxial tree node plots.pdf",width=10, height=10)
    cex <- 1.6
    tree33<-optccastreeoutput$tree
    dD33=optccastreeoutput$puredata[[length(optccastreeoutput$puredata)]]
    inner <- party::nodes(tree33, c(1:max(unique(where(tree33))))[-unique(where(tree33))])
    ##layout(matrix(1:length(inner), ncol = round(sqrt(length(inner)),0)))
    ###par(mfrow=c(round(s/2)+1,2))
    
    out <- sapply(inner, function(i) {
       splitstat <- i$psplit$splitstatistic
        w1=which(i$weights==1)
        x <- dD33[[i$psplit$variableName]][w1]
        y<-dD33[[ylabel]][w1]
        tryCatch({
            if (bandwidth.nrd(x)>0&bandwidth.nrd(y)>0) {
                xydens <- kde2d(x,y,n=25,lims = c(range(x), range(y)))
            }else {
                print("Bandwidths=0")
                print(c(bandwidth.nrd(x), bandwidth.nrd(y)))
                xydens <- kde2d(x,y,h=0.5,n=25,lims = c(range(x), range(y)))
            }

        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        plot(x,y, main = paste("Node", i$nodeID),
        xlab = i$psplit$variableName, ylab = ylabel, ylim = c(0, max(y)+5),
        cex.axis = cex, cex.lab = cex, cex.main = cex,cex=0.2, pch=20,col="dark blue")
        abline(v = i$psplit$splitpoint, lwd=4,lty = 3,col="red")
        #####filled.contour(xydens, color = terrain.colors, asp = 1) # simple
        contour(xydens,nlevels=10,add=T,col="dark orange")
     })
    
    dev.off();

}
