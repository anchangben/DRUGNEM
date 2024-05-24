#' plot graph of nested effects model, the marginal likelihood distribution or the posterior position of the effected proteins
#'
#' @param  x nem object to plot
#' @param  what (i), "graph", (ii) "mLL" = likelihood distribution, (iii) "pos" = posterior position of effected proteins
#' @param remove.singletons remove unconnected nodes from the graph plot
#' @param PDF output as PDF-file
#' @param filename filename of PDF-file
#' @param thresh if x has a real valued adjacency matrix (weight matrix), don't plot edges with |weight| <= thresh
#' @param transitiveReduction plot a transitively reduced graph
#' @param plot.probs plot edge weights/probabilities.
#' @param SCC plot the strongly connected components graph
#' @param D Visualize the nested subset structure of the dataset via plotEffects along with the graph and show the linking of E-proteins to S-proteins in the dataset. Should only be used for small networks. Default: Just plot the graph
#' @param draw.lines If the nested subset structure is shown, should additionally lines connecting S-genes and their associated E-genes be drawn? WARNING: For larger datasets than e.g. 5 S-genes this most probably does not work, because the nested subset structure picture then partially overlaps with the graph picture. Default: Do not draw these lines
#' @param palette olor palette to use: either 'BlueRed' (default) or 'Grey'
#' @param ... other arguments to be passed to the Rgraphviz plot function or to the graphics 'image' function.
#' @return 	None
#' @seealso \code{\link{fitmnemdown}} which this function wraps
#' @export
#' @examples
#' plot.nem2(rr,filename=paste(patient,"DRUGMNEM network and Heatmap for all states",infer,type,".pdf",sep=""),main="DRUGMNEM Network",PDF=TRUE,what = "graph",D=cdata1,draw.lines = FALSE)
#########
plot.nem2<-
function (x, what = "graph", remove.singletons = FALSE, PDF = FALSE, 
    filename = "nemplot.pdf", thresh = 0, transitiveReduction = FALSE, 
    plot.probs = FALSE, SCC = TRUE, D = NULL, draw.lines = FALSE, 
    palette = "BlueRed", ...) 
{
    if (!(what %in% c("graph", "mLL", "pos"))) 
        stop("\nnem> invalid plotting type: plot either 'graph', 'mLL', or 'pos'")
    if (what == "graph") {
        gR = x$graph
        if (numEdges(gR) == 0) 
            stop("Graph contains no edges - nothing to draw!")
        M = as(gR, "matrix")
        toremove = which((abs(M) <= thresh) & (abs(M) > 0), arr.ind = TRUE)
        if (nrow(toremove) > 0) {
            for (i in 1:nrow(toremove)) gR = removeEdge(from = graph:::nodes(gR)[toremove[i,
                1]], to = graph:::nodes(gR)[toremove[i, 2]], gR)
        }
        if (SCC) {
            gR = SCCgraph(gR)$graph
            M = as(gR, "matrix")
        }
        if (numEdges(gR) == 0) 
            edgeattr = list()
        else {
            if (transitiveReduction) 
                M = transitive.reduction(M)
            eDDn <- names(edgeDataDefaults(gR))
            if (!"weight" %in% eDDn) 
                edgeDataDefaults(gR, "weight") <- 1
            if (!"label" %in% eDDn) 
                edgeDataDefaults(gR, "label") <- 1
            if (!"arrowhead" %in% eDDn) 
                edgeDataDefaults(gR, "arrowhead") = "normal"
            if (!"style" %in% eDDn) 
                edgeDataDefaults(gR, "style") = "bold"
            nodes <- colnames(M)
            nodenames = vector("character", length(M[abs(M) > 
                0]))
            probs = double(length(nodenames))
            arr = character(length(probs))
            penwidth = rep("bold", length(probs))
            k = 1
            for (i in 1:ncol(M)) {
                for (j in 1:nrow(M)) {
                  if (M[i, j] != 0) {
                    if (class(x) != "dynoNEM") 
                      probs[k] = signif(ifelse(abs(M[i, j]) > 
                        1, abs(M[i, j]) - 1, abs(M[i, j])), 2)
                    else probs[k] = signif(M[i, j], 2)
                    edgeData(gR, from = nodes[i], to = nodes[j], 
                      attr = "style") = "bold"
                    edgeData(gR, from = nodes[i], to = nodes[j], 
                      attr = "label") = probs[k]
                    edgeData(gR, from = nodes[i], to = nodes[j], 
                      attr = "weight") = M[i, j]
                    if ((M[i, j] > 0) & (M[i, j] <= 1) || class(x) == 
                      "dynoNEM") {
                      edgeData(gR, from = nodes[i], to = nodes[j], 
                        attr = "arrowhead") = "normal"
                      arr[k] = "normal"
                    }
                    else if (M[i, j] > 1 & class(x) != "dynoNEM") {
                      edgeData(gR, from = nodes[i], to = nodes[j], 
                        attr = "arrowhead") = "vee"
                      arr[k] = "vee"
                    }
                    else {
                      edgeData(gR, from = nodes[i], to = nodes[j], 
                        attr = "arrowhead") = "tee"
                      arr[k] = "tee"
                    }
                    nodenames[k] <- paste(nodes[i], "~", nodes[j], 
                      sep = "")
                    k = k + 1
                  }
                  else {
                    if (nodes[i] %in% unlist(inEdges(nodes[j], 
                      gR))) 
                      gR = removeEdge(from = nodes[i], to = nodes[j], 
                        gR)
                  }
                }
            }
            names(arr) = nodenames
            names(probs) = nodenames
            names(penwidth) = nodenames
            fontcol = arr
            fontcol[arr == "tee"] = "blue"
            fontcol[arr == "normal"] = "black"
            fontcol[arr == "vee"] = "red"
            if (plot.probs) 
                edgeattr = list(label = probs, arrowhead = arr, 
                  fontcolor = fontcol, color = fontcol, style = penwidth)
            else edgeattr = list(arrowhead = arr, fontcolor = fontcol, 
                color = fontcol, style = penwidth)
        }
        edge.type = ifelse(plot.probs | class(x) == "dynoNEM", 
            "distinct", "combined")
        el = buildEdgeList(gR, recipEdges = edge.type, edgeAttrs = edgeattr)
        nodeattr = list(color = rep("white", length(graph:::nodes(gR))),
            penwidth = rep(0, length(graph:::nodes(gR))), fontsize = rep(14,
                length(graph:::nodes(gR))))
        names(nodeattr$color) = graph:::nodes(gR)
        names(nodeattr$penwidth) = graph:::nodes(gR)
        names(nodeattr$fontsize) = graph:::nodes(gR)
        args = list(...)
        if ("nodeAttrs" %in% names(args)) 
            nodeattr = c(nodeattr, args[[match("nodeAttrs", names(args))]])
        if ("edgeAttrs" %in% names(args)) 
            edgeattr = c(edgeattr, args[[match("edgeAttrs", names(args))]])
        main = NULL
        if ("main" %in% names(args)) 
            main = args[["main"]]
        G = agopen(gR, name = "test", edges = el, edgeAttrs = edgeattr, 
            nodeAttrs = nodeattr, recipEdges = edge.type)
        if (PDF) 
            pdf(file = filename)
        par(cex.main = 2, cex = 1)
        if (is.null(D)) 
            plot(G, main = main)
        else {
            zlim = NULL
            if ("zlim" %in% names(args)) 
                zlim = args[["zlim"]]
            plotnem2(D, G, x, SCC = SCC, main = main, zlim = zlim,
                draw.lines = draw.lines, palette = palette)
        }
        if (PDF) 
            dev.off()
        save(gR, file = paste(unlist(strsplit(filename, ".pdf")), 
            ".rda", sep = ""))
        toDotR(gR, paste(unlist(strsplit(filename, ".pdf")), 
            ".dot", sep = ""))
    }
    if (what == "mLL") {
        if (PDF) 
            pdf(file = filename)
        par(cex = 1.3)
        ss <- sort(unique(x$mLL), decreasing = TRUE)[1:min(30, 
            length(x$mLL))]
        plot(x = 1:length(ss), y = ss, pch = 19, main = "Score distribution", 
            xlab = paste(length(ss), "top ranked models"), ylab = "Marginal log-likelihood", 
            ylim = c(ss[length(ss)] - 10, ss[1] + 10))
        points(1, max(unique(x$mLL)), pch = 21, cex = 1.7, lwd = 2)
        if (PDF) 
            dev.off()
    }
    if (what == "pos") {
        if (length(x$mLL) > 1 & class(x) != "dynoNEM") {
            winner <- which.max(x$mLL)
            pos <- x$pos[[winner]]
            effects <- rownames(x$pos[[winner]])
        }
        else {
            pos <- x$pos
            effects <- rownames(x$pos)
        }
        pos[is.na(pos)] = 0
        if (PDF) 
            pdf(file = filename)
        par(las = 2, mgp = c(5.5, 1, 0), mar = c(6.7, 7, 4, 1), 
            cex.lab = 1.3, cex.main = 1.7)
        image(x = 1:ncol(pos), y = 1:nrow(pos), z = t(pos), main = "Posterior effect positions", 
            xlab = "Perturbations", xaxt = "n", ylab = "Effect reporters", 
            yaxt = "n", col = gray(seq(0.95, 0, length = 10)))
        abline(v = (1:(ncol(pos) - 1)) + 0.5)
        axis(1, 1:ncol(pos), colnames(pos))
        axis(2, 1:length(effects), effects)
        if (PDF) 
            dev.off()
    }
}


plotEffects2<-
function (D, nem, border = TRUE, legend = TRUE, order = NULL, 
    orderSCC = TRUE, palette = "BlueRed", ...) 
{
    if (!(class(D) %in% c("matrix", "data.frame"))) 
        stop("First argument has to be the data matrix and second the nem object!")
    sccg <- SCCgraph(nem$graph, name = TRUE)
    if (numEdges(sccg$graph)==0)  {
         topo.order <- graph:::nodes(sccg$graph)
   } else {
   topo.order <- tsort(sccg$graph)
    }
    if (is.null(order)) 
        myorder = topo.order
    else {
        if (!orderSCC) 
            myorder = rev(unique(sccg$which.scc[match(order, 
                names(sccg$which.scc))]))
        else myorder = rev(names(sccg$scc[order]))
    }
    if (class(nem) == "score") {
        mappos = nem$mappos[[which.max(nem$mLL)]]
    }
    else {
        mappos = nem$mappos
        if (length(mappos) == 1) 
            mappos = mappos[[1]]
    }
    selected = nem$selected
    null.genes = unique(unlist(mappos["null"], use.names = FALSE))
    if (length(null.genes) > 0) {
        sccg$scc[["null"]] = "null"
        myorder = c("null", myorder)
    }
    v <- list()
    nr <- list()
    for (i in myorder) {
        w = unique(unlist(mappos[sccg$scc[[i]]], use.names = FALSE))
        if (length(w) == 0) {
            v[[i]] <- NA
            nr[[i]] <- 1
        }
        if (length(w) == 1) {
            v[[i]] <- w
            nr[[i]] <- 1
        }
        if (length(w) > 1) {
            d <- dist(D[w, , drop = FALSE], method = "manhattan")
            v[[i]] <- w[hclust(d)$order]
            nr[[i]] <- length(w)
        }
    }
    v <- unique(rev(unlist(v)))
    nr <- rev(unlist(nr))
    D2 <- matrix(0, nrow = length(v), ncol = ncol(D))
    D2[which(!is.na(v)), ] <- D[v[which(!is.na(v))], ]
    dimnames(D2) <- list(rownames(D)[match(v,rownames(D))], colnames(D))
    colorder = unlist(sccg$scc[topo.order], use.names = FALSE)
    D2 = D2[, colorder]
    cs <- cumsum(nr)[-length(nr)]
    nrcolors = 200
    half = 1 + nrcolors/2
    if (palette == "BlueRed") 
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9, 
            "OrRd")[1:9])
    else if (palette == "Grey") 
        colpal = brewer.pal(9, "Greys")[9:1]
    else stop("Unknown palette!")
    allcolors = colorRampPalette(colpal)(nrcolors)
    maxwrittenrows = 140
    if (legend) {
        nf = split.screen(rbind(c(0, 0.75, 0, 1), c(0.75, 1, 
            0, 1)))
        erase.screen(1)
    }
    if (nrow(D2) < maxwrittenrows) 
        par(las = 2, mgp = c(5.5, 1, 0), mar = c(5, 1.235, 0, 
            0), cex.lab = 1.7, cex.main = 2, lwd = 2, oma = c(8, 
            0, 0, 0))
    else par(las = 2, mgp = c(5.5, 1, 0), mar = c(0, 1.235, 0, 
        0), cex.lab = 1.7, cex.main = 2, lwd = 2, oma = c(8, 
        0, 0, 0))
    args = list(...)
    if ("zlim" %in% names(args)) 
        rangeall = args[["zlim"]]
    else {
        r = c(quantile(D2[D2 < 0], 0.95), quantile(D2[D2 > 0], 
            0.95))
        rangeall = c(-max(abs(r), na.rm = TRUE), max(abs(r), 
            na.rm = TRUE))
    }
    D2[D2 < rangeall[1]] = rangeall[1]
    D2[D2 > rangeall[2]] = rangeall[2]
    a = (length(allcolors) - 1)/(rangeall[2] - rangeall[1])
    b = 1 - a * rangeall[1]
    colors = allcolors[round(a * min(D2) + b):round(a * max(D2) + 
        b)]
    image(x = 1:nrow(D2), y = 1:ncol(D2), z = D2, xaxt = "n",
        yaxt = "n", xlab = "", ylab = "", col = colors, ...)
    if (nrow(D2) < maxwrittenrows) 
        axis(side = 1, at = 1:nrow(D2), labels = rownames(D2), 
            las = 2, cex.axis = 0.25, family="sans",tick = !is.na(v))
    axis(side = 4, at = 1:ncol(D2), labels = colnames(D2), tick = FALSE, 
        las = 1, cex.axis = 0.7,family="sans")
    axis(side = 1, at = c(0, cs) + nr/2, labels = rev(myorder), 
        tick = FALSE, cex.axis = 0.7, family="sans",las = 2, outer = TRUE)
    
    box()
    if (border) 
        abline(v = cs + 0.5, col = "black", lwd = 2)
    if (legend) {
        erase.screen(2)
        screen(2)
        par(mar = c(0, 0, 0.1, 1))
        if (palette == "BlueRed") 
            color.legend(0.5, 0.1, 1, 1, signif(seq(rangeall[1], 
                rangeall[2], length.out = 10), digits = 1), rect.col = allcolors, 
                gradient = "y", cex = 0.75)
        else if (palette == "Grey") 
            color.legend(0.5, 0.1, 1, 1, signif(seq(rangeall[1], 
                rangeall[2], length.out = 10), digits = 1), rect.col = allcolors, 
                gradient = "y", cex = 0.75)
        close.screen(c(1, 2), all.screens = TRUE)
    }
    return(v)
}

#####
plotEffects3<-
function (D, nem, border = TRUE, legend = TRUE, order = NULL, 
    orderSCC = TRUE, palette = "BlueRed", ...) 
{
    if (!(class(D) %in% c("matrix", "data.frame"))) 
        stop("First argument has to be the data matrix and second the nem object!")
    sccg <- SCCgraph(nem$graph, name = TRUE)
    if (numEdges(sccg$graph)==0)  {
         topo.order <- graph:::nodes(sccg$graph)
   } else {
   topo.order <- tsort(sccg$graph)
    }
    if (is.null(order)) 
        myorder = topo.order
    else {
        if (!orderSCC) 
            myorder = rev(unique(sccg$which.scc[match(order, 
                names(sccg$which.scc))]))
        else myorder = rev(names(sccg$scc[order]))
    }
    if (class(nem) == "score") {
        mappos = nem$mappos[[which.max(nem$mLL)]]
    }
    else {
        mappos = nem$mappos
        if (length(mappos) == 1) 
            mappos = mappos[[1]]
    }
    selected = nem$selected
    null.genes = unique(unlist(mappos["null"], use.names = FALSE))
    if (length(null.genes) > 0) {
        sccg$scc[["null"]] = "null"
        myorder = c("null", myorder)
    }
    v <- list()
    nr <- list()
    for (i in myorder) {
        w = unique(unlist(mappos[sccg$scc[[i]]], use.names = FALSE))
        if (length(w) == 0) {
            v[[i]] <- NA
            nr[[i]] <- 1
        }
        if (length(w) == 1) {
            v[[i]] <- w
            nr[[i]] <- 1
        }
        if (length(w) > 1) {
            d <- dist(D[w, , drop = FALSE], method = "manhattan")
            v[[i]] <- w[hclust(d)$order]
            nr[[i]] <- length(w)
        }
    }
    v <- unique(rev(unlist(v)))
    nr <- rev(unlist(nr))
    D2 <- matrix(0, nrow = length(v), ncol = ncol(D))
    D2[which(!is.na(v)), ] <- D[v[which(!is.na(v))], ]
    dimnames(D2) <- list(rownames(D)[match(v,rownames(D))], colnames(D))
    colorder = unlist(sccg$scc[topo.order], use.names = FALSE)
    D2 = D2[, colorder]
    cs <- cumsum(nr)[-length(nr)]
    nrcolors = 200
    half = 1 + nrcolors/2
    if (palette == "BlueRed") 
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9, 
            "OrRd")[1:9])
    else if (palette == "Grey") 
        colpal = brewer.pal(9, "Greys")[9:1]
    else stop("Unknown palette!")
    allcolors = colorRampPalette(colpal)(nrcolors)
    maxwrittenrows = 60
    if (legend) {
        nf = split.screen(rbind(c(0, 0.75, 0, 1), c(0.75, 1, 
            0, 1)))
        erase.screen(1)
    }
    if (nrow(D2) < maxwrittenrows) 
        par(las = 2, mgp = c(5.5, 1, 0), mar = c(5, 1.235, 0, 
            0), cex.lab = 1.7, cex.main = 2, lwd = 2, oma = c(8, 
            0, 0, 0))
    else par(las = 2, mgp = c(5.5, 1, 0), mar = c(0, 1.235, 0, 
        0), cex.lab = 1.7, cex.main = 2, lwd = 2, oma = c(8, 
        0, 0, 0))
    args = list(...)
    if ("zlim" %in% names(args)) 
        rangeall = args[["zlim"]]
    else {
        r = c(quantile(D2[D2 < 0], 0.95), quantile(D2[D2 > 0], 
            0.95))
        rangeall = c(-max(abs(r), na.rm = TRUE), max(abs(r), 
            na.rm = TRUE))
    }
    D2[D2 < rangeall[1]] = rangeall[1]
    D2[D2 > rangeall[2]] = rangeall[2]
    a = (length(allcolors) - 1)/(rangeall[2] - rangeall[1])
    b = 1 - a * rangeall[1]
    colors = allcolors[round(a * min(D2) + b):round(a * max(D2) + 
        b)]
    image(x = 1:nrow(D2), y = 1:ncol(D2), z = D2, xaxt = "n",
        yaxt = "n", xlab = "", ylab = "", col = colors, ...)
    if (nrow(D2) < maxwrittenrows) 
        axis(side = 1, at = 1:nrow(D2), labels = rownames(D2), 
            las = 2, cex.axis = 0.75, tick = !is.na(v))
    axis(side = 4, at = 1:ncol(D2), labels = colnames(D2), tick = FALSE, 
        las = 1, cex.axis = 0.75)
    axis(side = 1, at = c(0, cs) + nr/2, labels = rev(myorder), 
        tick = FALSE, cex.axis = 1.0, las = 2, outer = TRUE)
    
    box()
    if (border) 
        abline(v = cs + 0.5, col = "black", lwd = 2)
    if (legend) {
        erase.screen(2)
        screen(2)
        par(mar = c(0, 0, 0.1, 1))
        if (palette == "BlueRed") 
            color.legend(0.5, 0.1, 1, 1, signif(seq(rangeall[1], 
                rangeall[2], length.out = 10), digits = 1), rect.col = allcolors, 
                gradient = "y", cex = 0.75)
        else if (palette == "Grey") 
            color.legend(0.5, 0.1, 1, 1, signif(seq(rangeall[1], 
                rangeall[2], length.out = 10), digits = 1), rect.col = allcolors, 
                gradient = "y", cex = 0.75)
        close.screen(c(1, 2), all.screens = TRUE)
    }
    return(v)
}


#######
plotnem2<-
function (D, G, x, SCC, main = NULL, zlim = NULL, draw.lines = FALSE, 
    palette = "BlueRed") 
{
    if (is(D, "list")) 
        stop("Plotting of effects matrix is only implemented for one dataset, not for several.")
    if (length(x$mLL) > 1 & class(x) != "dynoNEM") {
        winner <- which.max(x$mLL)
        mappos <- x$mappos[[winner]]
    }
    else {
        mappos <- x$mappos
        if (length(mappos) == 1) 
            mappos = mappos[[1]]
    }
    nf = split.screen(rbind(c(0, 1, 0.4, 1), c(0, 0.8, 0, 0.4), 
        c(0.9, 1, 0, 0.4)))
    erase.screen(2)
    screen(2)
    mynodes = AgNode(G)
    nodenames = sapply(mynodes, name)
    xy = getNodeXY(G)
    left = xy$x <= max(xy$x) * 0.5
    right = xy$x > max(xy$x) * 0.5
    nodenames = c(nodenames[left][order(xy$x[left] + xy$y[left])], 
        nodenames[right][order(xy$x[right] + xy$y[right], decreasing = TRUE)])
    if (!is.null(zlim)) 
        ord = plotEffects2(D, x, legend = FALSE, order = nodenames, 
            orderSCC = SCC, zlim = zlim, palette = palette)
    else ord = plotEffects2(D, x, legend = FALSE, order = nodenames, 
        orderSCC = SCC, palette = palette)
    erase.screen(3)
    screen(3)
    nrcolors = 200
    half = 1 + nrcolors/2
    if (palette == "BlueRed") 
        colpal = c(brewer.pal(9, "Blues")[9:1], brewer.pal(9, 
            "OrRd")[1:9])
    else if (palette == "Grey") 
        colpal = brewer.pal(9, "Greys")
    else stop("Unknown palette!")
    allcolors = colorRampPalette(colpal)(nrcolors)[101:200]
    if (is.null(zlim)) {
        r = c(quantile(D[D < 0], 0.95), quantile(D[D > 0], 0.95))
    ##rangeall = c(-max(abs(r), na.rm = TRUE), max(abs(r),na.rm = TRUE))
        rangeall = c(0,1)
    }
    else rangeall = zlim
    color.legend(0.4, 0.1, 1, 1, signif(seq(rangeall[1], rangeall[2], 
        length.out = 5), digits = 1), rect.col = allcolors, gradient = "y", 
        cex = 0.75)
    erase.screen(1)
    screen(1)
    par(mar = c(0, 0, 0, 0), cex = 1, cex.main = 2)
    plot(G, main = main)
    if (draw.lines) {
        c = 0.7912
        ma = c * par()$usr[2]
        mi = getX(botLeft(boundBox(G)))
        may = getY(upRight(boundBox(G)))
        miy = min(xy$y)
        if (miy - max(getNodeHeight(G)) > 0) 
            miy = 0
        xrescale = function(x, a2 = mi, b2 = ma) {
            alpha = (b2 - a2)/length(ord)
            beta = mi
            alpha * x + beta
        }
        for (i in 1:length(mynodes)) {
            xy = getNodeXY(mynodes[[i]])
            nodenames = strsplit(name(mynodes[[i]]), ":")[[1]]
            to = xrescale(match(unique(unlist(mappos[nodenames])), 
                ord) - 0.5)
            to = sort(to)
            for (j in 1:length(to)) {
                segments(xy$x, max(miy, xy$y - 7), to[j], miy, 
                  col = "grey", lty = 1, lwd = 1, type = "b", 
                  pch = 19, cex = 0.3)
            }
        }
        screen(1)
        plot(G, main = main)
    }
    close.screen(c(1, 2, 3), all.screens = TRUE)
}

