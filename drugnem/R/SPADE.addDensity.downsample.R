#' Function to perform density downsampling of a given FCS file
#'
#' @param infilename String of FCS file name that should be used as input
#' @param cols vector of column ids of lineage markers of interest to be used by ccast algorithm to define cell clusters.
#' @param arcsinh_cofactor Cofactor used in arcsinh transform asinh(data/arcsinh_cofactor) of data
#' @param kernel_mult Multiplier of the minimum median distance within which other observations are counted towards the density
#' @param apprx_mult Multiplier of the minimum median distance within which observations are approximated to have the same density
#' @param med_samples Number of observations used to estimate the minimum median distance
#' @param comp Apply compensation matrix if present in SPILL or SPILLOVER keywords
#' @param exclude_pctile Numeric value in [0,1]. Densities below this percentile will be excluded.
#' @param target_pctile Numeric value in [0,1]. Densities below this percentile, but above exclude_pctile will be retained. Only meaningful if desired_samples is NULL.
#' @param desired_samples Desired number of samples. If set to integer value, the target percentile will be set internally to downsample to approximately the desired number of samples.
#' @return The name of the written file is returned in addition to a dataframe of subsampled cells with their marker expression values.
#' @seealso \code{\link{fcsimport3}} which this function wraps
#' @export
#' @examples
#' 

##### Density downsampling ##############
SPADE.addDensity.downsample<-
function (infilename, cols = colid, arcsinh_cofactor = asinhp,
kernel_mult = 5, apprx_mult = 1.5, med_samples = 2000, comp = TRUE,exclude_pctile = 0.01, target_pctile = 0.05, desired_samples = subsamplesize)
{
	###require(flowCore)
    ##in_fcs <- spade:::SPADE.read.FCS(infilename, comp = comp)
    in_fcs <- read.FCS(infilename, transformation = FALSE)
    in_data1 <- exprs(in_fcs)
    rownames(in_data1)=paste(1:dim(in_data1)[1])
    in_data <- in_data1[,cols]
    if(desired_samples>dim(in_data1)[1]) {
        out_downsample2=in_data1
        params <- flowCore::parameters(in_fcs)
        desc <- description(in_fcs)
        pd <- pData(params)
        
        firstout <- paste(basename(infilename), ".downsample.fcs", sep="")
        out_frame <- flowFrame(out_downsample2, params, description = desc)
        
        write.FCS(out_frame,firstout)

    } else {
    ########Estimate density#############
    density <- spade:::SPADE.density(asinh(in_data1/arcsinh_cofactor),
    kernel_mult = kernel_mult, apprx_mult = apprx_mult, med_samples = med_samples)
    if (max(density) == 0)
    warning(paste(infilename, "has degenerate densities, possibly due to many identical observations",
    sep = " "))
    in_data <- cbind(in_data, density = density)
    ########## downsample#######
    d_idx <- match("density", colnames(in_data))
    if (is.na(d_idx)) {
        stop("No density parameter in FCS file")
    }
    boundary <- quantile(in_data[, d_idx], c(exclude_pctile,
    target_pctile), names = FALSE)
    out_data <- subset(in_data, in_data[, d_idx] > boundary[1])
    density <- out_data[, d_idx]
    if (is.null(desired_samples)) {
        boundary <- boundary[2]
        out_data <- subset(out_data, boundary/density > runif(nrow(out_data)))
    }
    else if (desired_samples < nrow(out_data)) {
        density_s <- sort(density)
        cdf <- rev(cumsum(1/rev(density_s)))
        boundary <- desired_samples/cdf[1]
        if (boundary > density_s[1]) {
            targets <- (desired_samples - 1:length(density_s))/cdf
            boundary <- targets[which.min(targets - density_s >
            0)]
        }
        
        out_data <- subset(out_data, boundary/density > runif(length(density)))
    }
    if(dim(out_data)[1]==0) {
        out_downsample=in_data1[sample(1:dim(in_data1)[1],desired_samples,replace = FALSE),]
    } else {
    out_downsample=out_data[,-dim(out_data)[2]]
    }
    out_downsample2=in_data1[as.numeric(rownames(out_downsample)),]
    params <- flowCore::parameters(in_fcs)
    desc <- description(in_fcs)
    pd <- pData(params)
        
    firstout <- paste(basename(infilename), ".downsample.fcs", sep="")
    out_frame <- flowFrame(out_downsample2, params, description = desc)
        
    write.FCS(out_frame,firstout)
    }
    ##out_downsample=out_data
    return(out_downsample2)
}

