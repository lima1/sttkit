#' deconvolute_spatial
#'
#' Simulate SpatialTranscriptomics data
#' @param obj_spatial Seurat object containing SpatialTranscriptomics data
#' @param sig Signature for deconvolution
#' @param write Write output data to text and PDF files
#' @param prefix Prefix of output files
#' @param ... Additional parameters passed to \code{\link{plot_spots}}
#' @export deconvolute_spatial
#' @examples
#' deconvolute_spatial()
deconvolute_spatial <- function(obj_spatial, sig, write = TRUE, prefix, ...) {

    if (!requireNamespace("CellMix", quietly = TRUE)) {
        stop("This function requires the CellMix package.")
    }

    spatial_scaled <- as.matrix(GetAssayData(obj_spatial, slot = "scale.data"))
    sig <- sig[rownames(sig) %in% rownames(spatial_scaled),]
    zz <- CellMix::ged(spatial_scaled, sig, log = FALSE)
    z <- CellMix::coef(zz)
    coords <- .parse_coords(obj, colnames(zz))
    d.f <- do.call(rbind, lapply(seq(nrow(z)), function(i) 
        data.frame( 
            x = coords[,1], 
            y = 1-coords[,2], 
            fraction = z[i,],
            cluster = colnames(sig)[i]
    )))
    if (write) {
        flog.info("Plotting clusters...") 
        ratio <- .get_image_ratio(nrow(z))
        width <- 10
             
        pdf(paste0(prefix, "_deconvolution_overview.pdf"), 
            width = width, height = width * ratio)
        plot_spots(d.f, ... ) 
        dev.off()
        write.csv(z, paste0(prefix, "_deconvolution.csv"))
    }
    z
}


#' find_nn_spatial
#'
#' Simulate SpatialTranscriptomics data
#' @param obj_spatial Seurat object containing SpatialTranscriptomics data
#' @param obj_ref Seurat object containing reference data
#' @param k k-nearest neighbors
#' @param write Write output data to text and PDF files
#' @param prefix Prefix of output files
#' @param ... Additional parameters passed to \code{\link{plot_features}}
#' @export find_nn_spatial
#' @importFrom RANN nn2
#' @examples
#' find_nn_spatial()
find_nn_spatial <- function(obj_spatial, obj_ref, k = 20, write = TRUE, prefix, ...) {
    nn <- nn2(data = t(as.matrix(GetAssayData(obj_spatial))), 
              query = t(as.matrix(GetAssayData(obj_ref))),
              k = k)
    counts <- matrix(0, nrow=length(levels(obj_ref)), ncol(obj_spatial))
    rownames(counts) <- levels(obj_ref)
    colnames(counts) <- colnames(obj_spatial)

    for (i in seq(ncol(counts))) {
        for (j in seq(ncol(nn$nn.idx))) {
            score <- ncol(nn$nn.idx) - j + 1
            idx <- which(nn$nn.idx[,j]==i)
            if (length(idx)) {
                counts[as.character(Idents(obj_ref))[idx],i] <- counts[as.character(Idents(obj_ref))[idx],i] + score
            }
        }    
    }    
    counts_norm <- t(t(counts)/apply(counts,2,sum))
    counts_norm[is.nan(counts_norm)] <- NA
    if (write) {
        flog.info("Plotting nearest neigbors...") 
        ratio <- .get_image_ratio(nrow(counts_norm))
        width <- 10
        pdf(paste0(prefix, "_nearest_neighbors_overview.pdf"), 
            width = width, height = width * ratio)
        plot_features(t(counts_norm), features = rownames(counts_norm), ...)
        dev.off()
        write.csv(counts_norm, paste0(prefix, "_nearest_neighbors.csv"))
    }
}    
