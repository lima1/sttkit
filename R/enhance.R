#' enhance_bayesspace
#'
#' Enhance BayesSpace data
#' @param x Object, converted by \code{\link{as_SingleCellExperiment}}
#' @param test_num_clusters Range of number of clusters to be tested (\code{qs}). 
#' @param num_clusters Optional picked number
#' @param force Recalculate, even when serialized objects are available
#' @param serialize Serialize output objects
#' @param prefix Prefix of output files
#' @param ... Additional paramters passed to \code{BayesSpace::spatialCluster}
#' @export enhance_bayesspace
#' @examples
#' #enhance_bayesspace()

enhance_bayesspace <- function(x, test_num_clusters = seq(2, 12),
                               num_clusters = NULL, ...) {
    flog.info("Enhancing %i clusters. This will take a while...", opt$num_clusters)
    ndata_enhanced <- spatialEnhance(ndata,
        q = opt$num_clusters,
        nrep = opt$num_iter,
        verbose = TRUE)
    flog.info("Writing R data structure to %s...", filename)
    saveRDS(ndata_enhanced, file = filename)
}    


#' enhance_deconvolve
#'
#' Enhances deconvoluted Spatial data
#' @param object Seurat object
#' @param num_subspots Split spot into n subspots
#' @param max_cell_types Too avoid plotting too many cell types, we can limit them here.
#' Set to -1 to turn off.
#' @param num_iter Number optimization iterations
#' @param stroke \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param crop \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param pt.size.factor \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param image.alpha \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param ... Additional parameters passed to \code{\link{Seurat::SingleSpatialPlot}}.
#' @export enhance_deconvolve
#' @examples
#' #enhance_deconvolve()

enhance_deconvolve <- function(object, num_subspots = 6, max_cell_types = 7,
    num_iter = 10, stroke = NA, crop = FALSE, pt.size.factor = 0.4, image.alpha = 1, ...) { 
    Y <- as.matrix(GetAssayData(object, assay = "predictions"))
    Y <- t(Y[which(rownames(Y) != "max"), ])
    if (max_cell_types > 0 && nrow(Y) > max_cell_types) {
        features <- names(head(sort(apply(Y, 2, function(x) sum(x > 1 / num_subspots)),
            decreasing = TRUE), max_cell_types))
        other <- 1 - apply(Y[,features], 1, sum) 
        Y <- cbind(Y[, features], "other" = other[rownames(Y)])
        features <- c(features, "other")
    } else {
        features <- names(sort(apply(Y, 2, function(x) sum(x > 1 / num_subspots)),
            decreasing = TRUE))
    }     
    d <- ncol(Y)
    n0 <- nrow(Y)

    positions <- as.matrix(GetTissueCoordinates(object))[, 2:1]
    colnames(positions) <- c("x", "y")
    Y2 <- Y[rep(seq_len(n0), num_subspots), ]
    rows_cols <- object@images[[Images(object)]]@coordinates[, 2:3]
    positions <- cbind(positions, rows_cols[rownames(positions), ])
    positions2 <- positions[rep(seq_len(n0), num_subspots), ]
    shift <- BayesSpace:::.make_subspot_offsets(num_subspots)
    xdist <- coef(lm(x ~ row, data = positions))[2]
    ydist <- coef(lm(y ~ col, data = positions))[2]
    shift <- t(t(shift) * c(xdist, ydist))
    dist <- max(rowSums(abs(shift))) * 1.05
    shift_long <- shift[rep(seq_len(num_subspots), each = n0), ]
    positions2[, "x"] <- positions2[, "x"] + shift_long[, "Var1"]
    positions2[, "y"] <- positions2[, "y"] + shift_long[, "Var2"]
    n <- nrow(Y2)
    df_j <- BayesSpace:::find_neighbors(positions2, dist, "manhattan")
    colnames(positions2)[1:4] <- c("imagecol", "imagerow", "spot.row", "spot.col")
    n_spots <- ncol(object)
    n_subspots <- nrow(positions2)
    idxs <- seq_len(n_subspots)
    spot_idxs <- ((idxs - 1) %% n_spots) + 1
    subspot_idxs <- rep(seq_len(num_subspots), each = n_spots)
    positions2 <- as.data.frame(positions2)
    rownames(positions2) <- gsub("\\.", "-", rownames(positions2))
    positions2$spot.idx <- spot_idxs
    positions2$subspot.idx <- subspot_idxs
    ct <- apply(Y,1,function(x) { 
        r1 <- ceiling(x * num_subspots)
        r <- round(x * num_subspots)
        r[r < 1] <- r1[r < 1]
        i <- 1
        while(sum(r) != num_subspots) {
            if (sum(r) < num_subspots) r[order(x, decreasing = TRUE)[i]] <- r[order(x, decreasing = TRUE)[i]] + 1
            if (sum(r) > num_subspots) r[order(x, decreasing = FALSE)[i]] <- r[order(x, decreasing = FALSE)[i]] - 1
            r[r < 0] <- 0
            i <- i + 1
        }
        return(r)
    })
    positions2$cell.type <- NA
    for (i in seq(n_spots)) {
        cts <- unlist(sapply(seq_along(ct[, i]), function(j) if (ct[j,i]) return(rep(j, ct[j,i])) else(return(NULL))))
        positions2$cell.type[positions2$spot.idx == i] <- cts
    }
    positions2$cell.type <- colnames(Y)[positions2$cell.type]
    if (num_iter > 0) {
        positions2 <- .optimize_subspot_celltypes(positions2, df_j, num_iter = num_iter)
    }
    vertices <- .make_triangle_subspots(positions2, fill="cell.type")
    
    vertices$imagecol <- positions2[vertices$spot, "imagecol"]
    vertices$imagerow <- positions2[vertices$spot, "imagerow"]
    fit <- lm(imagerow~x.vertex, data = vertices)
    vertices$y.image.vertex <- predict(fit, vertices)
    fit <- lm(imagecol~y.vertex, data = vertices)
    vertices$x.image.vertex <- predict(fit, vertices)
    vertices$fill <- factor(vertices$fill, levels = features)
    datax <- vertices[, c("y.image.vertex", "x.image.vertex", "fill")]
    colnames(datax) <- c("imagerow", "imagecol", "Cell Type")
    splot <- SingleSpatialPlot(data = datax, image = object@images[[Images(object)]],
        pt.size.factor = pt.size.factor, col.by = "Cell Type", stroke = stroke, crop = crop, ...)
    if (requireNamespace("ggthemes", quietly = TRUE) &&
            length(levels(vertices$fill)) <= 8) {
        splot <- splot +  ggthemes::scale_fill_colorblind()
    }
    return(splot)
}

.optimize_subspot_celltypes <- function(positions2, df_j, num_iter) {
    positions2_orig <- positions2
    n_total <- 0
    n_total_best <- .Machine$integer.max

    for (iter in seq(0, num_iter)) {
        table(positions2$cell.type)
        n_total_before <- n_total
        n_total <- 0
        for (i in seq_along(df_j)) {
            pss1 <- positions2[df_j[[i]], ]
            n_cts_before <- length(table(pss1$cell.type))
            n_total <- n_total + n_cts_before
            if (!iter) next
            pss2 <- positions2[positions2$spot.idx %in% pss1$spot.idx, ]
            pss2 <- pss2[!rownames(pss2) %in% rownames(pss1),]
            for (j in unique(pss1$spot.idx)) {
                pss2x <- pss2[pss2$spot.idx == j, ]
                pss2x <- pss2x[pss2x$cell.type %in% pss1$cell.type, ]
                if (nrow(pss2x)) {
                    pss1x <- pss1[pss1$spot.idx == j, ]
                    r1 <- sample(rownames(pss1x), 1)
                    r2 <- sample(rownames(pss2x), 1)
                    c1 <- positions2[r1, "cell.type"]
                    c2 <- positions2[r2, "cell.type"]
                    positions2[r1, "cell.type"] <- c2
                    positions2[r2, "cell.type"] <- c1
                    if (length(table(positions2[df_j[[i]], "cell.type"])) >= n_cts_before) {
                        positions2[r1, "cell.type"] <- c1
                        positions2[r2, "cell.type"] <- c2
                    }
                }
            }
        }
        flog.debug("Iteration %i: %i (best %i)", iter, n_total, n_total_best)
        if (n_total < n_total_best) {
            positions2_best <- positions2
            n_total_best <- n_total
        }
    }
    return(positions2_best)
}

.make_triangle_subspots <- function (cdata, fill = "spatial.cluster") {
    spot_positions <- BayesSpace:::.select_subspot_positions(cdata, x = "spot.col", 
        y = "spot.row", fill = fill)
    spot_positions <- BayesSpace:::.adjust_hex_centers(spot_positions)
    r <- 1/3
    R <- (2/sqrt(3)) * r
    vertex_offsets <- do.call(rbind, list(data.frame(x.offset = c(0, 
        0, r), y.offset = c(0, -R, -R/2), subspot.idx = 3), data.frame(x.offset = c(0, 
        r, r), y.offset = c(0, -R/2, R/2), subspot.idx = 5), 
        data.frame(x.offset = c(0, r, 0), y.offset = c(0, R/2, 
            R), subspot.idx = 1), data.frame(x.offset = c(0, 
            0, -r), y.offset = c(0, R, R/2), subspot.idx = 2), 
        data.frame(x.offset = c(0, -r, -r), y.offset = c(0, R/2, 
            -R/2), subspot.idx = 6), data.frame(x.offset = c(0, 
            -r, 0), y.offset = c(0, -R/2, -R), subspot.idx = 4)))
    spot_vertices <- BayesSpace:::.make_spot_vertices(spot_positions, vertex_offsets)
    spot_vertices$y.vertex <- -spot_vertices$y.vertex
    spot_vertices
}

