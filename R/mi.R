#' Cluster neighborhoods 
#'
#' This function provides the neighborhood of clusters
#'
#' @param obj Seurat Object
#' @param max.dist Defines neighborhood as maximum distance in number
#' of spots
#' @param fun.aggregate Aggregate cluster value of spot values
#' @return output pairwise cluster neighborhood overlaps
#'
#' @examples
#'
#' @importFrom stats weighted.mean dist
#' @export find_cluster_neighborhoods
find_cluster_neighborhoods <- function(obj, max.dist = 3, fun.aggregate = mean) {
    d <- dist(obj@images[[1]]@coordinates[, c("row", "col")])
    m <- as.matrix(d)
    barcodes <- colnames(m)
    .get_cluster_mean_dists <- function(labels) {
        dist_stats <- t(sapply(barcodes, function(i) {
            x <- m[i,]
            xs <- split(x, labels)
            sapply(xs, function(y) sum(y > 0 & y < max.dist))
        }))
        dd <- lapply(levels(labels), function(i) dist_stats[which(labels == i),])
        sapply(dd, apply, 2, fun.aggregate)
    }
    labels <- Idents(obj[,barcodes])
    dd <- .get_cluster_mean_dists(labels)
    rownames(dd) <- levels(labels)
    colnames(dd) <- levels(labels)
    return(dd)
}
# a couple of helper functions that calculate the correlation of features
# by taking the 4 neighors into account where available

.find_nn <- function(obj) {
    xy <- .parse_coords(obj, colnames(obj))
    xy[, 2] <- 1 - xy[, 2]
    d <- dist(xy, method = "manhattan")
    dm <- as.matrix(d)
    nn <- lapply(seq(nrow(dm)), function(i) as.vector(which(dm[i, ] > 0 & dm[i, ] < 1.1)))
    nn <- do.call(rbind, lapply(nn, function(x) c(x, rep(NA, 4 - length(x)))))
}
.cor_nn_vector <- function(obj, nn, f1, f2, zero_offset, method, slot = "data") {
    nn <- cbind(1:nrow(nn), nn)
    weights <- apply(nn, 2, function(x) 1-sum(is.na(x)) / length(x))
    idx <- weights > 0
    nn <- nn[, idx]
    weights <- weights[idx]
    weighted.mean(apply(nn, 2, function(i) {
        p <- .fix_offset_in_pair(
            FetchData(obj, vars = f1, slot = slot)[, 1],
            FetchData(obj, vars = f2, slot = slot)[, 1],
            zero_offset)
        cor(p[, 1], p[, 2], use = "complete.obs", method = method)
        }),
        w = weights)
}
.fix_offset_in_pair <- function(x, y, zero_offset = NULL) {
    if (is.null(zero_offset)) return(cbind(x = x, y = y))
    min_offset <- zero_offset / 2
    m <- cbind(x = x, y = y)
    m[m < min_offset] <- m[m < min_offset] - zero_offset
    m
} 
.cor_vector <- function(obj, f1, f2, zero_offset, method) {
    p <- .fix_offset_in_pair(obj[[f1]], obj[[f2]], zero_offset)
    cor(p[,1], p[,2], use="complete.obs", method = method)
}
.cor_nn <- function(obj, features, method = "spearman", average_nn = TRUE, zero_offset = NULL, slot = "data") {
    if (!average_nn) {
        if (is.null(zero_offset)) return(cor(FetchData(obj, vars = features, slot = slot)))
        cor_nn <- sapply(seq_along(features),function(i) sapply(seq_along(features), function(j) .cor_vector(obj,features[i], features[j], zero_offset, method)))
    } else {
        nn <- .find_nn(obj)
        cor_nn <- sapply(seq_along(features),function(i) sapply(seq_along(features), function(j) .cor_nn_vector(obj, nn, features[i], features[j], zero_offset, method, slot = slot)))
    }
    rownames(cor_nn) <- features
    colnames(cor_nn) <- features
    cor_nn
}    


#' Spatial Correlation
#'
#' This function provides the spatial correlation of two features. If both are the
#' same, this corresponds to Moran's I.
#'
#' @param obj Seurat Object
#' @param image Passed to the GetTissueCoordinates function
#' @param slot Passed to the FetchData function
#' @param feature.1 First feature
#' @param feature.2 Second feature(s)
#' @return spatial correlation
#'
#' @examples
#'
#' @export calculate_spatial_correlation
calculate_spatial_correlation <- function(obj, image = NULL, slot = "data",
    feature.1, feature.2 = rownames(obj)) {

    coord <- GetTissueCoordinates(obj, image = image)
    pos.dist <- dist(coord)
    pos.dist.mat <- as.matrix(x = pos.dist)
    weights <- 1 / pos.dist.mat^2
    diag(x = weights) <- 0
    if(length(feature.1) > 1) {
        flog.warn("feature.1 contains multiple values, using only the first")
        feature.1 <- feature.1[1]
    }
    x_values <- as.numeric(FetchData(obj, feature.1, slot = slot)[rownames(coord), 1])
    x_values <- (x_values - mean(x_values)) / sqrt(var(x_values))
    y_values_all <- FetchData(obj, feature.2, slot = slot)

    .spatial_corr <- function(gene) {
        y_values <- as.numeric(y_values_all[rownames(coord), gene])
        y_values <- (y_values - mean(y_values)) / sqrt(var(y_values))
        corr <- sum(outer(x_values, y_values) * weights) / sum(weights)
        return(corr)
    }

    results <- sapply(feature.2, .spatial_corr)
    crosscorr <- as.numeric(results)
    names(crosscorr) <- feature.2
    crosscorr <- sort(crosscorr, decreasing = TRUE)
    return(crosscorr)
}

.get_binsize <- function(obj, image = NULL, max_bins = 5000) {
    if (ncol(obj) <= max_bins) return(list(x.cuts = NULL, y.cuts = NULL))
    tc <- GetTissueCoordinates(obj, image = image)
    r <- (max(tc[,1], na.rm = TRUE) - min(tc[,1], na.rm = TRUE)) / (max(tc[,2], na.rm = TRUE) - min(tc[,2], na.rm = TRUE))
    xy.cuts <- ceiling(sqrt(max_bins))
    if (r > 1) return(list(x.cuts = xy.cuts, y.cuts = round(xy.cuts / r)))
    return(list(x.cuts = round(xy.cuts) / r, y.cuts = xy.cuts))
}    
