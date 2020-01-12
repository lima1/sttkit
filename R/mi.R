# a couple of helper functions that calculate the correlation of features
# by taking the 4 neighors into account where available

.find_nn <- function(obj) {
    xy <- .parse_coords(obj, colnames(obj))
    xy[,2] <- 1- xy[,2]
    d <- dist(xy, method="manhattan")
    dm <- as.matrix(d)
    nn <- lapply(seq(nrow(dm)), function(i) as.vector(which(dm[i,] > 0 & dm[i, ] < 1.1)))
    nn <- do.call(rbind, lapply(nn, function(x) c(x, rep(NA, 4-length(x)))))
}
.cor_nn_vector <- function(obj, nn, f1, f2, zero_offset, method, slot = "data") {
    nn <- cbind(1:nrow(nn), nn)
    weights <- apply(nn, 2, function(x) 1-sum(is.na(x))/length(x))
    idx <- weights > 0
    nn <- nn[,idx]
    weights <- weights[idx]
    weighted.mean(apply(nn, 2, function(i) {
        p <- .fix_offset_in_pair(
            FetchData(obj, vars = f1, slot = slot)[,1], 
            FetchData(obj, vars = f2, slot = slot)[,1], 
            zero_offset)
        cor(p[,1], p[,2], use="complete.obs", method = method)
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
