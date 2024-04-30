#' enhance_bayesspace
#'
#' Enhance BayesSpace data
#' @param x Object, converted by \code{\link{as_SingleCellExperiment}}
#' @param test_num_clusters Range of number of clusters to be tested (\code{qs}). 
#' @param num_clusters Optional picked number
#' @param ... Additional paramters passed to \code{BayesSpace::spatialEnhance}
#' @export enhance_bayesspace
#' @examples
#' #enhance_bayesspace()

enhance_bayesspace <- function(x, test_num_clusters = seq(2, 12),
                               num_clusters = NULL, ...) {
    flog.info("Enhancing %i clusters. This will take a while...", opt$num_clusters)
    ndata_enhanced <- spatialEnhance(ndata,
        q = opt$num_clusters,
        nrep = opt$num_iter,
        verbose = TRUE, ...)
    flog.info("Writing R data structure to %s...", filename)
    saveRDS(ndata_enhanced, file = filename)
}    


#' enhance_deconvolve
#'
#' Enhances deconvoluted Spatial data
#' @param object Seurat object
#' @param num_subspots Split spot into n subspots
#' @param max_cell_types Too avoid plotting too many cell types, we can limit them here.
#' Set to -1 to turn off. Ordered by prevalence if \code{cell_types} is \code{NULL}.
#' @param cell_types Only plot the specified cell types. 
#' @param image Image to plot
#' @param num_iter Number optimization iterations
#' @param stroke \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param crop \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param pt.size.factor \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param pt.alpha \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param image.alpha \code{\link{Seurat::SingleSpatialPlot}} argument
#' @param verbose Verbose output
#' @importFrom stats coef
#' @export enhance_deconvolve
#' @examples
#' #enhance_deconvolve()

enhance_deconvolve <- function(object, num_subspots = 6, max_cell_types = 8, cell_types = NULL,
    image = NULL, num_iter = 10, stroke = NA, crop = FALSE, pt.size.factor = 0.4, pt.alpha = NULL,
    image.alpha = 1, verbose = FALSE) { 
    image <- image %||% Images(object = object)
    if (length(x = image) == 0) {
        image <- Images(object = object)
    }
    if (length(x = image) < 1) {
        stop("Could not find any spatial image information")
    }
    image.use <- object[[image]]

    Y <- as.matrix(GetAssayData(object, assay = "predictions"))
    Y <- t(Y[which(rownames(Y) != "max"), ])
    if (!is.null(cell_types)) {
        cell_types <- cell_types[cell_types %in% colnames(Y)]
    }
    if (max_cell_types > 0 && ncol(Y) > max_cell_types) {
        if (!is.null(cell_types)) {
            features <- head(cell_types, max_cell_types)
        } else {
            features <- names(head(sort(apply(Y, 2, function(x) sum(x > 1 / num_subspots)),
                decreasing = TRUE), max_cell_types))
        }
        other <- 1 - apply(Y[,features], 1, sum) 
        Y <- cbind(Y[, features], "other" = other[rownames(Y)])
        features <- c(features, "other")
    } else {
        if (!is.null(cell_types)) {
            features <- cell_types
        } else {
            features <- names(sort(apply(Y, 2, function(x) sum(x > 1 / num_subspots)),
                decreasing = TRUE))
        }
        y <- Y[, features]    
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
        positions2 <- .optimize_subspot_celltypes(positions2, df_j, num_iter = num_iter, verbose = verbose)
    }
    vertices <- .make_triangle_subspots(positions2, fill = "cell.type")
    vertices$imagecol <- positions2[vertices$spot, "imagecol"]
    vertices$imagerow <- positions2[vertices$spot, "imagerow"]
    fit <- lm(imagerow~x.vertex, data = vertices)
    vertices$y.image.vertex <- predict(fit, vertices)
    fit <- lm(imagecol~y.vertex, data = vertices)
    vertices$x.image.vertex <- predict(fit, vertices)
    vertices$fill <- factor(vertices$fill, levels = features)
    vertices$spot <- as.factor(vertices$spot)
    data <- vertices[, c("y.image.vertex", "x.image.vertex", "fill", "spot")]
    splot <- ggplot(data = data, aes_string(x = colnames(x = data)[2], 
          y = colnames(x = data)[1],
          fill = colnames(x = data)[3],
          group = colnames(x = data)[4]
          ))

    if (is.null(x = pt.alpha)) {
        splot <- splot + geom_spatial_polygon(point.size.factor = pt.size.factor, 
            data = data, image = image.use, image.alpha = image.alpha, 
            crop = crop, stroke = stroke, )
    } else {
        splot <- splot + geom_spatial_polygon(point.size.factor = pt.size.factor, 
            data = data, image = image.use, image.alpha = image.alpha, 
            crop = crop, stroke = stroke, alpha = pt.alpha)
    }
    splot <- splot + coord_fixed() + theme(aspect.ratio = 1) + labs(fill = "Cell Type", x = "", y = "") +
        guides(fill = guide_legend(override.aes = list(size = 3)))

    palette <- .get_colorblind_pal(length(levels(vertices$fill)))
    if (!is.null(palette)) {
        splot <- splot + scale_fill_manual(values = palette)
    }
    splot <- splot + NoAxes() + theme(panel.background = element_blank())
    return(splot)
}

.get_colorblind_pal <- function(n) {
    palette <- NULL
    if (requireNamespace("ggthemes", quietly = TRUE) && n <= 9) {
        # black color last, not first and use gray instead
        palette <- ggthemes::colorblind_pal()(min(8, n))
        palette[1] <- "#999999"
        palette <- c(palette, "#014f3a")
        palette <- head(palette[c(seq(2, length(palette)), 1)], n)
    }
    return(palette)    
}    
.optimize_subspot_celltypes <- function(positions2, df_j, num_iter, verbose = FALSE) {
    positions2_orig <- positions2
    n_total <- 0
    n_total_best <- .Machine$integer.max

    for (iter in seq(0, num_iter)) {
        if (verbose) print(table(positions2$cell.type))
        n_total <- 0
        for (i in seq_along(df_j)) {
            pss1 <- positions2[df_j[[i]], ]
            n_cts_before <- length(table(pss1$cell.type))
            n_total <- n_total + n_cts_before
            if (!iter) next
            pss2 <- positions2[positions2$spot.idx %in% pss1$spot.idx, ]
            pss2$nn <- rownames(pss2) %in% rownames(pss1)
            pss2x <- split(pss2, pss2$spot.idx)
            if (length(pss2x) != 2) {
            #    flog.info("hhh %i %i iter %i", i, length(pss2x), iter)
                next
            }
            r <- lapply(1:2, function(j) replicate(30, {
                     sample(pss2x[[j]]$cell.type, nrow(pss2x[[j]]), replace = FALSE)
            }))
            n_cts_after <- sapply(seq(30), function(i) length(table(c(r[[1]][pss2x[[1]]$nn, i], r[[2]][pss2x[[2]]$nn, i]))))
            idx <- which.min(n_cts_after)
            if (n_cts_after[idx] < n_cts_before) {
                pss2x[[1]]$cell.type <- r[[1]][, idx]
                pss2x[[2]]$cell.type <- r[[2]][, idx]
                pss2 <- rbind(pss2x[[1]], pss2x[[2]])
                positions2[rownames(pss2),"cell.type"] <- pss2$cell.type
            }
        }
        if (verbose) flog.info("Iteration %i: %i (best %i)", iter, n_total, n_total_best)
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
    r <- 0.5
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

# For plotting the tissue image
#' @importFrom ggplot2 ggproto Geom aes ggproto_parent alpha draw_key_point
#' @importFrom grid unit gpar editGrob polygonGrob viewport gTree addGrob grobName
#'
GeomSpatialPolygon <- ggproto(
  "GeomSpatialPolygon",
  Geom,
  required_aes = c("x", "y", "group"),
  extra_params = c("na.rm", "image", "image.alpha", "crop"),
  default_aes = aes(
    shape = 21,
    colour = NA,
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = NA
  ),
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    # We need to flip the image as the Y coordinates are reversed
    data$y = max(data$y) - data$y + min(data$y)
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(data, panel_scales, coord, image, image.alpha, crop) {
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    if (!crop) {
      y.transform <- c(0, nrow(x = image)) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, ncol(x = image))
      panel_scales$y$continuous_range <- c(0, nrow(x = image))
      panel_scales$y.range <- c(0, nrow(x = image))
      panel_scales$x.range <- c(0, ncol(x = image))
    }
    z <- coord$transform(
      data.frame(x = c(0, ncol(x = image)), y = c(0, nrow(x = image))),
      panel_scales
    )
    # Flip Y axis for image
    z$y <- -rev(z$y) + 1
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    img.grob <- GetImage(object = image)

    img <- editGrob(grob = img.grob, vp = vp)
    coords <- coord$transform(data, panel_scales)
    firsts <- seq(1, nrow(coords), by = 3)
    firsts <- coords[firsts,]
    pts <- polygonGrob(
      x = coords$x,
      y = coords$y,
      id = rep(seq_len(nrow(firsts)), each = 3),
      default.units = "native",
      gp = gpar(
        col = firsts$colour,
        fill = scales::alpha(firsts$fill, firsts$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    if (image.alpha > 0) {
      if (image.alpha != 1) {
        img$raster = as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(x = img$raster),
            ncol = ncol(x = img$raster),
            byrow = TRUE)
        )
      }
      gt <- addGrob(gTree = gt, child = img)
    }
    gt <- addGrob(gTree = gt, child = pts)
    # Replacement for ggname
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial_polygon')
    return(gt)
  }
)

# influenced by: https://stackoverflow.com/questions/49475201/adding-tables-to-ggplot2-with-facet-wrap-in-r
# https://ggplot2.tidyverse.org/articles/extending-ggplot2.html
#' @importFrom ggplot2 layer
#'
#'
geom_spatial_polygon <-  function(
  mapping = NULL,
  data = NULL,
  image = image,
  image.alpha = image.alpha,
  crop = crop,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
) {
  layer(
    geom = GeomSpatialPolygon,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, image.alpha = image.alpha, crop = crop, ...)
  )
}
