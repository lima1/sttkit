#' plot_features
#'
#' Plots features in a Seurat object on HE slides
#' @param obj Seurat object
#' @param features \code{character(n)} of features to be plotted
#' @param cells Plot only specified cells
#' @param undetected_NA Remove spots with low value by setting them to \code{NA}
#' @param plot_correlations Plot pairwise correlations of feature scores
#' @param plot_violin Plot distribution of feature scores across idents
#' @param plot_map Plots features in \code{facet_wrap} for colorblind creatures.
#' Requires that features are \code{factors}, like cluster labels.
#' @param zero_offset If zero values were set to a low value, specify the offset here.
#' @param slot Slot to pull feature data for
#' @param ... Arguments passed to \code{\link{plot_spots}}.
#' @export plot_features
#' @import ggplot2
#' @importFrom stats quantile
#' @examples
#' plot_features()

plot_features <- function(obj, features, cells = NULL, undetected_NA = FALSE, 
                          plot_correlations = TRUE, plot_violin = TRUE, plot_map = FALSE, 
                          zero_offset = NULL, slot = "data", ...)  {
    if (is.null(cells)) cells <- colnames(obj)
    limits <- NULL 
    if (is(obj, "Seurat") && 0 && length(Images(obj))) {
        print(SpatialFeaturePlot(obj[, cells], features, slot = slot))
    } else {    
        if (is(obj, "Seurat")) {
            z <- FetchData(obj, vars = features, slot = slot)
            features <- features[features %in% colnames(z)]
        } else {
            z <- obj
        }    
        z0 <- z
        if (!is.null(zero_offset)) {
            z0[z0 < zero_offset / 2] <- NA
        }
        if (!identical(cells, colnames(obj))) {
            # make sure ggplot knows the true range of the feature when subset is 
            # plotted
            limits <- quantile(z0[, features], probs = c(1 - 0.999, 0.999),
                na.rm = TRUE)
            z <- z[cells, , drop = FALSE]
            z0 <- z0[cells, , drop = FALSE]
        }
        z <- t(z)
        coords <- .parse_coords(obj, colnames(z))
        d.f <- do.call(rbind, lapply(seq(nrow(z)), function(i) 
            data.frame( 
                x = coords[,1], 
                y = 1-coords[,2], 
                fraction = z[i,],
                cluster = rownames(z)[i]
        )))
        # ggplot sets values outside limits as NA, we want to truncate them here
        if (!is.null(limits)) {
            if (is.null(zero_offset)) {
                d.f$fraction[d.f$fraction < limits[1]] <- limits[1]
            } else {
                d.f$fraction[d.f$fraction < limits[1] & d.f$fraction > zero_offset / 2] <- limits[1]
            }
            d.f$fraction[d.f$fraction > limits[2]] <- limits[2]
        }

        plot_spots(d.f, undetected_NA = undetected_NA, limits = limits, ...) 

        if (plot_map && is(d.f$fraction, "factor")) {
            gp <- ggplot(d.f, aes_string("x", "y", color = "fraction")) + 
                geom_point() +
                scale_color_discrete(name = "") +
                theme_void() + theme(legend.position = "none")
            if (nrow(z) > 1) {
                gp <- gp + facet_wrap(~cluster+fraction)
            } else {
                gp <- gp + facet_wrap(~fraction)
            }
            print(gp)    
        }
        if (length(features) > 1) {
            if (plot_violin && length(levels(obj)) > 1) {
                if (length(features) > 36) {
                    flog.warn("Too many features for violin plot.")
                } else {    
                    print(plot_violin(obj, features, cells, zero_offset))
                }
            }    
            if (plot_correlations && requireNamespace("GGally", quietly = TRUE)) {
                if (length(features) > 50) {
                    flog.warn("Too many features for correlation plot.")
                } else {    
                    print(GGally::ggcorr(data = NULL, cor_matrix = .cor_nn(obj, features, average_nn = FALSE, zero_offset = zero_offset), 
                        label = TRUE, label_size = 3, layout.exp = 3,
                        hjust = 1))
                }
                #print(GGally::ggcorr(data = NULL, cor_matrix = .cor_nn(obj, features, average_nn = TRUE, zero_offset = zero_offset), 
                #    label = TRUE, label_size = 3, layout.exp = 3,
                #    hjust = 1))
            }       
        }
    }    
    obj
}


#' plot_violin
#'
#' Violin plot of a Seurat object
#' @param obj Seurat object
#' @param features \code{character(n)} of features to be plotted
#' @param cells Plot only specified cells
#' @param zero_offset If zero values were set to a low value, specify the offset here.
#' @param pt_size Size of overlayed dots 
#' @param slot Slot to pull feature data for
#' @param ... Arguments passed to \code{Seurat::VlnPlot}.
#' @export plot_violin
#' @examples
#' plot_violin()
plot_violin <- function(obj, features, cells = NULL, zero_offset = NULL, pt_size = 0.25, slot = "data", ...) {
    if (is.null(cells)) cells <- colnames(obj)
    obj <- .remove_offset(obj, features, zero_offset, slot = slot)
    VlnPlot(obj[, cells], features = features, pt.size = pt_size, slot = slot, ...)
}

.remove_offset <- function(obj, features, zero_offset, set_na = FALSE, slot = "data") {
    if (!is.null(zero_offset)) {
        for (i in features) {
            x <- FetchData(obj, vars = i, slot = slot)
            if (set_na) {
                x[x < zero_offset / 2] <- NA
            } else {
                x[x < zero_offset / 2] <- x[x < zero_offset / 2] - zero_offset
            }    
            obj[[i]] <- x
        }
    }
    obj
}
     
#' plot_spots
#'
#' Plots spots in a Seurat object on HE slides
#' @param x \code{data.frame} with spot coordinates, value and cluster name.
#' @param labels transformation of the labels, e.g. \code{scales::percent}
#' @param labels_title Title, shown in the legend
#' @param undetected_NA Remove spots with low value by setting them to \code{NA}
#' @param undetected_cutoff Define low value by this cutoff
#' @param na.rm Hide spots with NA values
#' @param max_quantile Cap values at this quantile to avoid that outliers.
#' Ignored if \code{limits} is not \code{NULL}.     
#' define the color scheme
#' @param reorder_clusters Reorder clusters alphabetically
#' @param palette Color palette. Can be two colors in format "low:high",
#' brewer_single_hue_x (x = red, green, blue, orange, gray, purple) for
#' single hue palettes from colorbrewer2.org     
#' @param palette_inverse Flip the low and high colors in \code{palette}     
#' @param alpha Transparency of spots
#' @param size Size of spot dots     
#' @param trans Transform values, currently only log2 is implemented when this 
#' parameter is set to a different value than \code{NULL}.
#' @param limits Optional limits for continuous scores
#' @export plot_spots
#' @importFrom viridis scale_color_viridis
#' @importFrom scales percent log2_trans math_format trans_breaks
#' @import ggplot2
#' @examples
#' plot_spots()

plot_spots <- function(x, labels = scales::percent, 
                      labels_title = "Contribution", 
                      undetected_NA = TRUE, 
                      undetected_cutoff = 0.001, na.rm = TRUE,
                      max_quantile = 0.975,
                      reorder_clusters = TRUE, palette = "viridis", 
                      palette_inverse = FALSE,
                      alpha = 0.85, pt.size.factor = 0.9, trans = NULL, limits = NULL) {
    if (reorder_clusters) {
        x$cluster <- factor(as.character(x$cluster), 
            levels = .order_clusters(x$cluster))
    }
    palette <- .get_palette(palette)
    cutoffs <- NULL
    if (!is(x$fraction, "factor")) {
        if (undetected_NA) x$fraction[which(x$fraction < undetected_cutoff)] <- NA
        if (is.null(limits) && max_quantile < 1) {
            cutoffs <- quantile(x$fraction, probs = 
                c(1 - max_quantile, max_quantile), na.rm = TRUE)
            if (abs(cutoffs[1] - 1) < 0.0001 && 
                abs(min(x$fraction, na.rm = TRUE) - 0) < 0.0001) {
                cutoffs[1] <- -1
            }
            x$fraction[which(x$fraction < cutoffs[1])] <- cutoffs[1]
            x$fraction[which(x$fraction > cutoffs[2])] <- cutoffs[2]
        }        
    }
    if (na.rm) {
        x <- x[!is.na(x$fraction),]
    }
        
    fun_scale_color <- .get_scale_color_cont(x, palette, palette_inverse)

    gp <- ggplot(x, aes_string("x", "y", color = "fraction"))+
        geom_point() + theme_void()

    if (!is.null(trans)) {
        gp <- gp + fun_scale_color(labels = labels, name = labels_title,
             trans = scales::log2_trans(),
             breaks = scales::trans_breaks("log2", function(x) 2^x))
    }  else {
        gp <- gp + fun_scale_color(labels = labels, name = labels_title,
            limits = limits)
    }
    if (!is(x$fraction, "factor")) {
        gp <- gp + labs(caption = paste0("Min: ",
            round(min(x$fraction, na.rm = TRUE), digits = 2), "; Max: ", 
            round(max(x$fraction, na.rm = TRUE), digits = 2)))
    }
    if (length(levels(x$cluster)) <= 16) {
        print(gp + facet_wrap(~cluster))
    } else {
        flog.info("More than 16 clusters, paginate plot.")
        n_pages <- ceiling(length(levels(x$cluster))/16)
        for (i in seq_len(n_pages)) {
            print(gp + ggforce::facet_wrap_paginate(~cluster, ncol = 4, nrow = 4, page = i))
        }
    }
}

.order_clusters <- function(s) {
    s <- unique(as.character(s))
    if (length(s) < 2) return(s)
    s[order(gsub("_\\d+$","", s), as.numeric(gsub("^.*_","", s)))]
}

.get_image_slice <- function(obj) {
    Images(obj)[which.max(sapply(Images(obj), function(i) 
        length(intersect(Cells(obj[[i]]), Cells(obj)))))]
}

.parse_coords <- function(obj, n, cols=1:2) {
    # Visium
    obj <- obj[, n]
    if (length(Images(obj))) {
        image_idx <- .get_image_slice(obj)
        return(obj[[image_idx]]@coordinates[,c("imagecol", "imagerow")])
    }
    # ST
    n <- gsub("_\\d+$","", n)
    n <- gsub("^\\d+_","", n)
    apply(do.call(rbind, strsplit(n, split = "x")), 2, as.numeric)[,cols]
}  

#' plot_signatures
#'
#' Plots signatures defined in a GMT file and applied to a Seurat object on HE slides
#' @param obj_spatial Seurat object with SpatialTranscriptomics data
#' @param file Output PDF
#' @param gmt GMT file with gene signatures or gene signature read by
#' \code{\link{read_signatures}}
#' @param nbin Argument of \code{Seurat::AddModuleScore}
#' @param ctrl Argument of \code{Seurat::AddModuleScore}
#' @param method Either use \code{Seurat::AddModuleScore} or simple mean
#' @param width Plot width 
#' @param png Create, in addition to PDF, PNG files
#' @param cells Plot only specified cells
#' @param zero_cutoff Cutoff defining zero. Defaults to half the number of genes in the signature.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param ... Arguments passed to \code{\link{plot_spots}}.
#' @export plot_signatures
#' @examples
#' plot_signatures()

plot_signatures <- function(obj_spatial, file, gmt, nbin = 24, 
    ctrl = 30, method = c("seurat", "mean"), width = 10, png = FALSE,
    cells = NULL, zero_cutoff = NULL, assay = "Spatial", ...) {
    if (is(gmt, "character") && file.exists(gmt)) {
        sigs <- read_signatures(gmt, obj_spatial)
    } else { 
        sigs <- gmt
    }    
    sig_names <- .get_signature_names(obj_spatial, sigs)
    if (sum(is.na(sig_names))) {
        flog.info("Adding signature scores...")
        obj_spatial <- .add_module_score(obj_spatial, sigs, zero_offset = -1000,
            zero_cutoff = zero_cutoff, nbin = nbin, ctrl = ctrl, 
            name = names(sigs), method = method, assay = assay)
    } else {
        flog.info("obj_spatial contains signature scores.")
    }    
    sig_names <- .get_signature_names(obj_spatial, sigs)
    ratio <- .get_image_ratio(length(sig_names))
    obj_spatial <- .plot_spatial_with_image(file, obj_spatial, sig_names, width, ratio, cells = cells,
                  zero_offset = 0, png = png, ...)
    obj_spatial
}

.get_image_ratio <- function(l) {
    ratio <- 1
    if (( l > 9 &&  l < 13 )  || 
        ( l > 4 &&  l < 7 ) ) ratio <- 3/4
    if ( l == 3 ) ratio <- 1/3
    if ( l == 2 ) ratio <- 1/2
    ratio
}

.extract_nmf_obj <- function(obj, rank, randomize = FALSE) {
    r <- if (randomize) "_random" else ""
    if (length(rank) > 1) {
        nmf_obj <- obj@misc[[paste0("nmf_k_", min(rank), "_to_", max(rank), r)]]
    } else {
        nmf_obj <- obj@misc[[paste0("nmf_k_", rank, r)]]
    }       
    nmf_obj
}

#' plot_nmf
#'
#' Basic plots of a saved \code{NMF} object
#' @param obj Seurat object that contains a \code{NMF} object after
#' running \code{\link{cluster_nmf}}
#' @param libs Library ids, must be stored in \code{obj$library}
#' @param labels Optional \code{character(n)} with sample labels
#' @param rank Number of clusters 
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param width Output PDF width
#' @param png Create, in addition to PDF, PNG files
#' @param plot_he Plot for each library the HE jpeg
#' @param plot_ranks Plot diagnostics about all tested ranks
#' @param plot_qc Plot QC statistics about all tested ranks
#' @param assay Name of the assay corresponding to the initial input data.
#' @param ... Additional parameters passed to \code{\link{plot_features}}
#' @importFrom graphics plot
#' @importFrom grDevices png
#' @importFrom stats heatmap lm
#' @importFrom data.table melt
#' @export plot_nmf
#' @examples
#' plot_nmf()
plot_nmf <- function(obj, libs, labels = NULL, rank, prefix, 
                     subdir = "nmf", width = 10, png = FALSE,
                     plot_he = TRUE, plot_ranks = TRUE, plot_qc = TRUE, assay = "Spatial", ...) {
    nmf_obj <- .extract_nmf_obj(obj, rank)
    nmf_obj_random <- .extract_nmf_obj(obj, rank, randomize = TRUE)

    libs <- as.vector(libs)
    #create subdir if not exists
    .get_sub_path(prefix, subdir, "")

    for (k in rank) {
        features <- paste0("nmf_k_", k, "_", seq(1, k))
        obj <- .extract_nmf_r2(obj, rank, k)
        obj <- .extract_nmf_rss(obj, rank, k)
        nmf_obj_f <- if (is(nmf_obj, "NMFfit")) nmf_obj else nmf_obj$fit[[as.character(k)]]
        Idents(obj) <- predict(nmf_obj_f)

        #obj <- .rescale_features(obj, features)
        ratio <- .get_image_ratio(k)
        if (plot_he) {
            flog.info("Generating output plots for k = %i...", k)
            tmp <- .get_sub_path(prefix, file.path(subdir, "he"), "") # make sure that HE directory exists

            for (i in seq_along(libs)) {
                label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
                libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
                obj_split <- obj[,obj$library == libs[i]]
                obj_split@images <- obj_split@images[which(names(obj_split@images) %in% make.names(libs[i]))]

                filename <- .get_sub_path(prefix, file.path(subdir, "he", k), paste0("_he_nmf_cluster_", k, label, libs_label, ".pdf"))
                .plot_spatial_with_image(filename, obj_split, features, width, ratio, 
                              plot_correlations = TRUE, plot_violin = TRUE, png = png, ...)

                sd.plot <- SpatialDimPlot(obj_split, label = TRUE, image = 
                                          sttkit:::.get_image_slice(obj_split), label.size = 3, ...)
                filename <- .get_sub_path(prefix, file.path(subdir, "he", k),
                    paste0("_he_nmf_discrete_cluster_", k, label, libs_label, ".pdf"))
                pdf(filename, width = 4, height = 3.9)
                print(sd.plot)
                dev.off()
                if(png) {
                    png(sub(".pdf$", ".png", filename), width = 4, height = 3.9, units = "in", res = 150)
                    print(sd.plot)
                    dev.off()
                }
            }
        }
        if (length(features) > 2 & length(libs) > 1) {
            tmp <- .get_sub_path(prefix, file.path(subdir, "advanced"), "") # make sure that advanced directory exists
            filename <- .get_sub_path(prefix, file.path(subdir,"advanced", k), 
                paste0("_nmf_cluster_", k, "_correlations.pdf"))
            pdf(filename, width = width, height = width * ratio, onefile = FALSE)
            plot_correlation_heatmap(lapply(libs, function(i) obj[, obj$library == i]), features)
            dev.off()
            
            .plot_correlation_labels(obj, cluster_labels = predict(nmf_obj_f), 
                prefix = prefix, file.path(subdir, "advanced", k), 
                suffix = paste0("_nmf_cluster_", k, "_label_correlations.pdf")) 
        }
        write_nmf_features(obj, rank = rank, k = k, prefix = prefix)
        if (plot_qc) {
            x <- obj@meta.data
            xm <- melt(x[,grep("library|nmf|nFeature", colnames(x))], id.vars=c("library", paste0("nFeature_", assay)))
            xm <- xm[grep(paste0("nmf_k_", k, "_"), xm$variable),]
            tmp <- .get_sub_path(prefix, file.path(subdir, "qc"), "") # make sure that qc directory exists
            filename <- .get_sub_path(prefix, file.path(subdir, "qc", k), paste0("_nmf_cluster_", k, "_qc.pdf"))
            pdf(filename, width = width, height = width * ratio)
            if (length(unique(xm$library)) < 6) {
                gp <- ggplot(xm, aes_string(paste0("nFeature_", assay), "value", color = "library")) + 
                    geom_point(size=0.2) + 
                    facet_wrap(~variable) + 
                    ylab("Contribution") + 
                    theme(legend.position = "bottom")
            } else {
                gp <- ggplot(xm, aes_string(paste0("nFeature_", assay), "value")) + 
                    geom_point(size=0.2) + 
                    facet_wrap(~variable) + 
                    ylab("Contribution")
            }       
            if (requireNamespace("ggthemes", quietly = TRUE)) {
                gp <- gp + ggthemes::scale_colour_colorblind()
            }    
            print(gp)
            dev.off()
            if (png) {
                filename <- .get_sub_path(prefix, file.path(subdir, "qc", k), paste0("_nmf_cluster_", k, "_qc.png"))
                png(filename, width = width, height = width * ratio, units = "in", res = 150)
                print(gp)
                dev.off()
            }    
        }
        filename <- .get_sub_path(prefix, file.path(subdir, "advanced", k), paste0("_nmf_cluster_", k, "_coefmap.pdf"))
        pdf(filename, width = width, height = width)
        NMF::coefmap(nmf_obj_f)
        dev.off()
        if (png) {
            filename <- .get_sub_path(prefix, file.path(subdir, "advanced", k), paste0("_nmf_cluster_", k, "_coefmap.png"))
            png(filename, width = width, height = width, units = "in", res = 150)
            NMF::coefmap(nmf_obj_f)
            dev.off()
        }    
    }
    if (plot_qc) {
        .plot_nmf_r2(obj, libs, rank, prefix, file.path(subdir, "qc"), width, png, ...)
    }
    if (plot_ranks && length(rank) > 1) {
        #.plot_nmf_r2(obj, libs, rank, prefix, subdir, width, 
        #    png, feature_suffix = "rss", ...)
        filename <- .get_sub_path(prefix, file.path(subdir, "qc"), "_nmf_ranks.pdf")
        pdf(filename)
        if (is.null(nmf_obj_random)) {
            print(plot(nmf_obj))
        } else {
            print(plot(nmf_obj))
            print(plot(nmf_obj, nmf_obj_random))
        }    
        dev.off()
    }
}        

.rescale_features <- function(obj, features) {
    m <- apply(obj@meta.data[,features], 2, function(x) { 
            x[x < 1 / length(features)] = 0 
            x})
    m <- m / rowSums(m)
    obj@meta.data[,features] <- m
    obj
}

.plot_nmf_r2 <- function(obj, libs, rank, prefix, subdir, width, 
                         png, feature_suffix = "r2", ...) {

    features <- paste0("nmf_k_", feature_suffix, "_", rank)
    ratio <- .get_image_ratio(length(rank))

    labels_title <- toupper(feature_suffix)
    for (i in seq_along(libs)) {
        flog.info("Generating output %s plots for %s...", labels_title, libs[i])
        filename <- .get_sub_path(prefix, subdir,
            paste0("_he_nmf_cluster_", feature_suffix, "_", libs[i], ".pdf"))
        obj_split <- obj[,obj$library == libs[i]]
        obj_split@images <- obj_split@images[which(names(obj_split@images) %in% make.names(libs[i]))]
        .plot_spatial_with_image (filename, obj_split, features, width, ratio, plot_violin = TRUE, png = png, ...)
    #        labels = waiver(), labels_title = sprintf("%12s", labels_title), ...)
    }
}    

.plot_correlation_labels <- function(obj, cluster_labels, prefix, subdir, suffix) {
    if (!requireNamespace("corrplot")) {
        flog.warn("Package corrplot not installed.")
    } else if ("library" %in% colnames(obj@meta.data)) {
        sample_labels <- obj$library
        if ("label" %in% colnames(obj@meta.data)) {
            sample_labels <- obj$label
        }    
        chisq <- chisq.test( table(cluster_labels, sample_labels))
        contrib <- 100 * chisq$residuals^2 / chisq$statistic
        filename <- .get_sub_path(prefix, subdir, suffix)
        pdf(filename, height = 8, width = 7)
        par(mfrow = c(1, 2))
        corrplot::corrplot(contrib, is.cor = FALSE)
        corrplot::corrplot(chisq$residuals, is.cor = FALSE)
        dev.off()
    }
}

#
#' plot_correlation_heatmap
#'
#' Given a list of Seurat objects, plots the fraction of available genes
#' for each signature in a GMT file
#' @param objs List of Seurat objects
#' @param features Features to plot
#' @param sample_labels \code{character(n)} with sample labels 
#' @param feature_labels \code{character(n)} with feature labels. 
#' If \code{NULL}, just use indices.
#' @param average_nn Include the 4 nearest neighbors (top, bottom,
#' right, left) in the correlation
#' @export plot_correlation_heatmap
#' @examples
#' plot_correlation_heatmap
plot_correlation_heatmap <- function(objs, features, sample_labels = NULL, feature_labels = NULL, average_nn = FALSE) {
    if (length(features) < 3) {
        stop("Need at least 3 features.")
    }    
    mm <- do.call(cbind, lapply(objs, function(obj) { 
        x <- .cor_nn(obj, features, average_nn = average_nn)
        x[lower.tri(x)]
    }))
    if (!is.null(sample_labels)) {
        colnames(mm) <- sample_labels
    } else if (!is.null(objs[[1]]@meta.data$label)) {
        colnames(mm) <- sapply(objs, function(obj) obj$label[1])
    } else if (!is.null(objs[[1]]@meta.data$library)) {
        colnames(mm) <- sapply(objs, function(obj) obj$library[1])
    }
    if (is.null(feature_labels)) {    
        kk <- sapply(seq_along(features), function(i)
            sapply(seq_along(features), function(j) paste0(i,"_",j)))
    } else {
        kk <- sapply(seq_along(features), function(i) 
            sapply(seq_along(features), function(j) 
                paste0(feature_labels[i],"_",feature_labels[j])))
    }        
    rownames(mm) <-kk[lower.tri(kk)]
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(mm)
    } else {
        heatmap(mm)
    }   
}    

#' plot_gmt_availability
#'
#' Given a list of Seurat objects, plots the fraction of available genes
#' for each signature in a GMT file
#' @param objs Seurat object 
#' @param gmt GMT File with gene signatures, read by \code{\link{read_signatures}}
#' @param assay Name of the assay corresponding to the initial input data.
#' @export plot_gmt_availability
#' @examples
#' plot_gmt_availability
plot_gmt_availability <- function(objs, gmt, assay = "Spatial") {
    libs <- sapply(objs, function(x) x$library[1])
    names(objs) <- as.character(libs)
    gmt <- lapply(gmt, function(sig) sig <- sort(unique(sig)))
    x <- lapply(objs, function(obj) {
            check_non_zero <- min(apply(GetAssayData(obj, slot = "counts", assay = assay), 1, max)) < 1
            lapply(gmt, function(sig) {
                # gene available in obj
                idx <- sig %in% rownames(obj)
                if (check_non_zero) {
                    idx_2 <- apply(GetAssayData(obj[sig[idx], ], slot = "counts", assay = assay), 1, max) > 0
                    idx[idx] <- idx_2
                }
                idx
            })
        })       
    
    # get unavailable symbols in any library
    unavailable <- lapply(x, function(i) 
        lapply(seq_along(gmt), function(j) gmt[[j]][!i[[j]]]))
    # unavailable in all libraries 
    unavailable_all <- lapply(seq_along(gmt), function(i) 
        Reduce(intersect, lapply(unavailable, function(x) x[[i]])))
    for (i in seq_along(unavailable_all)) {
        if (length(unavailable_all[[i]])) {
            flog.warn("Features unavailable in all libraries of signature %s: %s.",
                names(gmt)[i], paste(unavailable_all[[i]], collapse = ","))
        }    
    }

    x_df <- do.call(rbind, 
        lapply(names(x), function(i) 
            do.call(rbind, lapply(names(x[[i]]), 
                function(j) data.frame(i, j, sum(x[[i]][[j]])/length(x[[i]][[j]])))
    )))
    colnames(x_df) <- c("Library", "Signature", "Available")
    x_df$Signature <- as.character(x_df$Signature)
    if (length(objs) > 1) {
        gp <- ggplot(x_df, aes_string("Library", "Available")) + 
            geom_bar(stat = "identity") + 
            facet_wrap(~Signature) + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    } else{
        gp <- ggplot(x_df, aes_string("Signature", "Available")) + 
            geom_bar(stat = "identity") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }    
    gp    
}    


#' plot_nmf_gse
#'
#' Plots the association of a gene signature with a NMF clustering
#' @param obj Results of \code{\link{calculate_nmf_gse}}
#' @param sig The name of the signature to plot. Must match a row name
#' in \code{obj}.
#' @param binwidth P-value histogram binwidth
#' @param rank Number of clusters 
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param method Method that was used. Only used in output filenames.
#' @param width Output PDF width
#' @param png Create, in addition to PDF, PNG files
#' @export plot_nmf_gse
#' @examples
#' plot_nmf_gse
plot_nmf_gse <- function(obj, sig, binwidth = 0.025, rank, prefix,
                         subdir = "nmf", method = "", width = 10, png = FALSE) {
    x_df <- data.frame(
        do.call(rbind, 
            lapply(seq_along(obj), function(i) 
                do.call(rbind, 
                    lapply(seq_along(obj[[i]]), function(j) 
                    c(K = rank[i], 
                      Cluster = j, 
                      PV = as.numeric(obj[[i]][[j]][sig, 1])
    ))))))
    ratio <- .get_image_ratio(length(rank))
    flog.info("Generating output plots for signature = %s", sig)
    if (method != "") method <- paste0(method, "_")
    filename <- .get_sub_path(prefix, subdir, paste0("_nmf_", method,
        make.names(sig), ".pdf"))

    gp <- ggplot(x_df, aes_string("PV")) + 
        geom_histogram(binwidth = binwidth) + 
        facet_wrap(~K) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        geom_vline(xintercept = 0.05, linetype="dashed", color = "red", size = 0.5) +
        ggtitle(sig)
    pdf(filename, width = width, height = width * ratio)
    print(gp)
    dev.off()
    if (png) {
        filename <- .get_sub_path(prefix, subdir, paste0("_nmf_", method, 
            make.names(sig), ".png"))
        png(filename, width = width, height = width * ratio, units = "in", res = 150)
        print(gp)
        dev.off()
    }
    gp
}    


#' plot_signatures_fake_bulk
#'
#' Plots signatures defined in a GMT file and applied to list of
#' Seurat objects, ignores spatial information. This is useful for
#' a 30,000ft look at the data, similar to bulk
#' @param objs Seurat objects with SpatialTranscriptomics data
#' @param gmt GMT file with gene signatures or gene signature read by
#' \code{\link{read_signatures}}
#' @param assay Assay to extract counts from
#' @param log_trans Transform counts into log-space.
#' @param labels Optional \code{character(n)} vector with labels
#' @param plot_pairs Plot pairwise correlations of all samples
#' @param plot_bar Plot median counts per M for each sample
#' @param plot_heatmap Plot heatmaps of all signatures in \code{gmt}
#' @param ... Arguments passed to \code{\link{plot_spots}}.
#' @export plot_signatures_fake_bulk
#' @examples
#' plot_signatures_fake_bulk()
plot_signatures_fake_bulk <- function(objs, gmt, assay = "Spatial", log_trans = TRUE, labels = NULL, 
    plot_pairs = TRUE, plot_bar = TRUE, plot_heatmaps = TRUE) {
    idx <- sapply(objs, is, "Seurat")
    objs <- objs[idx]
    if (!is.null(labels)) labels <- labels[idx]
    if(!requireNamespace("edgeR", quietly = TRUE)) {
        stop("Install the edgeR package.")
    }
    gg_data <- do.call(rbind, lapply(seq_along(objs), function(i) {
        obj <- objs[[i]]
        counts <- Matrix::rowSums(GetAssayData(obj, slot = "counts", assay = assay))
        # make sure we don't miss genes because of make.names
        if (log_trans) {
            # this silly line makes sure that we use log, not log2
            # as Seurat does and that we deal with the pseudo-counts
            counts <- log(2^edgeR::cpm(counts, log = TRUE))
        } else { 
            counts <- edgeR::cpm(counts)
        }    
        sample_id <- as.character(obj$library[1])
        if (!is.null(labels)) {
            sample_id <- labels[[i]]
        } else if (!is.null(obj@meta.data$label)) {
            sample_id <- obj@meta.data$label[1]
        }
        data.frame(Symbol = make.names(rownames(counts)),
                   Id = sample_id, 
                   Spatial = counts[, 1])
    }))
    universe <- make.names(Reduce(intersect, lapply(objs, rownames)))
    gg_data <- gg_data[gg_data[,1] %in% universe,]

    gg_data_casted <- data.frame(data.table::dcast(data=gg_data, 
        formula = Symbol~Id, value.var = "Spatial"), row.names = 1)

    if (!log_trans) { 
        gg_data_casted[is.na(gg_data_casted)] <- 0
    } else {
        gg_data_casted[is.na(gg_data_casted)] <- min(gg_data_casted, na.rm = TRUE)
    }

    if(plot_pairs && requireNamespace("GGally", quietly = TRUE)) {
        print(GGally::ggpairs(gg_data_casted, xlab = "log counts per M", ylab="log count per M"))
    }
    x_melt <- NULL
    if (plot_bar) {
        if (is.null(gmt)) stop("plot_bar requires gmt.")
        x <- sapply(gmt, function(x) find_signature_means(gg_data_casted, features = make.names(x), fun = median))
        x <- data.frame(Id=unique(gg_data$Id), x)
        x_melt <- melt(x)
        x_melt$Id <- factor(x_melt$Id, levels=sort(levels(x_melt$Id)))
        print(ggplot(x_melt,aes(Id, value))+
            geom_col()+
            facet_wrap(~variable, scale="free_y")+
            xlab("")+
            ylab("Median counts per M") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)))
            
    }
    if (plot_heatmaps && requireNamespace("pheatmap", quietly = TRUE) && ncol(gg_data_casted)>=2) {
        if (is.null(gmt)) stop("plot_heatmaps requires gmt.")
        gmt_f <- lapply(gmt, function(x) make.names(x)[make.names(x) %in% rownames(gg_data_casted)])
        gmt_f <- gmt_f[sapply(gmt_f, length)>1]
        for (i in seq_along(gmt_f)) {
            pheatmap::pheatmap(gg_data_casted[gmt_f[[i]],], main = names(gmt_f)[i], cluster_rows = length(gmt_f[[i]])>1)
        }
    }
    list(genes = gg_data_casted, gmt = x_melt)
}


#' plot_signatures_nmf
#'
#' Plot correlation of gene signatures with basis of an \code{NMF} object
#' @param obj Seurat object that contains a \code{NMF} object after
#' running \code{\link{cluster_nmf}}
#' @param gmt GMT file with gene signatures or gene signature read by
#' \code{\link{read_signatures}}
#' @param rank Number of clusters 
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param width Output PDF width
#' @param png Create, in addition to PDF, PNG files
#' @export plot_signatures_nmf
#' @examples
#' plot_signatures_nmf()
plot_signatures_nmf <- function(obj, gmt, gmt_name = NULL, rank, prefix, 
                     subdir = "nmf", width = 10) {
    nmf_obj <- .extract_nmf_obj(obj, rank)
    if (is.null(nmf_obj)) {
        flog.warn("Cannot find NMF run in provided obj.")
        return(NULL)
    }    
    nmf_obj_random <- .extract_nmf_obj(obj, rank, randomize = TRUE)
    if (is(gmt, "character")) {
        sigs <- read_signatures(gmt, obj_spatial)
        if (is.null(gmt_name)) gmt_name <- gsub(".gmt", "", basename(gmt))
    } else { 
        sigs <- gmt
        if (is.null(gmt_name)) gmt_name <- ""
    }    
    if (gmt_name != "") gmt_name <- paste0("_", gmt_name)
    d_f <- do.call(rbind, lapply(rank, function(k) {
        nmf_obj_f <- if (is(nmf_obj, "NMFfit")) nmf_obj else nmf_obj$fit[[as.character(k)]]
        b <- NMF::basis(nmf_obj_f)
        sigsf <- lapply(sigs, function(x) x[x %in% rownames(b)])
        sigsf <- sigsf[sapply(sigsf, length)>0]
        reshape2::melt(do.call(rbind, lapply(seq_along(sigsf), 
            function(i) data.frame(Signature = names(sigsf)[i], 
                                   Rank = k,  
                                   K = b[sigsf[[i]],,drop = FALSE]))),
            id.vars=c("Signature", "Rank"))
    }))
    ratio <- .get_image_ratio(length(sigs))
    for (k in rank) {
        filename <- .get_sub_path(prefix, file.path(subdir, k), paste0("_nmf_signature_", k, gmt_name, ".pdf"))
        pdf(filename, width = width, height = width * ratio)
        gp <- ggplot(d_f[d_f$Rank == k,], 
            aes(as.character(variable), value)) + geom_boxplot() + 
            xlab("") +
            ylab("NMF basis") +
            ggtitle(sig)
        if (length(sigs) <= 16) {
            print(gp + facet_wrap(~Signature, scales = "free_y"))
        } else {
            flog.info("More than 16 signatures, paginate plot.")
            n_pages <- ceiling(length(sigs) / 16)
            for (i in seq_len(n_pages)) {
                print(gp + ggforce::facet_wrap_paginate(~Signature, scales = "free_y", ncol = 4, nrow = 4, page = i))
            }
        }
        dev.off()
    }    
}    

.get_palette <- function(x) {
    switch(x,
        "brewer_single_hue_red" = "#fee0d2:#de2d26",
        "brewer_single_hue_blue" = "#deebf7:#3182bd",
        "brewer_single_hue_green" = "#e5f5e0:#31a354",
        "brewer_single_hue_gray" = "#f0f0f0:#636363",
        "brewer_single_hue_grey" = "#f0f0f0:#636363",
        "brewer_single_hue_orange" = "#fee6ce:#e6550d",
        "brewer_single_hue_purple" = "#efedf5:#756bb1",
        x)
}        
.get_scale_color_cont <- function(x, palette, palette_inverse) {
    fun_scale_color <- scale_color_continuous
    palette_split <- strsplit(palette, ":")[[1]]

    if (palette_inverse) palette_split <- rev(palette_split)
        
    if (is(x$fraction, "factor")) {
        fun_scale_color <- scale_color_discrete
    } else if (palette == "viridis") {
        fun_scale_color <- viridis::scale_color_viridis
    } else if (length(palette_split)==2) {
        fun_scale_color <- function(...) scale_color_gradient(low  = palette_split[1],
                                                high = palette_split[2], ...)
    }
    return(fun_scale_color)
}        

.plot_spatial_with_image <- function(filename, object, features, width = 10, ratio = 0.5,
                                 ncol = 4, nrow = 4, cells = NULL, zero_offset = NULL, 
                                 png = FALSE, plot_correlations = FALSE, plot_violin = FALSE, ...) {
    object_resized <- .resize_slice_images(object)
    if (is.null(cells)) cells <- colnames(object_resized)

    if (length(features) > 16) {
        glist <- SpatialFeaturePlot(object_resized, image = .get_image_slice(object_resized), 
            features = features, combine = FALSE, ...)
        glist <- lapply(glist, ggplotGrob)
        glist <- marrangeGrob(glist, ncol = ncol, nrow = nrow)
        ggsave(filename, glist,
               width = width, height = width * ratio)
        if (png) {
            invisible(mapply(ggsave, filename = paste0(gsub(".pdf$", "", filename), "_", seq_along(glist),".png"),
                   plot = glist, 
                   MoreArgs = list(width = width, height = width * ratio, units = "in", dpi = 150)))
        }    
        flog.warn("Too many features for violin plot. Skipping...")
        flog.warn("Too many features for correlation plot. Skipping...")
    } else {
        pdf(filename, width = width, height = width * ratio)
        print(SpatialFeaturePlot(object_resized, image = .get_image_slice(object_resized), 
              features = features, combine = TRUE, ...))
        if (length(features) > 1) {
            if (plot_violin && length(levels(object_resized)) > 1)
                print(plot_violin(object_resized, features, cells, zero_offset))
            if (plot_correlations && requireNamespace("GGally", quietly = TRUE))
                print(GGally::ggcorr(data = NULL, cor_matrix = .cor_nn(object_resized, 
                      features, average_nn = FALSE, zero_offset = zero_offset), 
                      label = TRUE, label_size = 3, layout.exp = 3, hjust = 1))
        }
        dev.off()

        if(png) {
            png(gsub(".pdf$", ".png", filename), width = width, height = width * ratio, units = "in", res = 150)
            print(SpatialFeaturePlot(object_resized, image = .get_image_slice(object_resized),
                features = features, combine = TRUE, ...))
            dev.off()
        }
    }
    object_resized
}

#' plot_spatially_variable
#'
#' Plots results of spatial differential expression analysis
#' @param obj Seurat object 
#' @param labels Optional \code{character(n)} vector with labels
#' @param spatial_features Results of\code{Seurat::SpatiallyVariableFeatures}.
#' @param method Method used to find spatially variable features
#' @param number_features Plot that many features.
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param width Output PDF width
#' @param ncol Number of columns to plot genes
#' @param nrow Number of rows to plot genes
#' @param ... Arguments passed to \code{\link{plot_features}}
#' @importFrom gridExtra marrangeGrob
#' @export plot_spatially_variable
#' @examples
#' plot_spatially_variable
plot_spatially_variable <- function(ndata, labels = NULL, spatial_features, method = "markvariogram", 
    number_features = 80, prefix, subdir = "spatial_variation", width = 10, ncol = 4, nrow = 4, ...) {
    top_features <- lapply(spatial_features, function(x) head(x, number_features))
    top_features <- lapply(top_features, .reorder_spatially_variable_features, ndata)
    ndata_split <- SplitObject(ndata, split.by = "library")
    libs <- as.character(sapply(ndata_split, function(x) x$library[1]))
    for (i in seq_along(ndata_split)) {
        flog.info("Generating output plots for %s ...", libs[i])
        ratio <- .get_image_ratio(length(top_features[i]))
        label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
        libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
        filename <- sttkit:::.get_sub_path(prefix, subdir, 
            paste0("_he_variable_", method, label, libs_label, ".pdf"))
            .plot_spatial_with_image(filename, ndata_split[[i]], top_features[[i]], 
                                     width, ratio, ncol, nrow, ...)
    }
}

.reorder_spatially_variable_features <- function(features, x) {
    m1 <- FetchData(x, features)
    idx <- colnames(x@meta.data)[grep("nmf", colnames(x@meta.data))]
    if (length(idx)) {
        m2 <- FetchData(x, idx)
        d <- dist(t(apply(m1, 2, function(x) cor(x, m2))))
    } else {
        d <- dist(t(m1))
    }
    hc <- hclust(d, method = "ward.D2")
    hc$height <- round(hc$height, 6) # sometimes encounter rounding error
    ct <- cutree(hc, h = median(hc$height))
    ct <- split(features, ct[features])
    ranking <- seq_along(features)
    names(ranking) <- features
    idx <- order(sapply(ct, function(x) mean(ranking[x])))
    as.character(unlist(ct[idx]))
}


.resize_slice_images <- function(obj, w = 300) {
    if (!requireNamespace("EBImage", quietly = TRUE)) return(obj)
    .resize_image <- function(k) {
        new_k <- paste0(k, "_scaled")
        obj@images[[new_k]] <- obj@images[[k]]
        obj@images[[new_k]]@image <- EBImage::resize( obj@images[[k]]@image, w = w)
        r <- w / nrow(obj@images[[k]]@image)
        obj@images[[k]] <- NULL
        obj@images[[new_k]]@scale.factors$lowres <- obj@images[[new_k]]@scale.factors$lowres * r
        obj
    }
    all_images <- Images(obj)
    for(i in all_images) {
        obj <- .resize_image(i)
    }    
    obj
}

