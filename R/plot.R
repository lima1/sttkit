#' plot_features
#'
#' Plots features in a Seurat object on HE slides. Supports multi-page
#' plot when lots of features are requested.
#' @param object Seurat object
#' @param features \code{character(n)} of features to be plotted
#' @param cells Plot only specified cells
#' @param labels transformation of the labels, e.g. \code{scales::percent}
#' @param labels_title Title, shown in the legend
#' @param palette Color palette. Can be two colors in format "low:high",
#' brewer_single_hue_x (x = red, green, blue, orange, gray, purple) for
#' single hue palettes from colorbrewer2.org
#' @param palette_inverse Flip the low and high colors in \code{palette}
#' @param trans Transform values, currently only log2 is implemented when this
#' parameter is set to a different value than \code{NULL}.
#' @param ggcode Extra code added to \code{ggplot} call.
#' @param prefix Output file prefix
#' @param suffix Output file suffix
#' @param subdir Put files in a subdirectory
#' @param width Output PDF width
#' @param height Output PDF height
#' @param ncol Number of columns to plot genes for multi-page plots
#' @param nrow Number of rows to plot genes for multi-page plots
#' @param pdf Create PDF image files
#' @param png Create PNG image files
#' @param plot_correlations Plot pairwise correlations of feature scores
#' @param plot_violin Plot distribution of feature scores across idents
#' @param slot Slot to pull feature data for
#' @param ... Arguments passed to \code{SpatialFeaturePlot}.
#' @export plot_features
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom patchwork wrap_plots
#' @examples
#' plot_features()
plot_features <- function(object, features, cells = NULL,
                          labels = NULL, labels_title = "",
                          palette = NULL, palette_inverse = FALSE,
                          trans = NULL, ggcode = NULL,
                          prefix, suffix,
                          subdir = "", width = 10, height = NULL, pdf = FALSE,
                          png = FALSE,
                          ncol = 4, nrow = 4,
                          plot_correlations = FALSE, plot_violin = FALSE,
                          slot = "data", ...)  {
    filename <- .get_sub_path(prefix, subdir, suffix)
    ratio <- if (!is.null(height)) height / width else .get_image_ratio(length(features))
    .plot_spatial_with_image(filename, object = object[, cells], features = features,
                  width = width, ratio = ratio,
                  ncol = ncol, nrow = nrow, cells = cells,
                  labels = labels, labels_title = labels_title,
                  palette = palette, palette_inverse = palette_inverse, trans = trans,
                  ggcode = ggcode, plot_correlations = TRUE, plot_violin = TRUE,
                  pdf = pdf, png = png, ...)
}


#' plot_violin
#'
#' Violin plot of a Seurat object
#' @param object Seurat object
#' @param features \code{character(n)} of features to be plotted
#' @param cells Plot only specified cells
#' @param pt_size Size of overlayed dots
#' @param slot Slot to pull feature data for
#' @param ... Arguments passed to \code{Seurat::VlnPlot}.
#' @export plot_violin
#' @examples
#' plot_violin()
plot_violin <- function(object, features, cells = NULL, pt_size = 0.25, slot = "data", ...) {
    if (is.null(cells)) cells <- colnames(object)
    VlnPlot(object[, cells], features = features, pt.size = pt_size, slot = slot, ...)
}

.order_clusters <- function(s) {
    s <- unique(as.character(s))
    if (length(s) < 2) return(s)
    s[order(gsub("_\\d+$", "", s), as.numeric(gsub("^.*_", "", s)))]
}

.get_image_slice <- function(object) {
    Images(object)[which.max(sapply(Images(object), function(i)
        length(intersect(Cells(object[[i]]), Cells(object)))))]
}

.parse_coords <- function(object, n, cols=1:2) {
    # Visium
    object <- object[, n]
    if (length(Images(object))) {
        image_idx <- .get_image_slice(object)
        return(object[[image_idx]]@coordinates[, c("imagecol", "imagerow")])
    }
    # ST
    n <- gsub("_\\d+$", "", n)
    n <- gsub("^\\d+_", "", n)
    apply(do.call(rbind, strsplit(n, split = "x")), 2, as.numeric)[, cols]
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
#' @param pdf Create PDF image files
#' @param png Create PNG image files
#' @param cells Plot only specified cells
#' @param zero_cutoff Cutoff defining zero. Defaults to half the number of genes in the signature.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param ... Arguments passed to \code{Seurat::SpatialPlot}.
#' @export plot_signatures
#' @examples
#' plot_signatures()

plot_signatures <- function(obj_spatial, file, gmt, nbin = 24,
    ctrl = 30, method = c("seurat", "mean"), width = 10, pdf = FALSE, png = FALSE,
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
                  pdf = pdf, png = png, ...)
    obj_spatial
}

.get_image_ratio <- function(l) {
    ratio <- 1
    if ((l > 9 &&  l < 13)  ||
        (l > 4 &&  l < 7)) ratio <- 3/4
    if (l == 3) ratio <- 1/3
    if (l == 2) ratio <- 1/2
    ratio
}

.extract_nmf_obj <- function(object, rank, randomize = FALSE) {
    r <- if (randomize) "_random" else ""
    if (length(rank) > 1) {
        nmf_obj <- object@misc[[paste0("nmf_k_", min(rank), "_to_", max(rank), r)]]
    } else {
        nmf_obj <- object@misc[[paste0("nmf_k_", rank, r)]]
    }       
    nmf_obj
}

#' plot_nmf
#'
#' Basic plots of a saved \code{NMF} object
#' @param object Seurat object that contains a \code{NMF} object after
#' running \code{\link{cluster_nmf}}
#' @param libs Library ids, must be stored in \code{object$library}
#' @param labels Optional \code{character(n)} with sample labels
#' @param rank Number of clusters 
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param width Output PDF width
#' @param pdf Create PDF image files
#' @param png Create PNG image files
#' @param plot_he Plot for each library the HE jpeg
#' @param plot_ranks Plot diagnostics about all tested ranks
#' @param plot_qc Plot QC statistics about all tested ranks
#' @param assay Name of the assay corresponding to the initial input data.
#' @param ... Additional parameters passed to \code{\link{plot_features}}
#' @importFrom graphics plot
#' @importFrom grDevices png
#' @importFrom stats heatmap lm chisq.test cutree hclust median
#' @importFrom data.table melt
#' @export plot_nmf
#' @examples
#' plot_nmf()
plot_nmf <- function(object, libs, labels = NULL, rank, prefix, 
                     subdir = "nmf", width = 10, pdf = FALSE, png = FALSE,
                     plot_he = TRUE, plot_ranks = TRUE, plot_qc = TRUE, assay = "Spatial", ...) {
    nmf_obj <- .extract_nmf_obj(object, rank)
    nmf_obj_random <- .extract_nmf_obj(object, rank, randomize = TRUE)

    libs <- as.vector(libs)
    #create subdir if not exists
    .get_sub_path(prefix, subdir, "")

    for (k in rank) {
        features <- paste0("nmf_k_", k, "_", seq(1, k))
        object <- .extract_nmf_r2(object, rank, k)
        object <- .extract_nmf_rss(object, rank, k)
        nmf_obj_f <- if (is(nmf_obj, "NMFfit")) nmf_obj else nmf_obj$fit[[as.character(k)]]
        object <- set_idents_nmf(object, k, rank)

        tmp <- .get_sub_path(prefix, file.path(subdir, "umap"), "") # make sure that umap directory exists
        tmp <- .get_sub_path(prefix, file.path(subdir, "heatmap"), "") # make sure that umap directory exists
        filename <- .get_sub_path(prefix, file.path(subdir, "umap", k),
            paste0("_spatial_cluster.pdf"))
        if (pdf) {
            pdf(filename, width = 10, height = 5)
            gp <- DimPlot(object, reduction = .get_reduction(object), label = TRUE)
            if (requireNamespace("ggthemes", quietly = TRUE) &&
                length(levels(Idents(object))) <= 8) {
                gp <- gp + ggthemes::scale_colour_colorblind()
            }
            print(gp)
            dev.off()
        }
        if ("label" %in% colnames(object@meta.data)) {
            .plot_cluster_library(object, field = "label", prefix = prefix, 
                subdir = file.path(subdir, "umap", k),
                reference_technology = "spatial", pdf = pdf, png = png)
        } else if ("library" %in% colnames(object@meta.data)) {
            .plot_cluster_library(object, field = "library", prefix = prefix, 
                subdir = file.path(subdir, "umap", k),
                reference_technology = "spatial", pdf = pdf, png = png)
        }

        #object <- .rescale_features(object, features)
        ratio <- .get_image_ratio(k)
        features_nmf <- write_nmf_features(object, rank = rank, k = k, prefix = prefix)
        features_nmf <- unique(features_nmf$Gene[features_nmf$Gene %in% rownames(object)])
        if (plot_he) {
            flog.info("Generating output plots for k = %i...", k)
            tmp <- .get_sub_path(prefix, file.path(subdir, "he"), "") # make sure that HE directory exists

            for (i in seq_along(libs)) {
                label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
                libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
                obj_split <- object[,object$library == libs[i]]
                #obj_split@images <- obj_split@images[which(names(obj_split@images) %in% make.names(libs[i]))]

                filename <- .get_sub_path(prefix, file.path(subdir, "he", k), paste0("_he_nmf_cluster_", k, label, libs_label, ".pdf"))
                .plot_spatial_with_image(filename, obj_split, features, width, ratio, 
                              plot_correlations = TRUE, plot_violin = TRUE, pdf = pdf, png = png, ...)

                sd.plot <- SpatialDimPlot(obj_split, label = TRUE, images = 
                                          .get_image_slice(obj_split), label.size = 3, ...)

                if (requireNamespace("ggthemes", quietly = TRUE) &&
                    length(levels(Idents(obj_split))) <= 8) {
                    sd.plot <- sd.plot + ggthemes::scale_fill_colorblind()
                }    
                filename <- .get_sub_path(prefix, file.path(subdir, "he", k),
                    paste0("_he_nmf_discrete_cluster_", k, label, libs_label, ".pdf"))
                if (pdf) {
                    pdf(filename, width = 4, height = 3.9)
                    print(sd.plot)
                    invisible(dev.off())
                }
                if(png) {
                    png(sub(".pdf$", ".png", filename), width = 4, height = 3.9, units = "in", res = 150)
                    print(sd.plot)
                    invisible(dev.off())
                }
                filename <- .get_sub_path(prefix, file.path(subdir, "heatmap", k),
                    paste0("_heatmap_nmf_discrete_cluster_", k, label, libs_label, ".pdf"))
                if (ncol(obj_split) > 5000) {
                    obj_split <- subset(obj_split, downsample = 300)
                }
                gp <- DoHeatmap(obj_split, features = features_nmf)
                if (pdf) {
                    pdf(filename)
                    print(gp)
                    invisible(dev.off())
                }
                if(png) {
                    png(sub(".pdf$", ".png", filename), width = 7, height = 7, units = "in", res = 150)
                    print(gp)
                    invisible(dev.off())
                }
            }
        }
        if (length(features) > 2 & length(libs) > 1) {
            tmp <- .get_sub_path(prefix, file.path(subdir, "advanced"), "") # make sure that advanced directory exists
            filename <- .get_sub_path(prefix, file.path(subdir,"advanced", k), 
                paste0("_nmf_cluster_", k, "_correlations.pdf"))
            if (pdf) {
                pdf(filename, width = width, height = width * ratio, onefile = FALSE)
                plot_correlation_heatmap(lapply(libs, function(i) object[, object$library == i]), features)
                dev.off()
            }   
            .plot_correlation_labels(object, cluster_labels = Idents(object), 
                prefix = prefix, file.path(subdir, "advanced", k), 
                suffix = paste0("_nmf_cluster_", k, "_label_correlations.pdf"), pdf = pdf, png = png) 
        }
        if (plot_qc) {
            x <- object@meta.data
            xm <- melt(x[,grep("library|nmf|nFeature", colnames(x))], id.vars=c("library", paste0("nFeature_", assay)))
            xm <- xm[grep(paste0("nmf_k_", k, "_"), xm$variable),]
            tmp <- .get_sub_path(prefix, file.path(subdir, "qc"), "") # make sure that qc directory exists
            filename <- .get_sub_path(prefix, file.path(subdir, "qc", k), paste0("_nmf_cluster_", k, "_qc.pdf"))
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
            if (pdf) {
                pdf(filename, width = width, height = width * ratio)
                print(gp)
                dev.off()
            }
            if (png) {
                filename <- .get_sub_path(prefix, file.path(subdir, "qc", k), paste0("_nmf_cluster_", k, "_qc.png"))
                png(filename, width = width, height = width * ratio, units = "in", res = 150)
                print(gp)
                dev.off()
            }    
        }
        filename <- .get_sub_path(prefix, file.path(subdir, "advanced", k), paste0("_nmf_cluster_", k, "_coefmap.pdf"))
        if (pdf) {
            pdf(filename, width = width, height = width)
            NMF::coefmap(nmf_obj_f)
            dev.off()
        }
        if (png) {
            filename <- .get_sub_path(prefix, file.path(subdir, "advanced", k), paste0("_nmf_cluster_", k, "_coefmap.png"))
            png(filename, width = width, height = width, units = "in", res = 150)
            NMF::coefmap(nmf_obj_f)
            dev.off()
        }    
    }
    if (plot_qc) {
        .plot_nmf_r2(object, libs, rank, prefix, file.path(subdir, "qc"), width, pdf, png, ...)
    }
    if (plot_ranks && length(rank) > 1 && pdf) {
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

.rescale_features <- function(object, features) {
    m <- apply(object@meta.data[,features], 2, function(x) { 
            x[x < 1 / length(features)] = 0 
            x})
    m <- m / rowSums(m)
    object@meta.data[,features] <- m
    object
}

.plot_nmf_r2 <- function(object, libs, rank, prefix, subdir, width, 
                         pdf, png, feature_suffix = "r2", ...) {

    features <- paste0("nmf_k_", feature_suffix, "_", rank)
    ratio <- .get_image_ratio(length(rank))

    labels_title <- toupper(feature_suffix)
    for (i in seq_along(libs)) {
        flog.info("Generating output %s plots for %s...", labels_title, libs[i])
        filename <- .get_sub_path(prefix, subdir,
            paste0("_he_nmf_cluster_", feature_suffix, "_", libs[i], ".pdf"))
        obj_split <- object[,object$library == libs[i]]
        #obj_split@images <- obj_split@images[which(names(obj_split@images) %in% make.names(libs[i]))]
        .plot_spatial_with_image(filename, obj_split, features, width, ratio, plot_violin = TRUE, pdf = pdf, png = png, ...)
    #        labels = waiver(), labels_title = sprintf("%12s", labels_title), ...)
    }
}    

.plot_correlation_labels <- function(object, cluster_labels, prefix, subdir, suffix, pdf = FALSE, png = FALSE) {
    if (!requireNamespace("corrplot")) {
        flog.warn("Package corrplot not installed.")
    } else if ("library" %in% colnames(object@meta.data)) {
        sample_labels <- object$library
        if ("label" %in% colnames(object@meta.data)) {
            sample_labels <- object$label
        }    
        chisq <- chisq.test( table(cluster_labels, sample_labels))
        contrib <- 100 * chisq$residuals^2 / chisq$statistic
        filename <- .get_sub_path(prefix, subdir, suffix)
        if (pdf) {
            pdf(filename, height = 8, width = 7)
            par(mfrow = c(1, 2))
            corrplot::corrplot(contrib, is.cor = FALSE, col = "black", cl.pos="n")
            corrplot::corrplot(chisq$residuals, is.cor = FALSE, col = Seurat:::FeaturePalettes$Spatial)
            dev.off()
        }
        if (png) {
            png(paste0(gsub(".pdf$", "", filename), ".png"), width = 7, height = 8, units = "in", res = 150)
            par(mfrow = c(1, 2))
            corrplot::corrplot(contrib, is.cor = FALSE, col = "black", cl.pos="n")
            corrplot::corrplot(chisq$residuals, is.cor = FALSE, col = Seurat:::FeaturePalettes$Spatial)
            dev.off()
        }
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
#' @importFrom stats var
#' @export plot_correlation_heatmap
#' @examples
#' plot_correlation_heatmap
plot_correlation_heatmap <- function(objs, features, sample_labels = NULL, feature_labels = NULL, average_nn = FALSE) {
    if (length(features) < 3) {
        stop("Need at least 3 features.")
    }    
    mm <- do.call(cbind, lapply(objs, function(object) { 
        x <- .cor_nn(object, features, average_nn = average_nn)
        x[lower.tri(x)]
    }))
    if (!is.null(sample_labels)) {
        colnames(mm) <- sample_labels
    } else if (!is.null(objs[[1]]@meta.data$label)) {
        colnames(mm) <- sapply(objs, function(object) object$label[1])
    } else if (!is.null(objs[[1]]@meta.data$library)) {
        colnames(mm) <- sapply(objs, function(object) object$library[1])
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
    x <- lapply(objs, function(object) {
            check_non_zero <- min(apply(GetAssayData(object, slot = "counts", assay = assay), 1, max)) < 1
            lapply(gmt, function(sig) {
                # gene available in object
                idx <- sig %in% rownames(object)
                if (check_non_zero) {
                    idx_2 <- apply(GetAssayData(object[sig[idx], ], slot = "counts", assay = assay), 1, max) > 0
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
        x_df$Signature <- factor(x_df$Signature,
            levels = unique(x_df$Signature[order(x_df$Available)]))
        gp <- ggplot(x_df, aes_string("Signature", "Available")) + 
            geom_bar(stat = "identity") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }    
    gp    
}    


#' plot_nmf_gse
#'
#' Plots the association of a gene signature with a NMF clustering
#' @param object Results of \code{\link{calculate_nmf_gse}}
#' @param sig The name of the signature to plot. Must match a row name
#' in \code{object}.
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
plot_nmf_gse <- function(object, sig, binwidth = 0.025, rank, prefix,
                         subdir = "nmf", method = "", width = 10, png = FALSE) {
    x_df <- data.frame(
        do.call(rbind, 
            lapply(seq_along(object), function(i) 
                do.call(rbind, 
                    lapply(seq_along(object[[i]]), function(j) 
                    c(K = rank[i], 
                      Cluster = j, 
                      PV = as.numeric(object[[i]][[j]][sig, 1])
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
#' @param plot_heatmaps Plot heatmaps of all signatures in \code{gmt}
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
        object <- objs[[i]]
        counts <- Matrix::rowSums(GetAssayData(object, slot = "counts", assay = assay))
        # make sure we don't miss genes because of make.names
        if (log_trans) {
            # this silly line makes sure that we use log, not log2
            # as Seurat does and that we deal with the pseudo-counts
            counts <- log(2^edgeR::cpm(counts, log = TRUE))
        } else { 
            counts <- edgeR::cpm(counts)
        }    
        sample_id <- as.character(object$library[1])
        if (!is.null(labels)) {
            sample_id <- labels[[i]]
        } else if (!is.null(object@meta.data$label)) {
            sample_id <- object@meta.data$label[1]
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
        x <- lapply(gmt, function(x) find_signature_means(gg_data_casted, features = make.names(x), fun = median))
        x <- do.call(rbind, x)
        x_melt <- melt(x)
        colnames(x_melt)[1:2] <- c("variable", "Id")
        x_melt$Id <- factor(x_melt$Id, levels=sort(levels(x_melt$Id)))
        gp <- ggplot(x_melt,aes(Id, value))+
            geom_col()+
            xlab("")+
            ylab("Median counts per M") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

        if (length(gmt) <= 16) {
            print(gp + facet_wrap(~variable, scales = "free_y"))
        } else {
            flog.info("More than 16 signatures, paginate plot.")
            n_pages <- ceiling(length(gmt) / 16)
            for (i in seq_len(n_pages)) {
                print(gp + ggforce::facet_wrap_paginate(~variable, scales = "free_y", ncol = 4, nrow = 4, page = i))
            }
        }
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
#' @param object Seurat object that contains a \code{NMF} object after
#' running \code{\link{cluster_nmf}}
#' @param gmt GMT file with gene signatures or gene signature read by
#' \code{\link{read_signatures}}
#' @param gmt_name Optional name of the signature set in \code{gmt}
#' @param rank Number of clusters 
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param width Output PDF width
#' @export plot_signatures_nmf
#' @examples
#' plot_signatures_nmf()
plot_signatures_nmf <- function(object, gmt, gmt_name = NULL, rank, prefix, 
                     subdir = "nmf", width = 10) {
    nmf_obj <- .extract_nmf_obj(object, rank)
    if (is.null(nmf_obj)) {
        flog.warn("Cannot find NMF run in provided obj.")
        return(NULL)
    }    
    nmf_obj_random <- .extract_nmf_obj(object, rank, randomize = TRUE)
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
        d_f$variable <- as.character(d_f$variable)
        gp <- ggplot(d_f[d_f$Rank == k,], 
            aes(variable, value)) + geom_boxplot() + 
            xlab("") +
            ylab("NMF basis") 
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
        object <- set_idents_nmf(object, k, rank)
        filename <- .get_sub_path(prefix, file.path(subdir, k), paste0("_nmf_signature_heatmap_", k, gmt_name, ".pdf"))
        pdf(filename, width = width, height = width)
        field <- if ("label" %in% colnames(object@meta.data)) "label" else "library"
        # collapse replicates in heatmap
        field_prefix <- sort(unique(gsub("_\\d+$","", object[[field]][,1])))
        cells <- lapply(field_prefix, function(i) 
            Cells(object)[grep(i, object[[field]][,1], fixed = TRUE)])
        names(cells) <- field_prefix

        for (i in seq_along(sigs)) {
            gp <- lapply(seq_along(cells), function(j) DoHeatmap(object, features = sigs[[i]], cells = cells[[j]]) + 
                theme(legend.position = "none") + 
                ggtitle(paste(names(sigs)[[i]], names(cells)[j]))
            )
            print(wrap_plots(gp, ncol = 1))
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
    fun_scale_color <- if (is(x, "factor")) scale_color_discrete else scale_color_continuous
    if (is.null(palette)) return(fun_scale_color)

    palette_split <- strsplit(palette, ":")[[1]]

    if (palette_inverse) palette_split <- rev(palette_split)
        
    if (palette == "viridis") {
        fun_scale_color <- viridis::scale_color_viridis
    } else if (length(palette_split) == 2) {
        fun_scale_color <- function(...) scale_color_gradient(low  = palette_split[1],
                                                high = palette_split[2], ...)
    }
    return(fun_scale_color)
}        

.plot_spatial_with_image <- function(filename, object, features, width = 10, ratio = 0.5,
                                 ncol = 4, nrow = 4, cells = NULL, 
                                 labels = NULL, labels_title = "",
                                 palette = NULL, palette_inverse = FALSE, trans = NULL, ggcode = NULL,
                                 pdf = FALSE, png = FALSE, plot_correlations = FALSE, plot_violin = FALSE, ...) {
    fun_scale_color <- .get_scale_color_cont(FALSE, palette, palette_inverse) 
    .format_gp <- function(gp) {
        if (!is.null(trans)) {
            gp <- gp + fun_scale_color(labels = labels, name = labels_title,
                 trans = scales::log2_trans(),
                 breaks = scales::trans_breaks("log2", function(x) 2^x))
        }  else {
            gp <- gp + fun_scale_color(labels = labels, name = labels_title)
        }
        gp + ggcode
    }
    object_resized <- .resize_slice_images(object)
    if (is.null(cells)) cells <- colnames(object_resized)

    glist <- SpatialPlot(object_resized, images = .get_image_slice(object_resized), 
        features = features, combine = FALSE, ...)
    glist <- lapply(glist, .format_gp)
    if (length(features) > ncol * nrow) {
        glist <- lapply(glist, ggplotGrob)
        glist <- marrangeGrob(glist, ncol = ncol, nrow = nrow)
        if (pdf) {
            ggsave(filename, glist,
                   width = width, height = width * ratio)
        }
        if (png) {
            invisible(mapply(ggsave, filename = paste0(gsub(".pdf$", "", filename), "_", seq_along(glist),".png"),
                   plot = glist, 
                   MoreArgs = list(width = width, height = width * ratio, units = "in", dpi = 150)))
        }    
        flog.warn("Too many features for violin plot. Skipping...")
        flog.warn("Too many features for correlation plot. Skipping...")
    } else {
        if (pdf) {
            pdf(filename, width = width, height = width * ratio)
            print(wrap_plots(glist))
            if (length(features) > 1) {
                if (plot_violin && length(levels(object_resized)) > 1)
                    print(plot_violin(object_resized, features, cells))
            }
            dev.off()
        }
        if(png) {
            png(gsub(".pdf$", ".png", filename), width = width, height = width * ratio, units = "in", res = 150)
            print(wrap_plots(glist))
            dev.off()
        }
    }
    object_resized
}

#' plot_spatially_variable
#'
#' Plots results of spatial differential expression analysis
#' @param object Seurat object 
#' @param labels Optional \code{character(n)} vector with labels
#' @param spatial_features Results of\code{Seurat::SpatiallyVariableFeatures}.
#' @param method Method used to find spatially variable features
#' @param number_features Plot that many features.
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param ... Arguments passed to \code{\link{plot_features}}
#' @importFrom gridExtra marrangeGrob
#' @export plot_spatially_variable
#' @examples
#' plot_spatially_variable
plot_spatially_variable <- function(object, labels = NULL, spatial_features, method = "markvariogram", 
    number_features = 80, prefix, subdir = "spatial_variation", ...) {
    top_features <- lapply(spatial_features, function(x) head(x, number_features))
    top_features <- lapply(top_features, .reorder_spatially_variable_features, object)
    object_split <- SplitObject(object, split.by = "library")
    libs <- as.character(sapply(object_split, function(x) x$library[1]))
    for (i in seq_along(object_split)) {
        flog.info("Generating output plots for %s ...", libs[i])
        label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
        libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
        plot_features(object = object_split[[i]], 
            features = top_features[[i]],
            prefix = prefix,
            suffix = paste0("_he_variable_", method, label, libs_label, ".pdf"),
            subdir = subdir,
            ...)
    }
}

.reorder_spatially_variable_features <- function(features, x) {
    m1 <- FetchData(x, features)
    idx <- colnames(x@meta.data)[grep("^nmf_k", colnames(x@meta.data))]
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


.resize_slice_images <- function(object, w = 300) {
    if (!requireNamespace("EBImage", quietly = TRUE)) return(object)
    .resize_image <- function(k) {
        new_k <- paste0(k, "_scaled")
        object@images[[new_k]] <- object@images[[k]]
        object@images[[new_k]]@image <- EBImage::resize(object@images[[k]]@image, w = w)
        r <- w / nrow(object@images[[k]]@image)
        object@images[[k]] <- NULL
        object@images[[new_k]]@scale.factors$lowres <- object@images[[new_k]]@scale.factors$lowres * r
        object
    }
    all_images <- Images(object)
    for(i in all_images) {
        object <- .resize_image(i)
    }    
    object
}

#' plot_qc_read
#'
#' Plot QC after reading a Spatial read counts
#' @param object Seurat object
#' @param prefix Prefix of output files
#' @param assay Name of the assay corresponding to the initial input data.
#' @export plot_qc_read
#' @examples 
#' plot_qc_read
plot_qc_read <- function(object, prefix, assay) {
    pdf(paste0(prefix, "_qc.pdf"))
    assays <- assay
    for (assay in assays) {
        print(VlnPlot(object = object,
            features = c(paste0("nFeature_", assay), 
                         paste0("nCount_", assay), 
                         "percent.mito", "percent.ribo"),
            ncol = 2))
        par(mfrow = c(2, 2))
        print(FeatureScatter(object = object, 
            feature1 = paste0("nFeature_", assay), 
            feature2 = paste0("nCount_", assay)))
        print(FeatureScatter(object = object, 
            feature1 = paste0("nFeature_", assay),, 
            feature2 = "percent.ribo"))
        if (length(unique(object$percent.mito))>1) { 
            print(FeatureScatter(object = object, 
                feature1 = paste0("nFeature_", assay),
                feature2 = "percent.mito"))
            print(FeatureScatter(object = object, 
                feature1 = "percent.mito", 
                feature2 = "percent.ribo"))

        }
    }
    if ("S.Score" %in% colnames(object@meta.data)) {
        par(mfrow = c(1, 1))
        print(VlnPlot(object = object,
            features = c("S.Score", "G2M.Score")))
        if (length(unique(object$percent.mito))>1) { 
            par(mfrow = c(2, 2))
            print(FeatureScatter(object = object, 
                feature1 = "S.Score", 
                feature2 = "percent.ribo"))
            print(FeatureScatter(object = object, 
                feature1 = "G2M.Score", 
                feature2 = "percent.ribo"))
            print(FeatureScatter(object = object, 
                feature1 = "S.Score", 
                feature2 = "percent.mito"))
            print(FeatureScatter(object = object, 
                feature1 = "G2M.Score", 
                feature2 = "percent.mito"))
        }
    }    
    dev.off()
    .write_qc_stats(object, prefix)
}    

#' plot_sc3
#'
#' Basic plots of saved \code{SC3} information
#' @param object Seurat object that contains \code{SC3} information after
#' running \code{\link{cluster_sc3}}
#' @param libs Library ids, must be stored in \code{object$library}
#' @param labels Optional \code{character(n)} with sample labels
#' @param rank Number of clusters 
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param width Output PDF width
#' @param pdf Create PDF image files
#' @param png Create PNG image files
#' @param plot_he Plot for each library the HE jpeg
#' @param assay Name of the assay corresponding to the initial input data.
#' @param ... Additional parameters passed to \code{\link{plot_features}}
#' @export plot_sc3
#' @examples
#' plot_sc3()
plot_sc3 <- function(object, libs, labels = NULL, rank, prefix, 
                     subdir = "sc3", width = 10, pdf = FALSE, png = FALSE,
                     plot_he = TRUE, assay = "Spatial", ...) {

    libs <- as.vector(libs)
    #create subdir if not exists
    .get_sub_path(prefix, subdir, "")

    for (k in rank) {
        object <- set_idents_sc3(object, k, rank)

        tmp <- .get_sub_path(prefix, file.path(subdir, "umap"), "") # make sure that umap directory exists
        filename <- .get_sub_path(prefix, file.path(subdir, "umap", k),
            paste0("_spatial_cluster.pdf"))
        pdf(filename, width = 10, height = 5)
        gp <- DimPlot(object, reduction = "umap", label = TRUE)
        if (requireNamespace("ggthemes", quietly = TRUE) &&
            length(levels(Idents(object))) <= 8) {
            gp <- gp + ggthemes::scale_colour_colorblind()
        }
        print(gp)
        dev.off()
        if ("label" %in% colnames(object@meta.data)) {
            .plot_cluster_library(object, field = "label", prefix = prefix, 
                subdir = file.path(subdir, "umap", k),
                reference_technology = "spatial", pdf = pdf, png = png)
        } else if ("library" %in% colnames(object@meta.data)) {
            .plot_cluster_library(object, field = "library", prefix = prefix, 
                subdir = file.path(subdir, "umap", k),
                reference_technology = "spatial", pdf = pdf, png = png)
        }

        ratio <- .get_image_ratio(k)
        if (plot_he) {
            flog.info("Generating output plots for k = %i...", k)
            tmp <- .get_sub_path(prefix, file.path(subdir, "he"), "") # make sure that HE directory exists

            for (i in seq_along(libs)) {
                label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
                libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
                obj_split <- object[,object$library == libs[i]]

                sd.plot <- SpatialDimPlot(obj_split, label = TRUE, images = 
                                          .get_image_slice(obj_split), label.size = 3, ...)

                if (requireNamespace("ggthemes", quietly = TRUE) &&
                    length(levels(Idents(obj_split))) <= 8) {
                    sd.plot <- sd.plot + ggthemes::scale_fill_colorblind()
                }    
                filename <- .get_sub_path(prefix, file.path(subdir, "he", k),
                    paste0("_he_sc3_discrete_cluster_", k, label, libs_label, ".pdf"))
                if (pdf) {
                    pdf(filename, width = 4, height = 3.9)
                    print(sd.plot)
                    invisible(dev.off())
                }
                if(png) {
                    png(sub(".pdf$", ".png", filename), width = 4, height = 3.9, units = "in", res = 150)
                    print(sd.plot)
                    invisible(dev.off())
                }
            }
        }
    }        
}        

#' plot_scoloc
#'
#' Co-localization plot for cell types from \code{celltrek} integration
#' @param celltrek_object \code{celltrek} object to plot co-localization for
#' @param directed Should the graph be directed or not (default: FALSE)
#' @export plot_scoloc
#' @examples
#' plot_scoloc(celltrek)
plot_scoloc <- function(celltrek_object, directed=FALSE) {
  plot_cell <- names(sort(table(celltrek_object$cell_type), decreasing=TRUE)[1:8])
  names(plot_cell) <- make.names(plot_cell)
  ct_plot <- subset(celltrek_object, subset=cell_type %in% plot_cell)
  ct_plot$cell_type <- factor(ct_plot$cell_type, levels=plot_cell)

  sgraph_KL <- CellTrek::scoloc(ct_plot, col_cell='cell_type', use_method='KL', eps=1e-50)

  ## We extract the minimum spanning tree (MST) result from the graph
  sgraph_KL_mst_cons <- sgraph_KL$mst_cons
  rownames(sgraph_KL_mst_cons) <- colnames(sgraph_KL_mst_cons) <- plot_cell[colnames(sgraph_KL_mst_cons)]

  ## We then extract the metadata (including cell types and their frequencies)
  cell_class <- ct_plot@meta.data %>% dplyr::select(id=cell_type) %>% unique
  ct_count <- data.frame(freq = table(celltrek_object$cell_type))
  cell_counts <- merge(cell_class, ct_count, by.x ="id", by.y = "freq.Var1")

  write.csv(sgraph_KL_mst_cons, paste0(fileprefix, '_sgraph_KL_mst_cons', label, ".csv"), row.names=FALSE)
  write.csv(cell_counts, paste0(fileprefix, '_cell_class', label, ".csv"), row.names=FALSE)

  mst_cons_am <- sgraph_KL_mst_cons

  if (!directed) mst_cons_am[upper.tri(mst_cons_am, diag = T)] <- NA

  mst_cons_am <- data.frame(id=rownames(mst_cons_am), mst_cons_am, check.names=F)
  mst_cons_table <- reshape2::melt(mst_cons_am) %>% na.omit() %>% magrittr::set_colnames(c('source', 'destination', 'weight'))
  mst_cons_table$destination <- as.character(mst_cons_table$destination)
  mst_cons_table <- tibble(mst_cons_table)

  node_col='id'
  node_size='freq.Freq'

  if (!is.null(cell_counts)) {
    if (node_col!='None') {
      froms <- mst_cons_table %>% distinct(source) %>% rename(label=source)
      tos <- mst_cons_table %>% distinct(destination) %>% rename(label=destination)

      mst_cons_node <- full_join(froms, tos, by="label")
      mst_cons_node <- mst_cons_node %>% mutate(id=1:nrow(mst_cons_node)) %>% select(id, everything())

      mst_cons_edge <- mst_cons_table %>% left_join(mst_cons_node, by=c("source"="label")) %>% rename(from=id)
      mst_cons_edge <- mst_cons_edge %>% left_join(mst_cons_node, by=c("destination"="label")) %>% rename(to=id)
      mst_cons_edge <- select(mst_cons_edge, from, to, weight)
      mst_cons_edge <- mst_cons_edge[which(mst_cons_edge$weight>0),]


      node_cols <- ggpubr::get_palette('Set1', length(unique(mst_cons_node$label)))
      names(node_cols) <- unique(mst_cons_node$label)
      mst_cons_node$color <- node_cols[match(mst_cons_node$label, names(node_cols))]
    }
    if (node_size!='None') {
      size_colmn <- as.character(node_size)
      size_df <- data.frame(id=cell_counts$id, value=as.numeric(cell_counts[, size_colmn])/sum(as.numeric(cell_counts[, size_colmn])))
      mst_cons_node$size <- size_df$value[match(mst_cons_node$label, size_df$id)]
    }
      library(ggraph)
      library(tidygraph)
      net.tidy <- tbl_graph(nodes=mst_cons_node, edges=mst_cons_edge, directed=directed)
      ggraph(net.tidy, layout = "graphopt") + geom_node_point(aes(colour=color, size=size)) +
          geom_edge_link(aes(width = weight), alpha = 0.5) + scale_edge_width(range = c(0.5, 5)) +
          scale_size_continuous(range=c(5,25)) + geom_node_text(aes(label = label), repel = TRUE) 
  }
}



#' plot_predictions
#'
#' Plots predictions in a Seurat object on HE slides. Supports multi-page
#' plot when lots of features are requested.
#' @param object Seurat object
#' @param predictions object with transfer predictions
#' @param ignore Features to be ignored
#' @param label Label of single cell dataset
#' @param label_integration_method Label of integration method
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param pdf Create PDF image files
#' @param png Create PNG image files
#' @param ... Arguments passed to \code{plot_features}.
#' @export plot_predictions
#' @examples
#' plot_predictions()
plot_predictions <- function(object, predictions = NULL, ignore = c("max", "unassigned"),
                          label = NULL, label_integration_method = "default", prefix, 
                          subdir = "he", pdf = FALSE, png = FALSE,
                          ...) {
    if (is.null(predictions) && is.null(object[["predictions"]])) {
        stop("object needs a predictions object if not provided.")
    }
    if (!is.null(predictions)) {
        object$predictions <- predictions
        DefaultAssay(object) <- "predictions"
    }
    Idents(object) <- GetTransferPredictions(object)
    features <- names(Matrix::rowSums(GetAssayData(object)) > 0)
    features <- features[!features %in% ignore]
    label <- if (is.null(label)) "" else paste0("_", label)
    flog.info("Generating output plots for %s ...", label)
    if (length(Images(object)) > 1 && "library" %in% colnames(object@meta.data)) {
        x_split <- SplitObject(object, split.by = "library")
        libs <- sapply(x_split, function(y) y$library[1])
        libs_label <- rep("", length(libs))
        field <- "library"
        if ("label" %in% colnames(object@meta.data)) {
            libs_label <- paste0("_", sapply(x_split, function(y) y$label[1]))
            field <- "label"
        }
        for (j in seq_along(libs)) {
            plot_features(object = x_split[[j]], features = features,
                prefix = prefix, subdir = subdir,
                suffix = paste0("_he_labels", label, "_", libs[j], libs_label[j], "_", label_integration_method,".pdf"),
                pdf = pdf, png = png,
                ...)
            filename <- .get_sub_path(prefix, "he",
                    suffix = paste0("_he_labels_call", label, "_", libs[j], libs_label[j], "_", label_integration_method, ".pdf"))
            gp <- SpatialDimPlot(x_split[[j]], label = TRUE,
                images = .get_image_slice(x_split[[j]]),
                label.size = 3)
            if (pdf) {
                pdf(filename, width = 4, height = 3.9)
                print(gp)
                invisible(dev.off())
            }
            if (png) {
                png(gsub("pdf$", "png", filename), width = 4,
                    height = 3.9, units = "in", res = 150)
                print(gp)
                invisible(dev.off())
            }
        }
        filename <- .get_sub_path(prefix, "advanced",
                suffix = paste0("_labels", label, "_", libs[j], libs_label[j], "_", label_integration_method, ".pdf"))
        ratio <- .get_image_ratio(min(6, length(features)))
        glist <- VlnPlot(object, features = features, group.by = field, pt.size = 0.25, combine = FALSE)
        glist <- lapply(glist, function(p) ggplotGrob(p + theme(legend.position = "none")))
        if (length(features) > 6) {
            glist <- gridExtra::marrangeGrob(glist, ncol = 3, nrow = 2)
        } else {
            glist <- patchwork::wrap_plots(glist)
        }
        if (pdf) {
            ggsave(filename, glist,
                   width = 10, height = 10 * ratio)
        }
    } else {
        plot_features(object = object, features = features,
            prefix = prefix, subdir = subdir,
            suffix = paste0("_he_labels", label, "_", label_integration_method, ".pdf"),
            pdf = pdf, png = png, ...)
    }
}

#' plot_prediction_summary
#'
#' Plots a categorical prediction summary in a Seurat object on HE slides. 
#' @param object Seurat object
#' @param feature Name of the feature that contains the prediction summary
#' @param label Label of single cell dataset
#' @param label_integration_method Label of integration method
#' @param prefix Output file prefix
#' @param subdir Put files in a subdirectory
#' @param pdf Create PDF image files
#' @param png Create PNG image files
#' @param ... Arguments passed to \code{Seurat::SpatialDimPlot}.
#' @export plot_prediction_summary
#' @examples
#' plot_prediction_summary()
plot_prediction_summary <- function(object, feature = "spot_class",
                          label = NULL, label_integration_method = "default", prefix, 
                          subdir = "advanced", pdf = FALSE, png = FALSE,
                          ...) {
    Idents(object) <- FetchData(object, feature)[[feature]]

    label <- if (is.null(label)) "" else paste0("_", label)
    flog.info("Generating output plots for %s ...", label)
    if (length(Images(object)) > 1 && "library" %in% colnames(object@meta.data)) {
        x_split <- SplitObject(object, split.by = "library")
        libs <- sapply(x_split, function(y) y$library[1])
        libs_label <- rep("", length(libs))
        field <- "library"
        if ("label" %in% colnames(object@meta.data)) {
            libs_label <- paste0("_", sapply(x_split, function(y) y$label[1]))
            field <- "label"
        }
        for (j in seq_along(libs)) {
            filename <- .get_sub_path(prefix, subdir,
                    suffix = paste0("_he_labels_summary", label, "_", libs[j], libs_label[j], "_", label_integration_method, ".pdf"))
            gp <- SpatialDimPlot(x_split[[j]], label = TRUE,
                images = .get_image_slice(x_split[[j]]),
                label.size = 3, ...)
            if (requireNamespace("ggthemes", quietly = TRUE) &&
                length(levels(Idents(object))) <= 8) {
                gp <- gp + ggthemes::scale_fill_colorblind()
            }
            gp2 <- DimPlot(x_split[[j]], reduction = .get_reduction(x_split[[j]]), label = FALSE)
            if (requireNamespace("ggthemes", quietly = TRUE) &&
                length(levels(Idents(object))) <= 8) {
                gp2 <- gp2 + ggthemes::scale_colour_colorblind()
            }
            if (pdf) {
                pdf(filename, width = 7, height = 3.9)
                print(gp)
                print(gp2)
                invisible(dev.off())
            }
            if (png) {
                filename <- gsub("pdf$", "png", filename)
                png(filename, width = 4, height = 3.9, units = "in", res = 150)
                print(gp)
                invisible(dev.off())
                filename <- gsub("_he_labels_summary", "_umap_labels_summary", filename)
                png(filename, width = 7, height = 3.9, units = "in", res = 150)
                print(gp2)
                invisible(dev.off())
            }
        }
    } else {
        filename <- .get_sub_path(prefix, subdir,
                suffix = paste0("_he_labels_summary", label, "_", label_integration_method, ".pdf"))
        gp <- SpatialDimPlot(object, label = TRUE, label.size = 3, ...)
        if (requireNamespace("ggthemes", quietly = TRUE) &&
            length(levels(Idents(object))) <= 8) {
            gp <- gp + ggthemes::scale_fill_colorblind()
        }
        gp2 <- DimPlot(object, reduction = .get_reduction(object), label = FALSE)
        if (requireNamespace("ggthemes", quietly = TRUE) &&
            length(levels(Idents(object))) <= 8) {
            gp2 <- gp2 + ggthemes::scale_colour_colorblind()
        }
        if (pdf) {
            pdf(filename, width = 7, height = 3.9)
            print(gp)
            print(gp2)
            invisible(dev.off())
        }
        if (png) {
            filename <- gsub("pdf$", "png", filename)
            png(filename, width = 4, height = 3.9, units = "in", res = 150)
            print(gp)
            invisible(dev.off())
            filename <- gsub("_he_labels_summary", "_umap_labels_summary", filename)
            png(filename, width = 7, height = 3.9, units = "in", res = 150)
            print(gp2)
            invisible(dev.off())
        }
    }
}
