#' integrate_spatial
#'
#' Integrate SpatialTranscriptomics data with (matched) scRNA-seq data
#' @param obj_spatial Seurat object containing SpatialTranscriptomics data.
#' Can be \code{NULL}; in that case only references will be integrated.
#' @param references List of Seurat reference datasets 
#' @param features Features to be used for integration
#' @param min_features Remove cells or spots with small number of detected genes
#' @param min_spots Remove or merge small samples with fewer spots
#' @param scale Scale features 
#' @param max_percent_mito Remove cells or spots with high mitochondria expression.
#' @param min_max_counts Remove cells or spots with low maximum counts.
#' @param min_detected Remove features detected in fewer than \code{min_detected}
#' \code{references} (including \code{obj_spatial}). 
#' @param reference_technology Used to fill the \code{meta.data@technology} field
#' @param skip_alternative_batch_corrections By default, harmony and fastmnn will
#' be run when installed. Provide methods to skip here as array.
#' @param force Recalculate, even when serialized objects are available
#' @param serialize Serialize output objects
#' @param prefix Prefix of output files
#' @param verbose Verbose Seurat output
#' @export integrate_spatial
#' @examples
#' integrate_spatial()
integrate_spatial <- function(obj_spatial = NULL, references, features = 2000,
                            min_features = 400, min_spots = 200, scale = TRUE, 
                            max_percent_mito = 33, 
                            min_max_counts = 3,
                            min_detected = 2,
                            reference_technology = "single_cell",
                            skip_alternative_batch_corrections = c(),
                            force, 
                            serialize = TRUE, prefix, verbose = FALSE) {
    if (serialize || !force) { 
        filename <- .get_serialize_path(prefix, 
            paste0("_", reference_technology, "_integrated.rds"))
        if (!dir.exists(dirname(filename))) dir.create(dirname(filename))
    }
    if (!force && file.exists(filename)) {
        flog.warn("%s exists. Skipping alignment. Use --force to overwrite.", filename)
        return(readRDS(filename))
    }
    references <- lapply(references, function(x) { 
        x[["technology"]] <- reference_technology
        x[["reference"]] <- TRUE
        .filter_object(x, max_percent_mito, min_features, min_max_counts)

    })
    if (!is.null(obj_spatial)) {
        obj_spatial[["technology"]] <- "spatial"
        obj_spatial[["reference"]] <- FALSE
        obj_spatial <- .filter_object(obj_spatial, max_percent_mito, min_features, 0)
        references <- c(list(obj_spatial), references)
    }
    references <- .remove_undetected_features(references, min_detected)

    all_features <- Reduce(intersect, lapply(references, rownames))

    references <- lapply(references, function(x) x[rownames(x) %in% all_features,])

    
    flog.info("Integrating libraries...")

    if (is.null(features)) {
        flog.warn("No features provided, using number of variable features of first reference.")
        features <- length(VariableFeatures(references[[1]]))
    } else {
        if (length(features) > 1 && any(!features %in% all_features)) {
            flog.warn("Some requested features not available: %s",
                paste(features[!features %in% all_features], collapse = ","))
            features <- features[features %in% all_features]
        }
    }
    num_spots <- sapply(references, ncol)
    .plot_pre_integration_umap(references, scale, features = features,
        skip_alternative_batch_corrections = skip_alternative_batch_corrections,
        force = force, prefix = prefix)
    if (min(num_spots) < min_spots) {
        poor_libs <- sapply(references[num_spots < min_spots], function(x) x$library[1])
        flog.warn("Samples with low number of cells/spots: %s",
            paste(poor_libs, collapse = ", "))
        if (length(poor_libs) > 1) references <- .merge_poor_libs(references, min_spots)
    } else {
        flog.info("Minimum number of cells/spots: %i (%s)", min(num_spots),
            references[[which.min(num_spots)]]$library[1])
    }
    if (references[[1]]@active.assay == "SCT") {
        integrated <- .integrate_sct(references, features, scale, verbose)
    } else {
        integrated <- .integrate_default(references, features, scale, verbose)
    }
    DefaultAssay(object = integrated) <- "integrated"
    if (serialize) {
        flog.info("Writing R data structure to %s...", filename)
        saveRDS(integrated, filename)
        flog.info("Done writing R data structure.")
    }
    integrated
}

.integrate_default <- function(references, features, scale, verbose) {
    flog.info("Using default integration...")
    anchors <- FindIntegrationAnchors(object.list = references, dims = 1:30, 
        anchor.features = features, scale = scale, verbose = verbose)
    IntegrateData(anchorset = anchors, dims = 1:30, verbose = verbose)
}

.integrate_sct <- function(references, features, scale, verbose) {
    if (scale) {
        flog.warn("Not scaling because references are SCTransform normalized.")
        scale = FALSE
    }     
    flog.info("Using SCT integration...")
    references <- PrepSCTIntegration(object.list = references, 
        anchor.features = features, verbose = verbose)
    k_filter_max <- min(sapply(references, ncol))
    k_filter <- 200
    if (k_filter_max < k_filter) {
        k_filter <- k_filter_max - 5
        flog.warn("Have to lower the number of Nearest Neighbors to %i.", k_filter)
    }
    anchors <- FindIntegrationAnchors(object.list = references, dims = 1:30, 
        anchor.features = features, normalization.method = "SCT", k.filter = k_filter,
        verbose = verbose)
    integrated <- IntegrateData(anchorset = anchors, 
        normalization.method = "SCT", dims = 1:30, verbose = verbose)
}
    
.merge_poor_libs <- function(references, min_spots) {
    num_spots <- sapply(references, ncol)
    references <- c(references[num_spots >= min_spots], Reduce(merge, references[num_spots < min_spots]))
    references[order(sapply(references, function(x) x$library[1]))]
}

.filter_object <- function(x, max_percent_mito, min_features, min_max_counts, assay = "Spatial") {
    nbefore <- ncol(x)
    if ("percent.mito" %in% colnames(x@meta.data)) {
        x <- x[, is.na(x[["percent.mito"]]) | 
                 x[["percent.mito"]] <= max_percent_mito]
    }
    library_label <- paste("library", x$library[1])
    if (length(unique(x$library)) > 1) library_label <- "libraries"
    if (ncol(x) < nbefore) {
        flog.info("Percent Mito: Removing %i cells from reference %s.", 
            nbefore - ncol(x), library_label)
    }
    nbefore <- ncol(x)
    x <- x[ , x[[paste0("nFeature_", assay)]] >= min_features]
    if (ncol(x) < nbefore) {
        flog.info("Min Genes: Removing %i cells from reference library %s.", 
            nbefore - ncol(x), x$library[1])
    }
    if (min_max_counts > 0) {
        nbefore <- nrow(x) 
        max_counts <- apply(GetAssayData(x, slot = "counts"),1, max) 
        x <- x[names(which(max_counts >= min_max_counts)),]
        if (nrow(x) < nbefore) {
            flog.info("Max Counts: Removing %i features from reference library %s.", 
                nbefore - ncol(x), x$library[1])
        }
    }
    x
}

#' cluster_integrated
#'
#' Cluster integrated data (obtained by \code{\link{integrate_spatial}}.
#' @param integrated Seurat object containing integrated 
#' SpatialTranscriptomics and scRNA-seq data
#' @param regressout Regressout these features.
#' @param scale Scale data. If \code{NULL} use recommendation.
#' @param force Recalculate, even when serialized objects are available
#' @param plot_umap Generate UMAP plot of integrated data
#' @param serialize Serialize output objects
#' @param prefix Prefix of output files
#' @param umap_suffix File suffix of UMAP PDF plot
#' @param verbose Verbose Seurat output
#' @export cluster_integrated
#' @examples
#' cluster_integrated()
cluster_integrated <- function(integrated, regressout, scale = NULL,
                               force, plot_umap = TRUE,
                               serialize = TRUE, prefix, 
                               umap_suffix = "_umap_integration_overview.pdf",
                               verbose = FALSE) {
    reference_technology <- integrated$technology[integrated$reference][1]
    if (serialize || !force) {
        filename <- .get_serialize_path(prefix,
            paste0("_", reference_technology, "_integrated_scaled.rds"))
    }
    if (!force && file.exists(filename)) {
        flog.warn("%s exists. Skipping scaling. Use --force to overwrite.",
            filename)
        return(readRDS(filename))
    }
    if (is.null(scale)) scale <- ! "SCT" %in% names(integrated@assays)
    if (scale) {   
        flog.info("Scaling data...")
        regressout <- .check_regressout(integrated, regressout)
        integrated <- ScaleData(object = integrated, verbose = verbose,
            vars.to.regress = regressout)
    } else {
        flog.info("Skipping scaling.")
    }
    flog.info("Running PCA...")
    integrated <- RunPCA(object = integrated, npcs = 30, verbose = verbose)
    flog.info("Running UMAP...")
    integrated <- RunUMAP(object = integrated, reduction = .get_reduction(integrated, type = "pca"),
        dims = 1:30, verbose = verbose)
    if (plot_umap) {
        flog.info("Plotting UMAP...")
        pdf(.get_advanced_path(prefix, umap_suffix), width = 10, height = 5)
        if ("call" %in% colnames(integrated@meta.data)) {
            print(DimPlot(integrated, reduction = .get_reduction(integrated), split.by = "technology",
                          group.by = "call"))
        } else {
            print(DimPlot(integrated, reduction = .get_reduction(integrated), split.by = "technology",
                          group.by = "library"))
        }    
        dev.off()
    }
    if (serialize) {
        flog.info("Writing R data structure to %s...", filename)
        saveRDS(integrated, filename)
    }
    integrated
}

#' cluster_reference
#'
#' Cluster the scRNA-seq reference in the integrated data 
#' (obtained by \code{\link{integrate_spatial}}.
#' @param integrated Seurat object containing integrated 
#' SpatialTranscriptomics and scRNA-seq data
#' @param resolution Cluster resolution
#' @param sub Flag indicating a subset of the data was provided
#' @param idents Use idents from a serialized Seurat object or as provided
#' @param features Use only the specified features for clustering
#' @param force Recalculate, even when serialized objects are available
#' @param plot_umap Output UMAP plot 
#' @param serialize Serialize output objects
#' @param prefix Prefix of output files
#' @param verbose Verbose Seurat output
#' @param ... Additional parameters passed to \code{Seurat::FindClusters}
#' @export cluster_reference
#' @examples
#' cluster_reference()
cluster_reference <- function(integrated, resolution = 0.8, sub = FALSE,
                         idents = NULL, features = NULL, 
                         force, plot_umap = TRUE, serialize = TRUE, prefix, 
                         verbose = FALSE, ...) {
    sub <- ifelse(sub, "sub_","")
    features <- ifelse(!is.null(features), "features_", "")
    reference_technology <- integrated$technology[integrated$reference][1]
    resolution_label <-  paste0(resolution, collapse="_")
    if (serialize || !force) {
        filename <- .get_serialize_path(prefix, paste0("_", reference_technology,
            "_cluster_", sub, features, resolution_label, ".rds"))
    }
    if (!force && file.exists(filename)) {
        flog.warn("%s exists. Skipping scRNA clustering. Use --force to overwrite.", filename)
        obj_ref <- readRDS(filename)
    } else {
        obj_ref <- integrated[, integrated$reference]
        if (!is.null(idents)) {
            flog.info("Loading existing idents from %s...", idents)
            if (is(idents, "character")) {
                idents <- Idents(readRDS(idents))
            }
            Idents(obj_ref) <- idents
            flog.info("Found %i idents.", length(levels(Idents(obj_ref))))
        } else {
            flog.info("Finding neighbors of integrated reference data...")
            obj_ref <- FindNeighbors(obj_ref, verbose = verbose, features = features)
            flog.info("Finding clusters of integrated reference data...")
            obj_ref <- FindClusters(obj_ref, resolution = resolution, verbose = verbose, ...)
            flog.info("Generating UMAP plots of integrated reference data...")
            label <- "predicted.id" %in% colnames(obj_ref@meta.data)
            if (label) {
               ids <- sapply(levels(obj_ref), .get_ident_label, obj_ref)
               names(x = ids) <- levels(x = obj_ref)
               obj_ref <- RenameIdents(object = obj_ref, ids)
               if (plot_umap) {
                   pdf(.get_advanced_path(prefix, paste0("_", reference_technology, "_cluster_", sub, features,
                        resolution_label, "predictedid.pdf")), width = 10, height = 10)
                   print(DimPlot(obj_ref, reduction = .get_reduction(obj_ref), split.by = "predicted.id"))
                   dev.off()
               }
            }
        }
        obj_ref[["new.idents"]] <-  Idents(obj_ref)
        if (serialize) {
            flog.info("Writing R data structure to %s...", filename)
            saveRDS(obj_ref, filename)
        }
    }
    obj_ref
}    


#' find_nearest_neighbors
#'
#' Shows the correlation of spots across multiple samples,
#' after integration
#' @param object Seurat object containing integrated 
#' spatial data
#' @param split.by Split \code{object} by this feature.
#' @export find_nearest_neighbors
#' @examples
#' find_nearest_neighbors()
find_nearest_neighbors <- function(object, split.by = "library") {
    x <- SplitObject(object, split.by = split.by)
    if (length(x) < 2) {
        flog.info("Cannot find multiple samples in this object.")
        return(object)
    }    
    idx_ref <- which.max(sapply(x, ncol))
    m_ref <- Matrix::as.matrix(GetAssayData(x[[idx_ref]]))
    pos <- x[[idx_ref]]@images[[.get_image_slice(x[[idx_ref]])]]@coordinates[,c("row", "col")]
    x[[idx_ref]]$nn <- rank(apply(pos^2, 1, sum))
    for (i in seq_along(x)[-idx_ref]) {
        m <- Matrix::as.matrix(GetAssayData(x[[i]]))
        nn <- apply(m,2,function(x) which.max(cor(x, m_ref)))
        x[[i]]$nn <- x[[idx_ref]]$nn[nn]
    }

    md <- do.call(c, c(lapply(x, function(y) y$nn), list(use.names=FALSE)))
    names(md) <- Cells(object)

    AddMetaData(object,
        metadata = md,
        col.name = "int.nearest.neighbor"
    )
}

.get_ident_label <- function(i, obj_ref) {
    md <- obj_ref@meta.data[which(Idents(obj_ref) == i),]
    label <- names(sort(table(md$predicted.id[ md$prediction.score.max>0.9]),decreasing=TRUE)[1])
    if (is.null(label)) { 
        call <- names(sort(table(md$call), decreasing=TRUE)[1])
        label <- paste0("unk_", call)
    }
    paste0(label, "_", i)    
}

.merge_safely <- function(obj) {
    if (any(duplicated(unlist(sapply(obj, Cells))))) {
        return(
            Reduce(merge, lapply(seq_along(obj), function(i) {
                RenameCells(obj[[i]], add.cell.id = i)}))
        )
    }
    Reduce(merge, obj)
}

.remove_undetected_features <- function(references, min_detected) {
    # find number of undetected references for all features undetected
    # in at least one reference
    tbl <- table(unlist(lapply(references, function(x) { 
        y <- Matrix::rowSums(GetAssayData(x, slot = "counts"))
        names(y[which(y <= 0)])
        })))
    tbl_detected <- length(references) - tbl
    undetected_features <- names(tbl_detected[tbl_detected < min_detected])
    if (!length(undetected_features)) return(references)
    flog.warn("Removing features detected in fewer than %i references: %s", 
        min_detected, paste(undetected_features, collapse = ","))    
    lapply(references, function(x) x[!rownames(x) %in% undetected_features, ])
}
    

.merge_by_harmony <- function(references, force, prefix) {
    if (is(references, "Seurat")) {
        merged <- references
    } else {
        merged <- Reduce(merge, references)
        VariableFeatures(merged) <- features
    }
    if (requireNamespace("harmony", quietly = TRUE)) {
        filename <- .get_serialize_path(prefix, "_harmony.rds")
        if (!force && file.exists(filename)) {
            flog.warn("%s exists. Skipping Harmony. Use --force to overwrite.", filename)
            merged <- readRDS(filename)
        } else {
            assay.use <- "SCT"
            if (!assay.use %in% Assays(merged)) assay.use <- "Spatial"
            if (!assay.use %in% Assays(merged)) assay.use <- "RNA"
            flog.info("Running Harmony on assay %s...", assay.use)
            merged <- RunPCA(object = merged, npcs = 30, verbose = TRUE,
                assay = assay.use)
            merged <- harmony::RunHarmony(merged, "library", assay.use = assay.use)
            merged <- RunUMAP(object = merged, assay = assay.use,
                reduction = "harmony", dims = 1:30)
            flog.info("Writing R data structure to %s...", filename)
            saveRDS(merged, filename)
        }
        flog.info("Plotting UMAP...")
        umap_suffix <- "_umap_post_harmony_overview.pdf"
        pdf(.get_advanced_path(prefix, umap_suffix), width = 10, height = 5)
        print(DimPlot(merged, reduction = .get_reduction(merged), split.by = "technology",
                          group.by = "library"))
        dev.off()
    } else {
        flog.warn("Install harmony package for Harmony integration.")
    }
    merged
}        
.merge_by_fastmnn <- function(references, force, prefix, features = 3000, cluster = FALSE) {
    if (requireNamespace("SeuratWrappers", quietly = TRUE) &&
        requireNamespace("batchelor", quietly = TRUE)) {
        assay.use <- "SCT"
        if (!assay.use %in% Assays(references[[1]])) assay.use <- "Spatial"
        if (!assay.use %in% Assays(references[[1]])) assay.use <- "RNA"
        filename <- .get_serialize_path(prefix, "_fastmnn.rds")
        if (!force && file.exists(filename)) {
            flog.warn("%s exists. Skipping fastMNN. Use --force to overwrite.", filename)
            merged <- readRDS(filename)
        } else {
            flog.info("Running fastMNN on assay %s...", assay.use)
            merged <- SeuratWrappers::RunFastMNN(Seurat:::CheckDuplicateCellNames(references),
                assay = assay.use, features = features)
            VariableFeatures(merged) <- features
            merged <- RunUMAP(merged, reduction = "mnn", dims = 1:30, assay = assay.use)
            if (cluster) {
                merged <- FindNeighbors(merged, reduction = "mnn", dims = 1:30, assay = assay.use)
                merged <- FindClusters(merged, assay = assay.use)
            }
            flog.info("Writing R data structure to %s...", filename)
            saveRDS(merged, filename)
        }
        flog.info("Plotting UMAP...")
        umap_suffix <- "_umap_post_fastmnn_overview.pdf"
        pdf(.get_advanced_path(prefix, umap_suffix), width = 10, height = 5)
        DefaultAssay(merged) <- assay.use
        print(DimPlot(merged, reduction = .get_reduction(merged), split.by = "technology",
                          group.by = "library"))
        dev.off()
    } else {
        flog.warn("Install SeuratWrappers and batchelor packages for fastMNN integration.")
    }
}        
.plot_pre_integration_umap <- function(references, scale, features,
    skip_alternative_batch_corrections, force, prefix) {
    flog.info("Plotting pre-integration UMAP...")
    merged <- Reduce(merge, references)
    VariableFeatures(merged) <- features
    if (! "harmony" %in% skip_alternative_batch_corrections) {
        invisible(.merge_by_harmony(merged, force, prefix))
    } else {
        flog.info("Skipping harmony as requested.")
    }    
    if (! "fastmnn" %in% skip_alternative_batch_corrections) {
        invisible(.merge_by_fastmnn(references, force, prefix, features))
    } else {
        flog.info("Skipping fastmnn as requested.")
    }
    merged <- cluster_integrated(merged, regressout = NULL,
                               force = TRUE, plot_umap = TRUE,
                               scale = scale,
                               serialize = FALSE, prefix,
                               umap_suffix = "_umap_pre_integration_overview.pdf",
                               verbose = FALSE)
    return(1)
}

.harmonize_assayobjects <- function(x, ignore_unassigned = TRUE) {
    # first make sure all have the same set of features
    features <- Reduce(union, lapply(x, function(y) names(which(rowSums(y) > 0))))
    features <- features[features != "max"]
    if (ignore_unassigned) {
        features <- features[features != "unassigned"]
    }
    x <- lapply(x, function(y) {
        missing <- features[!features %in% rownames(y)]
        if (!length(missing)) return(y)
        m <- GetAssayData(y)
        m_missing <- matrix(0, nrow = length(missing), ncol = ncol(m))
        rownames(m_missing) <- missing
        if (is(m, "dgCMatrix")) {
            m_missing <- as(m_missing, "dgCMatrix")
        }
        m <- rbind(m, m_missing)
        CreateAssayObject(m)
   })
   return(x)  
}

#' find_assayobject_consensus
#'
#' Averages multiple cell type prediction into one AssayObject
#' @param x list of \code{AssayObject}s
#' @param min_fraction Set all probabilities smaller this value to 0
#' @param drop_zero Drop cell types never assigned
#' @param ignore_unassigned Ignore the unassigned category (not all methods
#' support them - currently only rctd_multi)
#' @param labels Labels of \code{x}, usually the method identifier
#' @param plot_correlations Plot heatmap of correlations across \code{x} 
#' @param plot_cor_method Method used in the \code{cor} function
#' @export find_assayobject_consensus
#' @examples
#' find_assayobject_consensus()
find_assayobject_consensus <- function(x, min_fraction = 0.05, drop_zero = TRUE,
    ignore_unassigned = TRUE, labels = names(x), plot_correlations = FALSE,
    plot_cor_method = "pearson") {
    names(x) <- labels
    x <- .harmonize_assayobjects(x, ignore_unassigned = ignore_unassigned)
    features <- Reduce(intersect, lapply(x, rownames))
    features <- features[features != "max"]
    if (plot_correlations) {
        print(.plot_assayobject_consensus(x, features, cor_method = plot_cor_method))
    }     
    m <- Reduce("+", lapply(x, function(y) GetAssayData(y)[features,] / length(x)))
    m <- sweep(m, 2, Matrix::colSums(m), "/")
    m[m < min_fraction] <- 0
    m <- sweep(m, 2, Matrix::colSums(m), "/")
    idx_na_cols <- is.na(Matrix::colSums(m))
    if (any(idx_na_cols)) {
         m[, idx_na_cols] <- 0
    }
    if (drop_zero) {
         m <- m[Matrix::rowSums(m) > 0, ]
    }    
    m <- rbind(m, max = apply(m, 2, max))
    m <- as(m, "sparseMatrix")
    CreateAssayObject(data = m)
}

#' combine_assayobjects
#'
#' Combines multiple cell type prediction into one AssayObject, keeps all predictions
#' @param x list of \code{AssayObject}s
#' @param ignore_unassigned Ignore the unassigned category (not all methods
#' support them - currently only rctd_multi)
#' @param labels Labels of \code{x}, usually the method identifier
#' @export combine_assayobjects
#' @examples
#' combine_assayobjects()
combine_assayobjects <- function(x, ignore_unassigned = TRUE, labels = names(x)) {
    if (is.null(labels)) {
        stop("need labels")
    }    
    names(x) <- labels
    x <- .harmonize_assayobjects(x, ignore_unassigned = ignore_unassigned)
    features <- head(rownames(x[[1]]), -1)

    xl <- lapply(x, function(y) GetAssayData(y)[features,])
    xl <- lapply(seq_along(xl), function(i) {
        rownames(xl[[i]]) <- paste(labels[i], make.names(rownames(xl[[i]])), sep = ".")
        return(xl[[i]])
    })
    m <- do.call(rbind, xl)
    rownames(m) <- gsub("_", "-", rownames(m))
    m <- as(m, "sparseMatrix")
    CreateAssayObject(data = m)
}

.plot_assayobject_consensus <- function(x, features, cor_method = "pearson") {
    m_cor <- lapply(features, function(f) suppressWarnings(cor(do.call(cbind,
        lapply(x, function(y) GetAssayData(y)[f, ])), method = cor_method)))
    names(m_cor) <- features
    dt_cor <- rbindlist(lapply(names(m_cor), function(i) melt(data.table(
        cell_type = i,
        method = rownames(m_cor[[i]]),
        m_cor[[i]]), id.vars = c("cell_type", "method"))))
    dt_cor <- dt_cor[dt_cor$cell_type != "unassigned",]
    cor_method_label <- paste(cor_method, "Cor")
    substr(cor_method_label, 1, 1) <- toupper(substr(cor_method_label, 1, 1))
    gp <- ggplot(dt_cor, aes(method, variable, fill = value)) +
        geom_tile() +
        facet_wrap(~cell_type) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        xlab("") + ylab("") + labs(fill = cor_method_label) +
        scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(7,"RdBu")), limits = c(-1, 1))
    return(gp)    
}

#'classify_cell_type_regions
#'
#' Used to classify cell-types as core, periphery, diffused, other
#' after integration
#' @param object Seurat object containing integrated
#' @param cell_type Cell-type of interest
#' @param cutoff Deconvolution cutoff for core membership
#' @export classify_cell_type_regions
#' @examples
#' classify_cell_type_regions()
classify_cell_type_regions <- function(object, cell_type = "Tumor", cutoff = NULL) {
    ct_score <- object[["predictions"]][cell_type, ]
    if (is.null(cutoff)) {
        if (!requireNamespace("mclust", quietly = TRUE)) {
            flog.info("Mclust not found and cutoff not specified, will use 0.5 as cutoff.")
            cutoff <- 0.5
        } else {
            p <- mclust::Mclust(ct_score, G = 3, verbose = FALSE)
            cutoff <- sapply(split(ct_score, p$classification), min)[[2]]
        }
    }
    giotto_object <- suppressMessages(as_GiottoObject(object))
    giotto_object <- suppressMessages(createSpatialNetwork(giotto_object))
    net <- giotto_object@spatial_network$cell$Delaunay_network@networkDT
    net$from_ct <- ct_score[1, as.character(net$from)]
    net$to_ct <- ct_score[1, as.character(net$to)]
    fraction_nn_above_cutoff <- sapply(colnames(ct_score), function(ct) {
        idx <- net$from == ct | net$to == ct
        sum(net[idx, "from_ct"] > cutoff & net[idx, "to_ct"] > cutoff) / sum(idx)
    })
    core_candidates <- names(which(fraction_nn_above_cutoff > 0.75))
    core_members <-  sapply(colnames(ct_score), function(ct) {
        idx <- net$from == ct | net$to == ct
        if (!any(idx)) return(FALSE)
        other_ct <- apply(net[idx,,drop = FALSE], 1, function(x) {
            col_ct_id <- "from"
            col_other_id <- "to"
            if (x["to"] == ct) {
                col_ct_id <- "to"
                col_other_id <- "from"
            }
            r1 <- as.numeric(x[paste0(col_other_id, "_ct")])
            r2 <- x[col_other_id] %in% core_candidates
            return(c(r1, r2))
        })
        # remove cases where region contains only one spot
        if (ct %in% core_candidates && sum(other_ct[2, ]) > 1) return(TRUE)
        if (ct %in% core_candidates && sum(other_ct[2, ]) <= 1) return(FALSE)
        return(all(other_ct[1, ] > cutoff))
    })
    core_members <- names(which(core_members))

    periphery_members <- sapply(colnames(ct_score), function(ct) {
        if (ct %in% core_members) return(FALSE)

        idx <- net$from == ct | net$to == ct
        if (!any(idx)) return(FALSE)
        other_ct <- apply(net[idx,,drop = FALSE], 1, function(x) {
            col_ct_id <- "from"
            col_other_id <- "to"
            if (x["to"] == ct) {
                col_ct_id <- "to"
                col_other_id <- "from"
            }
            r1 <- as.numeric(x[paste0(col_other_id, "_ct")])
            r2 <- x[col_other_id] %in% core_members
            return(c(r1, r2))
        })
        if (ct_score[1, ct] > 0.05 && sum(other_ct[2,]) > 0) return(TRUE)
        return(FALSE)
    })
    periphery_members <- names(which(periphery_members))
    
    diffuse_members <- sapply(colnames(ct_score), function(ct) {
        return(ct_score[1, ct] > cutoff && !ct %in% core_members && !ct %in% periphery_members)
    })
    diffuse_members <- names(which(diffuse_members))

    label <- sapply(colnames(ct_score), function(ct) {
        if (ct %in% core_members) return("core")
        if (ct %in% periphery_members) return("periphery")
        if (ct %in% diffuse_members) return("diffuse")
        return("other")
    })
    component <- NULL
    if (requireNamespace("igraph", quietly = TRUE)) {
        idx <- which(label == "core")
        core_barcodes <- names(label[idx])
        net_core <- net[which(net$from %in% core_barcodes &
                              net$to %in% core_barcodes), ]
        if (nrow(net_core)) {                      
            am <- matrix(0, nrow = length(core_barcodes), ncol = length(core_barcodes))
            colnames(am) <- core_barcodes
            rownames(am) <- core_barcodes
            for (i in seq(nrow(net_core))) {
                am[as.character(net_core$from[i]), as.character(net_core$to[i])] <- 1
                am[as.character(net_core$to[i]), as.character(net_core$from[i])] <- 1
            }
            g <- graph.adjacency(am)
            component <- components(g)$membership
            idx <- table(component)[component] == 1
            if (any(idx)) {
                component[idx] <- "Single"
            }
            component <- as.factor(component)
        } else {
            component <- label[label == "core"]
            component[component == "core"] <- "Single"
            component <- as.factor(component)
        }    
    }
    return(list(
        label = factor(label, levels = c("core", "periphery", "diffuse", "other")),
        component = component
    ))
}

