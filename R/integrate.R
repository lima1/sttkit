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
                            reference_technology = "single_cell", force, 
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
    if (min(num_spots) < min_spots) {
        poor_libs <- sapply(references[num_spots < min_spots], function(x) x$library[1])
        flog.warn("Samples with low number of cells/spots: %s", 
            paste(poor_libs, collapse=", "))
        if (length(poor_libs)>1) references <- .merge_poor_libs(references, min_spots)
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
        max_counts <- apply(GetAssayData(x, "counts"),1, max) 
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
#' @param verbose Verbose Seurat output
#' @export cluster_integrated
#' @examples
#' cluster_integrated()
cluster_integrated <- function(integrated, regressout, scale = NULL,
                               force, plot_umap = TRUE,
                               serialize = TRUE, prefix, verbose = FALSE) {
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
    integrated <- RunUMAP(object = integrated, reduction = "pca", 
        dims = 1:30, verbose = verbose)
    if (plot_umap) {
        flog.info("Plotting UMAP...")
        pdf(.get_advanced_path(prefix, "_integration_overview.pdf"), width = 10, height = 5)
        if ("call" %in% colnames(integrated@meta.data)) {
            print(DimPlot(integrated, reduction = "umap", split.by = "technology",
                          group.by = "call"))
        } else {
            print(DimPlot(integrated, reduction = "umap", split.by = "technology",
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
                        resolution_label, "predictedid.pdf")), width=10, height=10)
                   print(DimPlot(obj_ref, reduction = "umap", split.by = "predicted.id"))
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
    Reduce(merge, lapply(seq_along(obj), function(i) {
        RenameCells(obj[[i]], add.cell.id = i)
    }))
}

.remove_undetected_features <- function(references, min_detected) {
    # find number of undetected references for all features undetected
    # in at least one reference
    tbl <- table(unlist(lapply(references, function(x) { 
        y <- Matrix::rowSums(GetAssayData(x, "counts"))
        names(y[which(y <= 0)])
        })))
    tbl_detected <- length(references) - tbl
    undetected_features <- names(tbl_detected[tbl_detected < min_detected])
    if (!length(undetected_features)) return(references)
    flog.warn("Removing features detected in fewer than %i references: %s", 
        min_detected, paste(undetected_features, collapse = ","))    
    lapply(references, function(x) x[!rownames(x) %in% undetected_features, ])
}
    
