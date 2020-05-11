#' cluster_spatial
#'
#' Normalize spatial transcriptomics Seurat object
#' @param ndata Object, read by \code{\link{read_spatial}}.
#' @param resolution Seurat cluster resolution 
#' @param dims PCs to use
#' @param verbose Verbose Seurat output

#' @export cluster_spatial
#' @examples
#' cluster_spatial()

cluster_spatial <- function(ndata, resolution = 0.8, dims = 1:30, verbose = TRUE) {
    flog.info("Using resolution %f and %i PCs for clustering.", 
        resolution, length(dims))
    ndata <- RunPCA(object = ndata, npcs = min(ncol(ndata)-1, 50),
         verbose = verbose)
    ndata <- RunUMAP(object = ndata, dims = dims, verbose = verbose)
    ndata <- FindNeighbors(ndata, dims = dims, verbose = verbose)
    ndata <- FindClusters(ndata, resolution = resolution, verbose = verbose)
    ndata
}    

#' plot_clusters
#'
#' Plot clusters
#' @param obj Object, read by \code{\link{read_spatial}}.
#' @param prefix Prefix of output files

#' @export plot_clusters
#' @examples
#' plot_clusters()
plot_clusters <- function(obj, prefix) {
    reference_technology <- "spatial"
    if ("technology" %in% colnames(obj@meta.data)) {
        reference_technology <- obj$technology[obj$reference][1]
    }
    if (!"new.idents" %in% colnames(obj@meta.data)) {
        obj$new.idents <- Idents(obj)
    }
    obj$new.idents <- factor(as.character(obj$new.idents), 
        levels = .order_clusters(obj$new.idents))
    label <- "predicted.id" %in% colnames(obj@meta.data)
    # make sure than snn directory is created
    filename <- .get_sub_path(prefix, "snn", "")
    filename <- .get_sub_path(prefix, "snn/umap", paste0("_", reference_technology, "_cluster.pdf"))
    pdf(filename, width = 10, height = 5)
    if (label) { 
        flog.info("UMAP label...")
    } else {
        flog.info("UMAP idents...")
    }
    gp <- DimPlot(obj, reduction = "umap", label = label)
    if (requireNamespace("ggthemes", quietly = TRUE)) {
        gp <- gp + ggthemes::scale_colour_colorblind()
    }
    print(gp)
    dev.off()
    flog.info("UMAP splitted...")
    if ("call" %in% colnames(obj@meta.data)) {
        filename <- .get_sub_path(prefix, "snn/umap", paste0("_", reference_technology, "_cluster_call.pdf"))
        pdf(filename, width = 10, height = 5)
        print(DimPlot(obj, reduction = "umap", group.by = "call", label = label))
        dev.off()
        filename <- .get_sub_path(prefix, "snn/umap", paste0("_", reference_technology, "_cluster_splitted.pdf"))
        pdf(filename, width = 10, height = 10)
        print(DimPlot(obj, split.by = "new.idents", group.by = "call"))
        dev.off()
    }
    if ("label" %in% colnames(obj@meta.data)) {
        .plot_cluster_library(obj, field = "label", prefix = prefix, 
            reference_technology = reference_technology)
    } else if ("library" %in% colnames(obj@meta.data)) {
        .plot_cluster_library(obj, field = "library", prefix = prefix, 
            reference_technology = reference_technology)
    }
    if ("hg19" %in% colnames(obj@meta.data) && 
        "mm10" %in% colnames(obj@meta.data)) {
        flog.info("Violinplot hg19 vs mm10...")
        filename <- .get_sub_path(prefix, "snn/qc", paste0("_", reference_technology, "_cluster_violin_call.pdf"))
        pdf(filename, width=10, height=5)
        print(VlnPlot(obj, features = c("hg19", "mm10"), sort = TRUE))
        dev.off()
        filename <- .get_sub_path(prefix, "snn/qc", paste0("_", reference_technology, "_cluster_violin_qc.pdf"))
        pdf(filename, width = 10, height = 5)
        print(VlnPlot(obj, features = c("percent.mito", "percent.ribo"), sort = TRUE))
        dev.off()
    }
}

.find_technical_replicates <- function(labels) {
    labels <- gsub("_\\d+$","", labels)
    tbl <- table(labels)
    lapply(names(tbl[tbl>1]), function(x) which(labels == x)) 
}
    
.plot_cluster_library <- function(obj, field = "library", prefix, subdir = "snn/umap", reference_technology) {
    flog.info("UMAP %s...", field)
    filename <- .get_sub_path(prefix, subdir, paste0("_", reference_technology, "_cluster_", field, ".pdf"))

    pdf(filename, width = 10, height = 5)
    gp <- DimPlot(obj, reduction = "umap", split.by = field)
    if (requireNamespace("ggthemes", quietly = TRUE) &&
        length(levels(Idents(obj))) <= 8) {
        gp <- gp + ggthemes::scale_colour_colorblind()
    }
    print(gp)
    
    df <- melt(table(obj@meta.data[[field]], Idents(obj)))
    print(ggplot(df, aes_string("Var1", "value")) + geom_col() + 
        facet_wrap(~Var2) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ylab("Counts") + xlab("")
    )
    dfx <- data.frame(Label=obj@meta.data[[field]], Cluster=Idents(obj))
    num_labels <- length(unique(dfx$Label))
    if (num_labels > 1) {
        gp <- ggplot(dfx, aes_string("Cluster", fill = "Label")) +
            geom_bar(position="fill") + 
            scale_y_continuous(labels = scales::percent) + 
            ylab("Label")
        if (requireNamespace("ggthemes", quietly = TRUE) &&
            nrow(unique(obj[[field]])) <= 8) {
            gp <- gp + ggthemes::scale_fill_colorblind()
        }
        print(gp)
    }
    # collapse technical replicates
    dfx$Label <- gsub("_\\d+$","", dfx$Label)
    num_labels_2 <- length(unique(dfx$Label))
    if (num_labels_2 < num_labels && num_labels_2 > 1) {
        gp <- ggplot(dfx, aes_string("Cluster", fill = "Label")) +
            geom_bar(position="fill") + 
            scale_y_continuous(labels = scales::percent) + 
            ylab("Label")
        if (requireNamespace("ggthemes", quietly = TRUE) && 
            length(unique(dfx$Label)) <= 8) {
            gp <- gp + ggthemes::scale_fill_colorblind()
        }
        print(gp)
    }
    dev.off()
}

#' cluster_prediction_strength
#'
#' Calculates the cluster consistency across two technical 
#' replicates. Shows the fraction of spot pairs clustered together
#' for each spot in both samples.
#' @param obj1 Object, clustered by \code{\link{cluster_spatial}}.
#' @param obj2 Object, clustered by \code{\link{cluster_spatial}}.
#' @param col1 \code{meta.data} column containing cluster labels of 
#' \code{obj1}. Default use first available.
#' @param col2 \code{meta.data} column containing cluster labels of 
#' \code{obj2}. Default use first available.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param prefix Prefix of output files

#' @export cluster_prediction_strength
#' @examples
#' cluster_prediction_strength

cluster_prediction_strength <- function(obj1, obj2, col1 = NULL, col2 = NULL,
                                        assay = "Spatial",
                                        prefix) {

    .find_cluster_col <- function(obj, label) {
        col <- grep("snn_res", colnames(obj@meta.data))
        if (!length(col)) stop("cannot find cluster column for %s.", label)
        col <- colnames(obj1@meta.data)[col]
    }    
    if (is.null(col1)) col1 <- .find_cluster_col(obj1, "obj1")
    if (is.null(col2)) col2 <- .find_cluster_col(obj2, "obj2")
    .transfer_cluster_labels <- function(reference, query, refdata) {
        anchors <- FindTransferAnchors(reference = reference, query = query,
            dims = 1:30, reduction = "cca")
        predictions <- TransferData(anchorset = anchors, refdata = refdata,
            dims = 1:30, weight.reduction = "cca")           
        query <- AddMetaData(object = query, metadata = predictions)
    }    
    .calc_consistency <- function(obj1, obj2, col) {
        query <- .transfer_cluster_labels(obj1, obj2, obj1@meta.data[[col]])
        m <- .get_sample_consistency_matrix(query, "predicted.id", col)
        .plot_consistency_matrix(m, query, obj1, obj2, assay = assay, prefix = prefix)
    }
   .calc_consistency(obj1, obj2, col1)
   .calc_consistency(obj2, obj1, col2)
}
    

.get_sample_consistency_matrix <- function(query, col1, col2) {
    m <- matrix(0, ncol(query), ncol(query))
    for (i in seq(1, ncol(query)-1)) {
        for (j in seq(i + 1, ncol(query))) {
            if (query@meta.data[i, col1] == query@meta.data[j, col1] &&
                query@meta.data[i, col2] == query@meta.data[j, col2]) { 
                m[i, j] <- sum(query@meta.data[, col1] == query@meta.data[i, col1], na.rm=TRUE)
                m[j, i] <- sum(query@meta.data[, col2] == query@meta.data[i, col2], na.rm=TRUE)
            }    
        }
    }
    m
}

.plot_consistency_matrix <- function(m, query, obj1, obj2,  assay = "Spatial", prefix, suffix = "") {
    query$cluster_consistency <- apply(m, 1, function(x) sum(x > 0)) / 
        pmax(1, pmax(apply(m, 1, function(x) max(x)), 
                     apply(m, 2, function(x) max(x))))

    l1 <- if (!is.null(obj1@meta.data$library)) obj1$library[1] else "obj1"
    l2 <- if (!is.null(obj2@meta.data$library)) obj2$library[1] else "obj2"
    l2 <- if (l2 == l1)  "" else paste0("_", l2)
    pdf(paste0(prefix, "_cluster_consistency_", l1, l2, suffix, ".pdf"),
        width = 4, height = 3.9)
    print(SpatialFeaturePlot(query, "cluster_consistency"))         
    print(FeatureScatter(object = query, feature1 = "cluster_consistency",
                                  feature2 = paste0("nFeature_", assay)))
    print(FeatureScatter(object = query, feature1 = "cluster_consistency",
                                         feature2 = "percent.mito"))
    print(FeatureScatter(object = query, feature1 = "cluster_consistency",
                                         feature2 = "percent.ribo"))
    dev.off()
    query
}


#' cluster_nmf
#'
#' Performes NMF clustering 
#' @param obj Object, clustered by \code{\link{cluster_spatial}}.
#' @param rank Number of clusters 
#' @param randomize Randomize data, useful for diagnostics and picking rank
#' @param variable_features If \code{TRUE}, only use variable features
#' @param max_features Reduce runtime by only using the top 
#' \code{max_features} features
#' @param ... Additional parameters passed to the \code{nmf} function.
#' @export cluster_nmf
#' @examples
#' cluster_nmf
cluster_nmf <- function(obj, rank, randomize = FALSE, variable_features = TRUE, 
    max_features = NULL, ...) {
    if (!requireNamespace("NMF", quietly = TRUE)) {
        stop("This function requires the NMF library.")
    }
    if (length(rank) == 1) {
        key <- paste0("nmf_k_", rank)
    } else {
        key <- paste0("nmf_k_", min(rank), "_to_", max(rank))
    }
    if (!is.null(max_features)) {
        obj <- FindVariableFeatures(obj, nfeatures = max_features,
            selection.method = "disp")
    }    
    m <- GetAssayData(obj)
    if (variable_features) {
        m <- m[VariableFeatures(obj),]
    }
    flog.info("Performing clustering on %i samples and %i features.",
        ncol(m), nrow(m))
    
    m <- as.matrix(m) - min(m)
    if (randomize) {
        key <- paste0(key, "_random")
        m <- NMF::randomize(m)
    }
    if (!is.null(obj@misc[[key]])) {
        flog.warn("Object %s in misc slot already exists. Skipping...")
        return(obj)
    }        
    nmf_obj <- NMF::nmf(m, rank, ...) 
    obj@misc[[key]] <- nmf_obj
    # fill meta.data slot with cluster scores
    if (!randomize) {
        if (is(nmf_obj, "NMFfit")) {
            coef <- NMF::scoef(nmf_obj)
            rownames(coef) <- paste0(key, "_", seq(1,rank))
            obj@meta.data <- cbind(obj@meta.data, t(coef))
        } else {
            for (k in rank) {
                coef <- NMF::scoef(nmf_obj$fit[[as.character(k)]])
                rownames(coef) <- paste0("nmf_k_", k, "_", seq(1, k))
                obj@meta.data <- cbind(obj@meta.data, t(coef))
            }
        }
    }
    obj
}

#' write_nmf_features
#'
#' Output metagenes to a CSV file
#' @param obj Object, clustered by \code{\link{cluster_nmf}}.
#' @param rank Number of clusters (the one used in \code{\link{cluster_nmf}})
#' @param k Features of rank to be written (must be a single k, not a range)
#' @param min_features Minimum number of features for each cluster to be reported
#' in case \code{method} returns less.
#' @param method Parameter of \code{NMF::extractFeatures}
#' @param prefix Prefix of output files
#' @returns \code{data.frame} with features.
#' @export write_nmf_features
#' @examples
#' write_nmf_features
write_nmf_features <- function(obj, rank, k, min_features = 20, method = "kim", prefix) {
    nmf_obj <- .extract_nmf_obj(obj, rank)
    nmf_obj_f <- if (is(nmf_obj, "NMFfit")) nmf_obj else nmf_obj$fit[[as.character(k)]]

    features <- .extract_nmf_features(nmf_obj_f, method = method, min_features = min_features )
    idx <- sapply(features, function(x) !is.null(nrow(x)))

    features_all <- do.call(rbind, 
        lapply(seq_along(features)[idx], 
            function(i) data.frame(Gene = rownames(features[[i]]), 
                               features[[i]], K = i)))
    tmp <- .get_sub_path(prefix, file.path("nmf", "advanced"), "") # make sure that advanced directory exists
    filename <- .get_sub_path(prefix, "nmf", paste0("_nmf_cluster_", k, ".csv"))
    write.csv(features_all, file = filename, row.names = FALSE)
    filename <- .get_sub_path(prefix, file.path("nmf", "advanced", k), paste0("_nmf_cluster_", k, "_all_basis.csv"))
    write.csv(NMF::basis(nmf_obj_f), file = filename)
    filename <- .get_sub_path(prefix, file.path("nmf", "advanced", k), paste0("_nmf_cluster_", k, "_all_coef.csv"))
    write.csv(t(NMF::coef(nmf_obj_f)), file = filename)
    return(features_all)
}

#' export_nmf_loupe
#'
#' Output NMF clustering to a CSV file loadable in Loupe
#' @param obj Object, clustered by \code{\link{cluster_nmf}}.
#' @param rank Number of clusters (the one used in \code{\link{cluster_nmf}})
#' @param k Features of rank to be written (must be a single k, not a range)
#' @param libs Library ids, must be stored in \code{obj$library}
#' @param labels Optional \code{character(n)} with sample labels
#' @param prefix Prefix of output files
#' @export export_nmf_loupe
#' @examples
#' export_nmf_loupe
export_nmf_loupe <- function(obj, rank, k, libs, labels = NULL, prefix) {
    nmf_obj <- .extract_nmf_obj(obj, rank)
    nmf_obj_f <- if (is(nmf_obj, "NMFfit")) nmf_obj else nmf_obj$fit[[as.character(k)]]
    m <- t(NMF::scoef(nmf_obj$fit[[as.character(k)]]))
    d <- data.frame( "Barcode" =  .extract_barcode(obj), 
                     "NMF" = paste("Cluster", apply(m, 1, which.max)))
    colnames(d)[2] <- paste0(colnames(d)[2], "_", k)
    for (i in seq_along(libs)) {
        label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
        libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
        idx <- which(obj$library == libs[i]) 
        filename <- .get_sub_path(prefix, "nmf/loupe", 
            paste0("_nmf_cluster_loupe_", k, label, libs_label, ".csv"))
        write.csv(d[idx, , drop = FALSE], file = filename, row.names = FALSE)
    }    
}

.extract_barcode <- function(obj) {
    barcode <- colnames(obj)
    if ("barcode" %in% colnames(obj@meta.data)) {
        # Loupe wants the barcode with number suffix
        if (any(grepl("-1$", obj$barcode))) {
            barcode <- obj$barcode
        } else {
            barcode <- paste0(obj$barcode, "-1")
        }
    }
    return(barcode)
}
    
#' export_snn_loupe
#'
#' Output SNN clustering to a CSV file loadable in Loupe
#' @param obj Object, clustered by \code{\link{cluster_nmf}}.
#' @param libs Library ids, must be stored in \code{obj$library}
#' @param labels Optional \code{character(n)} with sample labels
#' @param prefix Prefix of output files
#' @export export_snn_loupe
#' @examples
#' export_snn_loupe
export_snn_loupe <- function(obj, libs, labels = NULL, prefix) {
    barcode <- .extract_barcode(obj)

    sids <- grep("snn_res", colnames(obj@meta.data))
    for (i in sids) {
        sid_us <- gsub("\\.","_", colnames(obj@meta.data)[i])
        sid_us <- gsub("snn_res_", "", sid_us)
        id <- as.numeric(obj@meta.data[, i]) 
        if (min(id, na.rm = TRUE) < 1) id <- id + 1
        d <- data.frame( "Barcode" =  barcode, 
                         "SNN" = paste("Cluster", id))
        colnames(d)[2] <- paste0(colnames(d)[2], "_", sid_us)

        for (i in seq_along(libs)) {
            label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
            libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
                
            filename <- sttkit:::.get_sub_path(prefix, "snn/loupe", 
                paste0("_snn_cluster_loupe_", sid_us, label, libs_label, ".csv"))
            idx <- which(obj$library == libs[i]) 
            write.csv(d[idx, , drop = FALSE], file = filename, row.names = FALSE)
        }    
    }
}
    
.extract_nmf_features <- function(nmf_obj_f, method = "kim", min_features = 20) {
    features_method <- NMF::extractFeatures(nmf_obj_f, nodups = FALSE, method = method)
    features_min <- NMF::extractFeatures(nmf_obj_f, nodups = FALSE, method = min_features)
    features <- lapply(seq_along(features_min), 
        function(i) NMF::basis(nmf_obj_f)[na.omit(unique(c(features_min[[i]], features_method[[i]]))),])
    return(features)
}

.extract_nmf_r2 <- function(obj, rank, k) {
    nmf_obj <- .extract_nmf_obj(obj, rank)
    nmf_obj_f <- if (is(nmf_obj, "NMFfit")) nmf_obj else nmf_obj$fit[[as.character(k)]]

    nmf_m <- NMF::fitted(nmf_obj_f)
    nmf_d <- as.matrix(GetAssayData(obj))
    idx <- intersect(rownames(nmf_d), rownames(nmf_m))
   
    r2 <- sapply(seq(ncol(nmf_d)), function(i) summary(lm(nmf_m[idx,i]~nmf_d[idx,i]))$r.squared)
    obj[[paste0("nmf_k_r2_", k)]] <- r2
    obj
}
.extract_nmf_rss <- function(obj, rank, k) {
    nmf_obj <- .extract_nmf_obj(obj, rank)
    nmf_obj_f <- if (is(nmf_obj, "NMFfit")) nmf_obj else nmf_obj$fit[[as.character(k)]]
    nmf_m <- NMF::fitted(nmf_obj_f)
    nmf_d <- as.matrix(GetAssayData(obj))
    idx <- intersect(rownames(nmf_d), rownames(nmf_m))
    nmf_x <- scale(nmf_d[idx,], scale = FALSE) - 
             scale(nmf_m[idx,], scale  = FALSE)

    rss <- apply(nmf_x, 2, function(x) sum(x^2))
    obj[[paste0("nmf_k_rss_", k)]] <- rss
    obj
}

init_nmf_seed <- function(obj, b) {
    m <- GetAssayData(obj)
    m <- as.matrix(m) - min(m)

}

#' set_idents_nmf 
#'
#' Sets \code{Idents} to a discrete NMF clustering
#' @param object Object, clustered by \code{\link{cluster_nmf}}.
#' @param k Features of rank to be written (must be a single k, not a range)
#' @param rank Number of clusters (the one used in \code{\link{cluster_nmf}})
#' If \code{NULL}, use first available.
#' @return object with new \code{Idents}
#' @export set_idents_nmf
#' @examples
#' set_idents_nmf
set_idents_nmf <- function(object, k, rank = NULL) {
    if (is.null(rank)) {
        available_misc <- names(object@misc)
        available_misc <- available_misc[grep("^nmf", available_misc)]
        available_misc <- available_misc[!grepl("_random", available_misc)]
        if (!length(available_misc)) return(object)
        rank <- as.numeric(names(object@misc[[available_misc]]$fit))
    }    
    nmf_obj <- .extract_nmf_obj(object, rank)
    if (k < min(rank)) {
        flog.warn("requested k not available.")
        k <- min(rank)
    } else if (k > max(rank)) {
        flog.warn("requested k not available.")
        k <- max(rank)
    }
    nmf_obj_f <- if (is(nmf_obj, "NMFfit")) nmf_obj else nmf_obj$fit[[as.character(k)]]
    Idents(object) <- predict(nmf_obj_f)
    return(object)
}    
