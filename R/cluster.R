#' cluster_spatial
#'
#' Normalize spatial transcriptomics Seurat object
#' @param ndata Object, read by \code{\link{read_spatial}}.
#' @param resolution Seurat cluster resolution 
#' @param dims PCs to use
#' @param sketch Number of sketch bins/cell. Ignored when < 1.
#' @param verbose Verbose Seurat output

#' @export cluster_spatial
#' @examples
#' cluster_spatial()

cluster_spatial <- function(ndata, resolution = 0.8, dims = 1:50, sketch = 0, verbose = TRUE) {
    if (sketch > max(dims) && ncol(ndata) > sketch) {
        if (sketch < 5000) flog.warn("Small sketch value, will continue anyway.")
        if (length(resolution) > 1) {
            resolution <- max(resolution)
            flog.warn("Only one resolution with sketching. Using %.2f", resolution)
        }
        def_assay <- DefaultAssay(ndata)
        flog.info("Sketching assay %s from %i to %i bins/cells.", def_assay, ncol(ndata), sketch)
        ndata <- SketchData(
          object = ndata,
          ncells = sketch,
          method = "LeverageScore",
          sketched.assay = "sketch",
          features = head(VariableFeatures(ndata), 5000)
        )
        DefaultAssay(ndata) <- "sketch"
        ndata <- FindVariableFeatures(ndata)
        ndata <- ScaleData(ndata)
        flog.info("Using resolution %f and %i PCs for clustering.", 
            resolution, length(dims))
        ndata <- RunPCA(object = ndata, assay = "sketch", npcs = min(ncol(ndata)-1, 50),
             verbose = verbose, reduction.name = "pca.sketch")
        ndata <- FindNeighbors(ndata, assay = "sketch", reduction = "pca.sketch", dims = dims, verbose = verbose)
        ndata <- FindClusters(ndata, cluster.name = "seurat_cluster.sketched", resolution = resolution, verbose = verbose)
        ndata <- RunUMAP(object = ndata, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = dims, verbose = verbose)
        flog.info("Projecting sketched data back to %s", def_assay)
        ndata <- ProjectData(
          object = ndata,
          assay = def_assay,
          full.reduction = "full.pca.sketch",
          sketched.assay = "sketch",
          sketched.reduction = "pca.sketch",
          umap.model = "umap.sketch",
          dims = dims,
          refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
        )
        DefaultAssay(ndata) <- def_assay
        Idents(ndata) <- "seurat_cluster.projected"

    } else {
        flog.info("Using resolution %f and %i PCs for clustering.", 
            resolution, length(dims))
        ndata <- RunPCA(object = ndata, npcs = min(ncol(ndata)-1, 50),
             verbose = verbose)
        # check for crappy samples with less than 30 spots/cells
        if (length(dims) > length(Seurat::Reductions(ndata, slot = .get_reduction(ndata, type = "pca")))) {
            flog.warn("Ignoring provided dims argument because ndata does not have sufficient dimension.")
            dims <- seq(1, ncol(ndata) - 1)
        }
        ndata <- RunUMAP(object = ndata, dims = dims, verbose = verbose)
        ndata <- FindNeighbors(ndata, dims = dims, verbose = verbose)
        ndata <- FindClusters(ndata, resolution = resolution, verbose = verbose)
    }
    ndata
}    

#' cluster_bayesspace
#'
#' Cluster BayesSpace data
#' @param x Object, converted by \code{\link{as_SingleCellExperiment}}
#' @param test_num_clusters Range of number of clusters to be tested (\code{qs}). 
#' @param num_clusters Optional picked number
#' @param ... Additional paramters passed to \code{BayesSpace::spatialCluster}
#' @export cluster_bayesspace
#' @examples
#' #cluster_bayesspace()

cluster_bayesspace <- function(x, test_num_clusters = seq(2, 12),
                               num_clusters = NULL, ...) {
    if (!requireNamespace("BayesSpace", quietly = TRUE)) {
        stop("Install BayesSpace package.")
    }    
    # check if already calculated
    ql <- attr(x, "q.loglik")
    if (is.null(ql) || !identical(ql$q, test_num_clusters)) {
         x <- BayesSpace::qTune(x, qs = test_num_clusters) 
         ql <- attr(x, "q.loglik")
    } else {
        flog.info("Found log likelihoods in x. Skipping qTune...")
    }
    eql <- .elbow(ql)
    attr(x, "q.auto.selected") <- eql$q_selected

    if (is.null(num_clusters)) {
        num_clusters <- eql$q_selected
    } else {
        if (length(num_clusters) != 1 || !num_clusters %in% test_num_clusters) {
            stop("Invalid num_clusters")
        }
    }
    flog.info("Running spatialCluster with q = %i...", num_clusters)
    x <- BayesSpace::spatialCluster(x, q = num_clusters, ...)

    return(x)
}
            
#' plot_clusters
#'
#' Plot clusters
#' @param obj Object, read by \code{\link{read_spatial}}.
#' @param prefix Prefix of output files
#' @param subdir Put files in a subdirectory
#' @param pdf Create PDF image files
#' @param png Create PNG image files
#' @export plot_clusters
#' @examples
#' plot_clusters()
plot_clusters <- function(obj, prefix, subdir = "snn", pdf = FALSE, png = FALSE) {
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
    filename <- .get_sub_path(prefix, subdir, "")
    filename <- .get_sub_path(prefix, file.path(subdir, "umap"),
        paste0("_", reference_technology, "_cluster.pdf"))
    if (label) { 
        flog.info("UMAP label...")
    } else {
        flog.info("UMAP idents...")
    }
    gp <- DimPlot(obj, reduction = .get_reduction(obj), label = label)
    if (requireNamespace("ggthemes", quietly = TRUE) &&
        length(levels(Idents(obj))) <= 8) {
        gp <- gp + ggthemes::scale_colour_colorblind()
    }
    if (pdf) {
        pdf(filename, width = 10, height = 5)
        print(gp)
        dev.off()
    }
    if (png) {
        png(gsub(".pdf$", ".png", filename), width = 10, height = 5,
            units = "in", res = 150)
        print(gp)
        dev.off()
    }    
    flog.info("UMAP splitted...")
    if ("call" %in% colnames(obj@meta.data)) {
        filename <- .get_sub_path(prefix, file.path(subdir, "umap"), paste0("_", reference_technology, "_cluster_call.pdf"))
        if (pdf) {
            pdf(filename, width = 10, height = 5)
            print(DimPlot(obj, reduction = .get_reduction(obj), group.by = "call", label = label))
            dev.off()
        }
        if (png) {
            png(gsub(".pdf$", ".png", filename), width = 10, height = 5,
                units = "in", res = 150)
            print(DimPlot(obj, reduction = .get_reduction(obj), group.by = "call", label = label))
            dev.off()
        }    
        filename <- .get_sub_path(prefix, file.path(subdir, "umap"), paste0("_", reference_technology, "_cluster_splitted.pdf"))
        if (pdf) {
            pdf(filename, width = 10, height = 10)
            print(DimPlot(obj, reduction = .get_reduction(obj), split.by = "new.idents", group.by = "call"))
            dev.off()
        }
        if (png) {
            png(gsub(".pdf$", ".png", filename), width = 10, height = 10, res = 150, units = "in")
            print(DimPlot(obj, reduction = .get_reduction(obj), split.by = "new.idents", group.by = "call"))
            dev.off()
        }
    }
    if ("label" %in% colnames(obj@meta.data)) {
        .plot_cluster_library(obj, field = "label", prefix = prefix,
            subdir = file.path(subdir, "umap"),
            reference_technology = reference_technology)
    } else if ("library" %in% colnames(obj@meta.data)) {
        .plot_cluster_library(obj, field = "library", prefix = prefix, 
            subdir = file.path(subdir, "umap"),
            reference_technology = reference_technology)
    }
    if ("hg19" %in% colnames(obj@meta.data) && 
        "mm10" %in% colnames(obj@meta.data)) {
        flog.info("Violinplot hg19 vs mm10...")
        filename <- .get_sub_path(prefix, file.path(subdir, "qc"),
            paste0("_", reference_technology, "_cluster_violin_call.pdf"))
        pdf(filename, width=10, height=5)
        print(VlnPlot(obj, features = c("hg19", "mm10"), sort = TRUE))
        dev.off()
        filename <- .get_sub_path(prefix, file.path(subdir, "qc"),
            paste0("_", reference_technology, "_cluster_violin_qc.pdf"))
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
    
.plot_cluster_library <- function(obj, field = "library", prefix, subdir = "snn/umap",
    reference_technology, pdf = FALSE, png = FALSE) {
    if (nrow(unique(obj[[field]])) < 2) return()
        
    flog.info("UMAP %s...", field)
    gp <- DimPlot(obj, reduction = .get_reduction(obj), split.by = field)
    if (requireNamespace("ggthemes", quietly = TRUE) &&
        length(levels(Idents(obj))) <= 8) {
        gp <- gp + ggthemes::scale_colour_colorblind()
    }
    filename <- .get_sub_path(prefix, subdir, paste0("_", reference_technology, "_cluster_", field, ".pdf"))
    if (png) {
        png(paste0(gsub(".pdf", "", filename), ".png"), width = 10, height = 5,
            res = 150, units = "in")
        print(gp)
        dev.off()
    }
    if (!pdf) return()
    
    pdf(filename, width = 10, height = 5)
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
#' \code{obj1}. Use \code{ident} to use spot identity classes.
#' @param col2 \code{meta.data} column containing cluster labels of 
#' \code{obj2}. Use \code{ident} to use spot identity classes.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param png Create, in addition to PDF, PNG files
#' @param prefix Prefix of output files

#' @export cluster_prediction_strength
#' @examples
#' cluster_prediction_strength

cluster_prediction_strength <- function(obj1, obj2, col1 = "ident", col2 = "ident",
                                        assay = "Spatial", png = FALSE,
                                        prefix) {

    .transfer_cluster_labels <- function(reference, query, refdata) {
        anchors <- FindTransferAnchors(reference = reference, query = query,
            dims = 1:30, reduction = "cca")
        predictions <- TransferData(anchorset = anchors, refdata = refdata,
            dims = 1:30, weight.reduction = "cca")           
        query <- AddMetaData(object = query, metadata = predictions)
    }    
    .calc_consistency <- function(obj1, obj2, col) {
        refdata <- FetchData(obj1, vars = col[1])[Cells(obj1),]
        query <- .transfer_cluster_labels(obj1, obj2, refdata)
        m <- .get_sample_consistency_matrix(query, "predicted.id", col)
        .plot_consistency_matrix(m, query, obj1, obj2, assay = assay, png = png, prefix = prefix)
    }
   .calc_consistency(obj1, obj2, col1)
   .calc_consistency(obj2, obj1, col2)
}
    

.get_sample_consistency_matrix <- function(query, col1, col2) {
    m <- matrix(0, ncol(query), ncol(query))
    md <- FetchData(query, c(col1, col2))
    for (i in seq(1, ncol(query)-1)) {
        for (j in seq(i + 1, ncol(query))) {
            if (md[i, 1] == md[j, 1] &&
                md[i, 2] == md[j, 2]) { 
                m[i, j] <- sum(md[, 1] == md[i, 1], na.rm=TRUE)
                m[j, i] <- sum(md[, 2] == md[i, 2], na.rm=TRUE)
            }    
        }
    }
    m
}

.plot_consistency_matrix <- function(m, query, obj1, obj2,  assay = "Spatial",
                                     png = FALSE, prefix, suffix = "") {
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

    loupe_col <- colnames(query@meta.data)[grep("^import_", colnames(query@meta.data))]
    if (length(loupe_col)) {
        print(VlnPlot(query, features = "cluster_consistency", group.by = loupe_col))
    }    
    dev.off()
    if (png) { 
        png(paste0(prefix, "_cluster_consistency_", l1, l2, suffix, ".png"),
        width = 4, height = 3.9, units = "in", res = 150)
        print(SpatialFeaturePlot(query, "cluster_consistency"))         
        dev.off()
    }    
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
#' @importFrom stats na.omit
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
    obj <- set_idents_nmf(obj, k, rank)
    d <- data.frame( "Barcode" =  .extract_barcode(obj), 
                     "NMF" = paste("Cluster", Idents(obj)))
    colnames(d)[2] <- paste0(colnames(d)[2], "_", k)
    for (i in seq_along(libs)) {
        label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
        libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
        idx <- which(obj$library == libs[i]) 
        filename <- .get_sub_path(prefix, "nmf/loupe", 
            paste0("_nmf_cluster_loupe_", k, label, libs_label, ".csv"))
        write.csv(d[idx, , drop = FALSE], file = filename, row.names = FALSE)
    }    
    filename <- .get_sub_path(prefix, "nmf/loupe", 
        paste0("_nmf_cluster_loupe_", k, "_all.csv"))
    d$Barcode <- .extract_barcode(obj, aggr = TRUE)
    write.csv(d, file = filename, row.names = FALSE)
}

#' import_loupe
#'
#' Impurt annotation from a CSV file loadable in Loupe
#' @param obj Object, clustered by \code{\link{cluster_nmf}}.
#' @param file CSV file loadable in Loupe 
#' @returns Seurat object with annotation added.
#' @export import_loupe
#' @examples
#' import_loupe
import_loupe <- function(obj, file) {
    anno <- read.csv(file, as.is = TRUE)
    label <- paste0("import_", tolower(colnames(anno)[2]))
    anno[,2] <- gsub("Cluster ", "", anno[,2])
    if (any(anno$Barcode %in% Cells(obj))) {
        obj <- AddMetaData(obj, metadata = data.frame(anno[,2], row.names=anno[,1]), col.name = label)
    }
    obj   
}
    
.extract_barcode <- function(obj, aggr = FALSE) {
    barcode <- Cells(obj)
    # get aggregated barcode
    if (aggr && any(grepl("-\\d+_\\d+$", barcode))) {
        return(gsub("-\\d+_","-", barcode))
    }
    if ("barcode" %in% colnames(obj@meta.data)) {
        # Loupe wants the barcode with number suffix
        if (any(grepl("-\\d+$", obj$barcode))) {
            barcode <- obj$barcode
        } else {
            barcode <- paste0(obj$barcode, "-1")
        }
    }
    return(barcode)
}

.extract_barcode_loupe <- function(obj, tissue_positions_list = NULL) {
    if (is.null(tissue_positions_list)) return(Cells(obj))
    barcodes1 <- strsplit(Cells(obj), "-")
    barcodes1 <- split(sapply(barcodes1, function(x) x[1]), sapply(barcodes1, function(x) x[2]))
    pos <- read.csv(tissue_positions_list, as.is = TRUE, header = FALSE)
    barcodes2 <- strsplit(pos[,1], "-")
    barcodes2 <- split(sapply(barcodes2, function(x) x[1]), sapply(barcodes2, function(x) x[2]))
    .sdiff <- function(x,y) length(intersect(x,y))/length(union(x,y))
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
                
            filename <- .get_sub_path(prefix, "snn/loupe", 
                paste0("_snn_cluster_loupe_", sid_us, label, libs_label, ".csv"))
            idx <- which(obj$library == libs[i]) 
            write.csv(d[idx, , drop = FALSE], file = filename, row.names = FALSE)
        }    
        filename <- .get_sub_path(prefix, "snn/loupe", 
            paste0("_snn_cluster_loupe_", sid_us, "_all.csv"))
        d$Barcode <- .extract_barcode(obj, aggr = TRUE)
        write.csv(d, file = filename, row.names = FALSE)
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
#' @param stop_if_unavail Dies when NMF clustering is not found.
#' If \code{NULL}, use first available.
#' @return object with new \code{Idents}
#' @importFrom stats predict
#' @export set_idents_nmf
#' @examples
#' set_idents_nmf
set_idents_nmf <- function(object, k, rank = NULL, stop_if_unavail = FALSE) {
    if (!requireNamespace("NMF", quietly = TRUE)) {
        stop("This function requires the NMF library.")
    }

    if (is.null(rank)) {
        available_misc <- names(object@misc)
        available_misc <- available_misc[grep("^nmf", available_misc)]
        available_misc <- available_misc[!grepl("_random", available_misc)]
        if (!length(available_misc)) {
            if (stop_if_unavail) stop("NMF clustering not found.")
            return(object)
        }    
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
    object[["orig.ident.before.nmf"]] <- Idents(object = object)
    Idents(object) <- NMF::predict(nmf_obj_f)
    return(object)
}    

.elbow <- function(data) {
  # This code is taken from   
  # https://github.com/ahasverus/elbow/blob/master/R/elbow.R


  if (!is.list(data)) {
    stop("`data` must be a two-columns data frame.")
  }

  if (!is.numeric(data[ , 1]) || !is.numeric(data[ , 2])) {
    stop("Non-numeric data detected.")
  }

  if (sum(is.na(data[ , 1])) + sum(is.na(data[ , 2]))) {
    stop("Missing values detected.")
  }

  ## Data transformation ----

  data <- data[ , 1:2]
  data <- data[order(data[ , 1]), ]

  ## Get constant increase/decrease in y ----

  constant <- data[c(1, nrow(data)), ]
  colnames(constant) <- c("x", "y")

  mod <- stats::lm(y ~ x, data = constant)

  data[ , "constant"] <- round(mod$coef[[1]] + mod$coef[[2]] * data[ , 1], 3)

  ## Detect inflection point ----
  pos <- round(nrow(data) / 2)

  if (data[pos, "constant"] < data[pos, 2]) { # Concave Down
    ymin <- min(data[ , 2])
    data[ , "benefits"] <- ymin + round(data[ , 2] - data[ , "constant"], 3)
    maxi <- data[which.max(data[ , "benefits"]), ]
  } else { # Concave Up
    ymax <- max(data[ , 2])
    data[ , "benefits"] <- ymax - round(data[ , "constant"] - data[ , 2], 3)
    maxi <- data[which.min(data[ , "benefits"]), ]
  }

  xxx <- list()
  xxx[[1]] <- maxi[1, 1]
  xxx[[2]] <- data
  names(xxx) <- c(paste(colnames(data)[1], "selected", sep = "_"), "data")
  return(xxx)
}

.merge_bayesspace <- function(x, y) {
# there must be a cleaner way, let me know in case there is!
    if (!identical(dim(x), dim(y))) return(NULL)
    features_x <- names(which(!is.na(rowData(x)$enhanceFeatures.rmse)))
    features_y <- names(which(!is.na(rowData(y)$enhanceFeatures.rmse)))
    if (length(intersect(features_x, features_y))) {
        flog.warn("Duplicated features in .merge_bayesspace")
        features_y <- features_y[!features_y %in% features_x]
    }
    if (!length(features_x)) return(y)
    if (!length(features_y)) return(x)
    m <- x
    SingleCellExperiment::logcounts(m)[features_y,] <- SingleCellExperiment::logcounts(y)[features_y,]
    rowData(m)[features_y,] <- rowData(y)[features_y,]
    return(m)
}    

#' cluster_sc3
#'
#' Performes SC3 clustering 
#' @param obj Object, clustered by \code{\link{cluster_spatial}}.
#' @param rank Number of clusters 
#' @param variable_features If \code{TRUE}, only use variable features
#' @param max_features Reduce runtime by only using the top 
#' \code{max_features} features
#' @param force Recalculate, even when serialized objects are available
#' @param serialize Serialize output objects
#' @param prefix Prefix of output files
#' @param ... Additional parameters passed to the \code{sc3} function.
#' @export cluster_sc3
#' @examples
#' cluster_sc3
cluster_sc3 <- function(obj, rank, variable_features = FALSE, 
    max_features = NULL, force, serialize = TRUE, prefix, ...) {
    if (!requireNamespace("SC3", quietly = TRUE)) {
        stop("This function requires the SC3 library.")
    }
    if (!requireNamespace("scater", quietly = TRUE)) {
        stop("This function requires the scater library.")
    }
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("This function requires the SingleCellExperiment library.")
    }
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("This function requires the SummarizedExperiment library.")
    }
    rank <- sort(unique(rank))

    if (length(rank) == 1) {
        key <- paste0("sc3_k_", rank)
    } else {
        key <- paste0("sc3_k_", min(rank), "_to_", max(rank))
    }
    if (!is.null(max_features)) {
        obj <- FindVariableFeatures(obj, nfeatures = max_features,
            selection.method = "disp")
    }    
    if (variable_features) {
        obj <- obj[VariableFeatures(obj),]
    }
    flog.info("Performing SC3 clustering on %i samples and %i features.",
        ncol(obj), nrow(obj))
    
    if (!is.null(obj@misc[[key]])) {
        flog.warn("Object %s in misc slot already exists. Skipping...")
        return(obj)
    }
    sce <- as.SingleCellExperiment(obj)
    SummarizedExperiment::rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sce[!duplicated(SummarizedExperiment::rowData(sce)$feature_symbol), ]
    sce <- scater::runPCA(sce)
    SingleCellExperiment::counts(sce) <- as.matrix(SingleCellExperiment::counts(sce))
    SingleCellExperiment::logcounts(sce) <- as.matrix(SingleCellExperiment::logcounts(sce))
    sce <- SC3::sc3(sce, ks = rank, ...) 
    if (serialize) {
        filename <- .get_serialize_path(prefix, paste0("_", key, ".rds"))
        if (!file.exists(filename) || force) {
            flog.info("Writing R data structure to %s...", filename)
            saveRDS(sce, filename)
        }    
    }
    obj@misc[[key]] <- list(
        rank = rank,
        sc3_rowData = SummarizedExperiment::rowData(sce),
        sc3_colData = SummarizedExperiment::colData(sce)
    )    
    obj
}


#' set_idents_sc3
#'
#' Sets \code{Idents} to SC3 clustering
#' @param obj Object, clustered by \code{\link{cluster_sc3}}.
#' @param k Features of rank to be written (must be a single k, not a range)
#' @param rank Number of clusters (the one used in \code{\link{cluster_sc3}})
#' @param stop_if_unavail Dies when SC3 clustering is not found.
#' If \code{NULL}, use first available.
#' @return obj with new \code{Idents}
#' @export set_idents_sc3
#' @examples
#' set_idents_sc3
set_idents_sc3 <- function(obj, k, rank = NULL, stop_if_unavail = FALSE) {
    if (!requireNamespace("SC3", quietly = TRUE)) {
        stop("This function requires the SC3 library.")
    }

    if (is.null(rank)) {
        available_misc <- names(obj@misc)
        available_misc <- available_misc[grep("^sc3", available_misc)]
        if (!length(available_misc)) {
            if (stop_if_unavail) stop("SC3 clustering not found.")
            return(obj)
        }    
        rank <- obj@misc[[available_misc]]$rank
    } else {
        if (length(rank) > 1) {
            available_misc <- paste0("sc3_k_", min(rank), "_to_", max(rank))
        } else {
            available_misc <- paste0("sc3_k_", rank)
        }       
    }    
    if (k < min(rank)) {
        flog.warn("requested k not available.")
        k <- min(rank)
    } else if (k > max(rank)) {
        flog.warn("requested k not available.")
        k <- max(rank)
    }
    obj[["orig.ident.before.sc3"]] <- Idents(object = obj)

    Idents(obj) <- obj@misc[[available_misc]]$sc3_colData[colnames(obj),paste0("sc3_", k, "_clusters")]
    return(obj)
}    

.get_reduction <- function(obj, type = "umap") {
    avail <- grep(type, Seurat::Reductions(obj), value = TRUE)
    if (type %in% avail) return(type)
    if (paste0(type, ".sketch") %in% avail) return(paste0(type, ".sketch"))
    return(avail[1])
}
