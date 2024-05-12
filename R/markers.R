
#' find_markers
#'
#' Integrate SpatialTranscriptomics data with (matched) scRNA-seq data
#' @param obj Seurat object containing reference data obtained by
#' \code{\link{cluster_reference}}
#' @param references List of Seurat scRNA-seq reference datasets 
#' @param resolution Cluster resolution
#' @param max_markers Use only top markers
#' @param merged Flag indicating whether clusters were already merged
#' @param write Write output text files and PDFs
#' @param force Recalculate, even when serialized objects are available
#' @param serialize Serialize output objects
#' @param prefix Prefix of output files
#' @param label Label of \code{obj}
#' @param verbose Verbose Seurat output
#' @param ... Additional parameters passed to \code{Seurat::FindMarkers}
#' @export find_markers
#' @examples
#' find_markers()
#' @import ggplot2
#' @importFrom utils head
find_markers <- function(obj, references = NULL, resolution, max_markers = NULL, merged = FALSE, write = TRUE, force = FALSE,
    serialize = TRUE, prefix, label = NULL, verbose = TRUE, ...) {
    merged <- ifelse(merged, "merged_", "")
    label <- if (is.null(label)) "" else paste0(label, "_")

    if (!"reference" %in% colnames(obj@meta.data)) {
        obj$reference <- TRUE
    }    
    if (!"technology" %in% colnames(obj@meta.data)) {
        obj$technology <- "single_cell"
    }    
    reference_technology <- obj$technology[obj$reference][1]
    if (serialize || !force) { 
        filename <- .get_serialize_path(prefix, paste0("_", reference_technology, 
            "_cluster_signature_", merged, resolution, ".rds"))
    }
    if (!force && file.exists(filename)) {
        flog.warn("%s exists. Skipping signature search. Use --force to overwrite.", filename)
        return(readRDS(filename))
    }
    flog.info("Finding all markers of integrated and clustered scRNA data...")
    all_markers_unfiltered <- .find_all_markers(obj, prefix = prefix,
        suffix = paste0("_", reference_technology, "_", merged, resolution, "_markers.rds"),
        only.pos = TRUE, verbose = verbose, ...)
    all_markers <- all_markers_unfiltered

    n_before <- nrow(all_markers)
    all_markers <- all_markers[!grepl(regex_mito(), all_markers$gene),]
    if (nrow(all_markers) < n_before) {
        flog.info("Removing %i mitochondria markers.", n_before - nrow(all_markers))
    }
         
    n_before <- nrow(all_markers)
    all_markers <- all_markers[!grepl(regex_ribo(), all_markers$gene),]
    if (nrow(all_markers) < n_before) {
        flog.info("Removing %i ribosomal markers.", n_before - nrow(all_markers))
    }
         
    n_before <- nrow(all_markers)
    if (is.null(references)) {
        all_markers <- .filter_markers_maxref(list(obj), all_markers)
    } else {    
        all_markers <- .filter_markers_maxref(references, all_markers)
    }    
    if (nrow(all_markers) < n_before) {
        flog.info("Removing %i markers with low counts in references.", n_before - nrow(all_markers))
    } 
    
    obj_sc_average <- AverageExpression(obj, verbose = verbose)
    all_markers <- all_markers[rownames(all_markers) 
        %in% rownames(obj_sc_average[[1]]),]

    all_markers <- .filter_markers_max(all_markers, max_markers)

    obj_sc_sig <- as.matrix(obj_sc_average[[1]][all_markers$gene,])
    rownames(obj_sc_sig) <-gsub("obj-", "", rownames(obj_sc_sig))
    attr(obj_sc_sig, "markers") <- all_markers
    if (serialize) {
        flog.info("Writing R data structure to %s...", filename)
        saveRDS(obj_sc_sig, filename)
    }    
    if (write) {
        filename <- .get_sub_path(prefix, "markers",
            paste0("_single_cell_cluster_signature_", merged, resolution, ".csv"))
        write.csv(all_markers, filename, row.names = FALSE)
        filename <- .get_sub_path(prefix, "markers",
            paste0("_single_cell_cluster_signature_unfiltered_", merged, resolution, ".csv"))
        write.csv(all_markers_unfiltered, filename, row.names = FALSE)
        filename <- .get_sub_path(prefix, "markers",
            paste0("_single_cell_cluster_signature_", merged, resolution, ".pdf"))
        flog.info("Plotting signature heatmap to %s...", filename)
        if (ncol(obj) > 5000) {
            obj <- subset(obj, downsample = 300)
        }    
        pdf(filename, width = 10, height = 5)
        p <- DoHeatmap(obj, features = rownames(obj_sc_sig), label = FALSE) + 
            theme(axis.text.y = element_blank())
        print(p)
        dev.off()
    }
    obj_sc_sig
}
.filter_markers_maxref <- function(references, all_markers, min_max = 20) {
     max_counts <- apply(do.call(cbind, 
        lapply(references, function(x) 
            apply(GetAssayData(x[all_markers$gene,], slot = "counts"),1,max))), 
     1, max)
     all_markers$max_counts <- max_counts[all_markers$gene]
     all_markers[which(all_markers$max_counts >= min_max),]
}

.filter_markers_cellmix <- function(obj_sc, all_markers) {
    if (requireNamespace("CellMix", quietly = TRUE)) {
        flog.info("Extracting markers...")
        ems <- lapply(c(0.01, 0.001, 0.0001), function(threshold) 
            CellMix::extractMarkers(as.matrix(GetAssayData(obj_sc)[unique(all_markers$gene),]), 
            Idents(obj_sc), threshold = threshold, lbase = exp(1)))
        
        filtered <- lapply(ems, function(em) lapply(names(CellMix::geneIds(em)), 
            function(i) paste(i, CellMix::geneIds(em)[[i]])))
        labels <- do.call(cbind, lapply(filtered, function(x) paste(all_markers$cluster, 
            all_markers$gene) %in% unlist(x)))
        colnames(labels) <- paste0("em", seq_along(filtered))
        
        all_markers <- cbind(all_markers, labels)
        all_markers <- all_markers[all_markers$em1,]     
    } else {
        flog.warn("Install CellMix package!")
    }        
    return(all_markers)
}

.filter_markers_max <- function(all_markers, max_markers) {
    if (is.null(max_markers)) return(all_markers)
    markers <- split(all_markers, all_markers$cluster)
    do.call(rbind, lapply(markers, function(x) 
        head(x[order(!x$em3, !x$em2, x$p_val),], max_markers)))
}        

.find_all_markers <- function(obj, prefix, suffix, ...) {
    filename_markers <- .get_serialize_path(opt$outprefix, suffix)
    if (!opt$force && file.exists(filename_markers)) {
        flog.warn("%s exists. Skipping differential expression analysis. Use --force to overwrite.", filename_markers)
        markers <- readRDS(filename_markers)
    } else {
        wrong_assays <- c("SCT", "integrated")
        preferred_assays <- c("Spatial", "RNA") 
        if (DefaultAssay(obj) %in% wrong_assays) {
            for (assay in Assays(obj)) {
                if (assay %in% preferred_assays) {
                    DefaultAssay(obj) <- assay
                    flog.info("Setting default assay to %s for FindMarkers...", assay)
                    break
                }
            }
        }
        markers <- FindAllMarkers(obj, ...)
        flog.info("Writing R data structure to %s...", filename_markers)
        .serialize(markers, opt$outprefix, suffix)
    }
    return(markers)
}

