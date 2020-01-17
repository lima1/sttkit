#' read_spatial
#'
#' Read spatial transcriptomics into Seurat object
#' @param file File or data containing Spatial read counts
#' @param sampleid \code{character(1)} with sampleid
#' @param mt_pattern Regex to match mitochondria genes
#' @param rp_pattern Regex to match ribosomal genes
#' @param min_features Minimum number of features to include spot
#' @param min_spots Minimum number of spots to include gene
#' @param required_features Add those features even if not present or
#' passing \code{min_spots}
#' @param transpose Transpose input matrix
#' @param barcodes Optionally, barcodes for all spots
#' @param plot_qc Generate QC plots
#' @param serialize Automatically serialize object
#' @param prefix Prefix of output files
#' @import Seurat
#' @importFrom futile.logger flog.info flog.warn
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom methods is
#' @importFrom utils read.delim write.csv write.table
#' @importFrom Matrix t rowSums
#' @export read_spatial
#' @examples
#' read_spatial()
read_spatial <- function(file, sampleid, mt_pattern = regex_mito(), 
                        rp_pattern = regex_ribo(), 
                        min_features = 300, min_spots = 2, required_features = NULL, 
                        transpose = FALSE, barcodes = NULL,
                        plot_qc = TRUE, serialize = TRUE, prefix) {
    if (is(file, "character")) {
        flog.info("Loading %s...", basename(file))
        raw_data <- Matrix::t(read.delim(file, row.names=1))
    } else {
        raw_data <- Matrix::t(file)
    }    
    if (transpose) raw_data <- Matrix::t(raw_data)
    spots_passing <- sum(Matrix::rowSums(raw_data) > min_features)    
    if (spots_passing < 100) {
        flog.warn("Less than 100 spots passing min_features %i (%i out of %i).", 
            min_features, spots_passing, ncol(raw_data))
    }
    if (sum(grepl("^mm|^hs", rownames(raw_data))) > nrow(raw_data)*0.9) { 
        rownames(raw_data) <- gsub("^mm", "mm10", rownames(raw_data))
        rownames(raw_data) <- gsub("^hs", "hg19", rownames(raw_data))
    }
    if (!is.null(required_features)) {
        required_features <- make.names(required_features)
        required_features <- required_features[!required_features %in% rownames(raw_data)]
        if (length(required_features)) {
            raw_data_required <- matrix(0, nrow = length(required_features),
                ncol = ncol(raw_data), 
                dimnames=list(required_features, colnames(raw_data)))
            raw_data <- rbind(raw_data, raw_data_required)
        }
        raw_data <- raw_data[Matrix::rowSums(raw_data) >= min_spots | 
                             rownames(raw_data) %in% required_features, ]
    } else {
        raw_data <- raw_data[Matrix::rowSums(raw_data) >= min_spots, ]
    }    
    ndata <- CreateSeuratObject(raw_data, min.cells = 0,
        min.features = min_features, project = sampleid)
    mito.features <- grep(pattern = mt_pattern, x = rownames(x = ndata), value = TRUE)
    ndata <- PercentageFeatureSet(object = ndata, pattern = mt_pattern, col.name = "percent.mito")
    ndata <- PercentageFeatureSet(object = ndata, pattern = rp_pattern, col.name = "percent.ribo")

    if (plot_qc) {
        pdf(paste0(prefix, "_qc.pdf"))
        print(VlnPlot(object = ndata, 
            features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"),
            ncol = 2))
        par(mfrow = c(2, 2))
        print(FeatureScatter(object = ndata, feature1 = "nFeature_RNA", feature2 = "nCount_RNA"))
        print(FeatureScatter(object = ndata, feature1 = "nFeature_RNA", feature2 = "percent.ribo"))
        if (length(unique(ndata$percent.mito))>1) { 
            print(FeatureScatter(object = ndata, feature1 = "nFeature_RNA", feature2 = "percent.mito"))
            print(FeatureScatter(object = ndata, feature1 = "percent.mito", feature2 = "percent.ribo"))

        }
        dev.off()
        .write_qc_stats(ndata, prefix)
    }    
    ndata@meta.data$library <- sampleid
    if (!is.null(barcodes)) {
        ndata <- AddMetaData(object = ndata, metadata = barcodes, col.name = "barcode")
    }    
    cnts <- GetAssayData(object = ndata, slot = 'counts')

    if (sum(grepl("^hg19", rownames(ndata)))>50) {
        hg19.features <- grep(pattern = "^hg19", x = rownames(x = ndata), value = TRUE)
        percent.hg19 <- Matrix::colSums(x = cnts[hg19.features, ]) / Matrix::colSums(x = cnts)
        ndata <- AddMetaData(object = ndata, metadata = percent.hg19, col.name = "hg19")
        ndata$nFeature_RNA_hg19 <- apply(as.matrix(cnts[grep("hg19", rownames(cnts)),]),2,function(x) length(which(x>0)))
    }    
    if (sum(grepl("^mm10", rownames(ndata)))>50) {
        mm10.features <- grep(pattern = "^mm10", x = rownames(x = ndata), value = TRUE)
        percent.mm10 <- Matrix::colSums(x = cnts[mm10.features, ]) / Matrix::colSums(x = cnts)
        ndata <- AddMetaData(object = ndata, metadata = percent.mm10, col.name = "mm10")
        ndata$nFeature_RNA_mm10 <- apply(as.matrix(cnts[grep("mm10", rownames(cnts)),]),2,function(x) length(which(x>0)))
    }
    if (serialize) .serialize(ndata, prefix, "_raw.rds")
    ndata
}

.write_qc_stats <- function(ndata, prefix) {
    stats_mean <- suppressWarnings(apply(ndata@meta.data, 2, 
        function(x) mean(as.numeric(x),na.rm=TRUE)))
    stats_median <- suppressWarnings(apply(ndata@meta.data, 2, 
        function(x) median(as.numeric(x),na.rm=TRUE)))
    names(stats_mean) <- paste0(names(stats_mean), ".mean")
    names(stats_median) <- paste0(names(stats_median), ".median")
    stats <- c(stats_mean, stats_median)
    stats <- stats[order(names(stats))]
    stats <- t(as.matrix(stats[is.finite(stats)]))
    write.csv(stats, file = paste0(prefix, "_qc.csv"), row.names = FALSE)
}

#' read_visium
#'
#' Read 10X SpaceRanger data into Seurat object
#' @param filtered_feature_bc_matrix_dir Path to SpaceRanger filtered matrix
#' @param spatial_dir Path to SpaceRanger \code{spatial} directory
#' @param tissue_positions_file File providing the coordindates
#' @param ... Arguments passed to \code{\link{read_spatial}}
#' @export read_visium
#' @examples
#' read_visium()
read_visium <- function(filtered_feature_bc_matrix_dir, spatial_dir, tissue_positions_file = "tissue_positions_list.csv", ...) {
    raw_data <- Read10X(filtered_feature_bc_matrix_dir)
    spatial <- read.csv(file.path(spatial_dir, tissue_positions_file),
        header = FALSE, as.is = TRUE)
    idx <- match(colnames(raw_data), gsub("-\\d+$", "", spatial[,1]))
    barcodes <- colnames(raw_data)
    colnames(raw_data) <- paste0(spatial$V6[idx], "x", spatial$V5[idx])
    names(barcodes) <- colnames(raw_data)
    read_spatial(Matrix::t(raw_data), barcodes = barcodes, ...)
} 

.serialize <- function(x, prefix, suffix) {
    filename <- .get_serialize_path(prefix, suffix)
    saveRDS(x, filename)     
}
.get_serialize_path <- function(prefix, suffix) {
    .get_sub_path(prefix, "serialize", suffix)
}
.get_advanced_path <- function(prefix, suffix) {
    .get_sub_path(prefix, "advanced", suffix)
}
.get_sub_path <- function(prefix, sub, suffix) {
    s_dir <- file.path(dirname(prefix), sub)
    if (!dir.exists(s_dir)) dir.create(s_dir)
    file.path(s_dir, paste0(basename(prefix), suffix))
}    
