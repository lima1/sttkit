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
#' @param image Optional \code{VisiumV1} object containing image information.
#' Can also be a JPEG for Spatial Transcriptomics data.
#' @param slice Optional name for the stored \code{image}
#' @param downsample_prob Downsample input matrix. Requires \code{DropletUtils}.
#' @param assay Name of the assay corresponding to the initial input data.
#' @param serialize Automatically serialize object
#' @param prefix Prefix of output files
#' @import Seurat
#' @importFrom futile.logger flog.info flog.warn
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom methods is new
#' @importFrom utils read.delim write.csv write.table data
#' @importFrom Matrix t rowSums
#' @export read_spatial
#' @examples
#' read_spatial()
read_spatial <- function(file, sampleid, mt_pattern = regex_mito(), 
                        rp_pattern = regex_ribo(), 
                        min_features = 300, min_spots = 2, required_features = NULL, 
                        transpose = FALSE, barcodes = NULL, image = NULL, slice = sampleid,
                        downsample_prob = NULL,
                        assay = "Spatial", 
                        serialize = TRUE, prefix) {
    if (is(file, "character")) {
        flog.info("Loading %s...", basename(file))
        raw_data <- Matrix::t(read.delim(file, row.names=1))
    } else {
        raw_data <- Matrix::t(file)
    }    
    if (transpose) raw_data <- Matrix::t(raw_data)
    if (!is.null(downsample_prob)) {
        if (!requireNamespace("DropletUtils")) {
            flog.warn("downsample_prob defined but DropletUtils not installed. Ignoring.")
        }
        flog.info("Downsample input to %.2f...", downsample_prob) 
        raw_data <- DropletUtils::downsampleMatrix(raw_data, prop = downsample_prob)
    }        
    spots_passing <- sum(Matrix::rowSums(raw_data) > min_features)    
    nfeatures_0 <- Matrix::colSums(x = raw_data > 0)

    if (spots_passing < 100) {
        flog.warn("Less than 100 spots passing min_features %i (%i out of %i).", 
            min_features, spots_passing, ncol(raw_data))
    }
    if (sum(grepl("^mm|^hs", rownames(raw_data))) > nrow(raw_data)*0.9) { 
        rownames(raw_data) <- gsub("^mm_|^mm-", "mm10-", rownames(raw_data))
        rownames(raw_data) <- gsub("^hs_|^hs-", "hg19-", rownames(raw_data))
    }
    if (!is.null(required_features)) {
        required_features <- required_features[!make.names(required_features) %in% make.names(rownames(raw_data))]
        if (length(required_features)) {
            raw_data_required <- matrix(0, nrow = length(required_features),
                ncol = ncol(raw_data), 
                dimnames=list(required_features, colnames(raw_data)))
            raw_data <- rbind(raw_data, raw_data_required)
        }
        raw_data <- raw_data[Matrix::rowSums(raw_data) >= min_spots | 
                             make.names(rownames(raw_data)) %in% make.names(required_features), ]
    } else {
        raw_data <- raw_data[Matrix::rowSums(raw_data) >= min_spots, ]
    }    
    ndata <- CreateSeuratObject(raw_data, min.cells = 0,
        min.features = min_features, project = sampleid, assay = assay)
    mito.features <- grep(pattern = mt_pattern, x = rownames(x = ndata), value = TRUE)
    ndata <- PercentageFeatureSet(object = ndata, pattern = mt_pattern, col.name = "percent.mito")
    ndata <- PercentageFeatureSet(object = ndata, pattern = rp_pattern, col.name = "percent.ribo")

    ndata@meta.data$library <- sampleid
    if (!is.null(barcodes)) {
        ndata <- AddMetaData(object = ndata, metadata = barcodes, col.name = "barcode")
    }    
    if (!is.null(image)) {
        if (class(image) == "character") {
            image <- .read_spatial_transcriptomics_image(ndata, image)
        }
        image <- image[Cells(x = ndata)]
        DefaultAssay(object = image) <- assay
        ndata[[slice]] <- image
    }
    cnts <- GetAssayData(object = ndata, slot = 'counts')

    if (sum(grepl("^hg19", rownames(ndata)))>50) {
        hg19.features <- grep(pattern = "^hg19", x = rownames(x = ndata), value = TRUE)
        percent.hg19 <- Matrix::colSums(x = cnts[hg19.features, ]) / Matrix::colSums(x = cnts)
        ndata <- AddMetaData(object = ndata, metadata = percent.hg19, col.name = "hg19")
        key <- paste0("nFeature_", DefaultAssay(ndata), "_hg19")
        ndata[[key]] <- apply(as.matrix(cnts[grep("hg19", rownames(cnts)),]),2,function(x) length(which(x>0)))
    }    
    if (sum(grepl("^mm10", rownames(ndata)))>50) {
        mm10.features <- grep(pattern = "^mm10", x = rownames(x = ndata), value = TRUE)
        percent.mm10 <- Matrix::colSums(x = cnts[mm10.features, ]) / Matrix::colSums(x = cnts)
        ndata <- AddMetaData(object = ndata, metadata = percent.mm10, col.name = "mm10")
        key <- paste0("nFeature_", DefaultAssay(ndata), "_mm10")
        ndata[[key]] <- apply(as.matrix(cnts[grep("mm10", rownames(cnts)),]),2,function(x) length(which(x>0)))
    }
    if (serialize) {
        flog.info("Writing R data structure to %s...", paste0(prefix, "_raw.rds"))
        .serialize(ndata, prefix, "_raw.rds")
    }
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
    stats <- c("total.spots" = ncol(ndata), "total.features" = nrow(ndata),
        stats)
    stats <- t(as.matrix(stats[is.finite(stats)]))
    write.csv(stats, file = paste0(prefix, "_qc.csv"), row.names = FALSE)
}

#' read_visium
#'
#' Read 10X SpaceRanger data into Seurat object
#' @param filtered_feature_bc_matrix_dir Path to SpaceRanger filtered matrix
#' @param spatial_dir Path to SpaceRanger \code{spatial} directory
#' @param assay Name of the assay corresponding to the initial input data.
#' @param ... Arguments passed to \code{\link{read_spatial}}
#' @export read_visium
#' @examples
#' read_visium()
read_visium <- function(filtered_feature_bc_matrix_dir, 
    spatial_dir = file.path(filtered_feature_bc_matrix_dir, "spatial"), 
    assay = "Spatial", ...) {
    requireNamespace("hdf5r")
    
    if (!dir.exists(spatial_dir)) {
        flog.warn("%s does not exist.", spatial_dir)
        if (file.exists(file.path(filtered_feature_bc_matrix_dir, 
            "scalefactors_json.json"))) spatial_dir <- filtered_feature_bc_matrix_dir
    }    
    filename <- file.path(filtered_feature_bc_matrix_dir, 
        "filtered_feature_bc_matrix.h5")
    raw_data <- Read10X_h5(filename = filename)
    lowres <- file.path(spatial_dir, "tissue_lowres_image.png")
    if (file.exists(lowres)) {
        image <- Read10X_Image(spatial_dir)
        DefaultAssay(object = image) <- assay
    } else {
        # check multi-sample
        lowres <- file.path(dir(spatial_dir, full.names = TRUE), "tissue_lowres_image.png")
        lowres <- lowres[file.exists(lowres)]
        image <- lapply(dirname(lowres), function(d) {
            image <- Read10X_Image(d)
            DefaultAssay(object = image) <- assay
        })    
    }    
    barcodes <- colnames(raw_data)
    names(barcodes) <- colnames(raw_data)

    read_spatial(Matrix::t(raw_data), barcodes = barcodes, image = image, ...)
} 

.read_spatial_transcriptomics_image <- function(obj, hejpeg) {
    if (!requireNamespace("jpeg", quietly = TRUE)) {
        stop("Please install packages jpeg SpatialTranscriptomics data.")
    }
    image <- jpeg::readJPEG(hejpeg)
    scale.factors <- data.frame(
        "fiducial_diameter_fullres" = 0.5,
        "tissue_hires_scalef" = 1,
        "tissue_lowres_scalef" = max(dim(x = image))/34
    )
    coords <- .parse_coords(obj, Cells(obj))
    tissue.positions <- data.frame( 
        "tissue" = 1,
        "row" = coords[,2], 
        "col" = coords[,1], 
        "imagerow" = coords[,2] - 1.5,
        "imagecol" = coords[,1] - 1.5,
        row.names = Cells(obj)
    )
    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * 
        scale.factors$tissue_lowres_scalef
    spot.radius <- unnormalized.radius/max(dim(x = image))
    return(new(Class = "VisiumV1", image = image, 
        scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
        fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, 
        scale.factors$tissue_lowres_scalef), coordinates = tissue.positions, 
        spot.radius = spot.radius))
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
.add_cc_score <- function(obj, ...) {
    data(cc.genes.all, envir = environment())
    prefix <- .check_for_symbol_prefix(obj)
    if (!is.null(prefix)) {
        if (grepl("hg19|hg38", prefix[1])) {
            cc.genes.all$Mixed <- lapply(cc.genes.all$Hs, function(x) c(paste0(prefix[1], "-", x)))
        }     
    }
    id <- which.max(sapply(cc.genes.all, function(x) sum(unlist(x) %in% rownames(obj))))
    s.features <- cc.genes.all[[id]]$s.genes
    g2m.features <- cc.genes.all[[id]]$g2m.genes
    
    if (!length(s.features) || !length(g2m.features)) {
        flog.warn("No cell cycle genes found.")
        return(obj)
    }
    flog.info("Adding cell cycle scoring...")
    objm <- try(CellCycleScoring(obj, s.features = s.features,
        g2m.features = cc.genes.all[[id]]$g2m.genes, set.ident = FALSE, ...))
    if (is(objm, "try-error")) return(obj)
    objm
}

.check_for_symbol_prefix <- function(obj) {
    ss <- strsplit(rownames(obj), "-")
    if (all(sapply(ss, length) > 1)) {
        return(sort(unique(sapply(ss, function(x) x[1]))))
    }
    return(NULL)
}        
