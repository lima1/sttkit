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
#' @param hybrid_reference_prefixes For dataset aligned to hybrid references, the
#' the feature name prefixes. First one reserved for human (for non-human,
#' hybrid references, this just means that some downstream tools might not work).
#' @param gtf Optional GTF for including feature meta data like alternative gene ids.
#' @param mane Optional NCBI MANE file for resolving ambiguous mapping in \code{gtf} 
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
#' @importFrom data.table fread
#' @export read_spatial
#' @examples
#' read_spatial()
read_spatial <- function(file, sampleid, mt_pattern = regex_mito(), 
                        rp_pattern = regex_ribo(), 
                        min_features = 300, min_spots = 2,
                        required_features = NULL, 
                        transpose = FALSE, barcodes = NULL, image = NULL, slice = sampleid,
                        downsample_prob = NULL,
                        hybrid_reference_prefixes = c("hg19", "mm10"),
                        gtf = NULL, mane = NULL, assay = "Spatial", 
                        serialize = TRUE, prefix) {
    ndata <- NULL
    if (is(file, "character")) {
        flog.info("Loading %s...", basename(file))
        raw_data <- Matrix::t(read.delim(file, row.names=1))
    } else if (is(file, "Seurat")) {    
        ndata <- file
    } else {
        raw_data <- Matrix::t(file)
    }
    if (is.null(ndata)) {
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
        # TODO: this hack is obsolete with hybrid_reference_prefixes
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
    }
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
    ndata <- .annotate_features_gtf(ndata, assay, gtf, mane)

    cnts <- GetAssayData(object = ndata, slot = 'counts')
    for (hrp in hybrid_reference_prefixes) {
        pattern <- paste0("^", hrp)
        if (sum(grepl(pattern, rownames(ndata))) > 50) {
            features <- grep(pattern = pattern, x = rownames(x = ndata), value = TRUE)
            percent <- Matrix::colSums(x = cnts[features, ]) / Matrix::colSums(x = cnts)
            ndata <- AddMetaData(object = ndata, metadata = percent, col.name = hrp)
            key <- paste0("nFeature_", DefaultAssay(ndata), "_", hrp)
            ndata[[key]] <- apply(as.matrix(cnts[grep(hrp, rownames(cnts)), ]),
                2, function(x) length(which(x > 0)))
        }
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
#' @param bin_size Bin size for VisiumHD
#' @param assay Name of the assay corresponding to the initial input data.
#' @param probe_set Path to SpaceRanger probe set file. Only useful
#' if SpaceRanger was run without probe set filtering.
#' @param ... Arguments passed to \code{\link{read_spatial}}
#' @export read_visium
#' @examples
#' read_visium()
read_visium <- function(filtered_feature_bc_matrix_dir,
    spatial_dir = file.path(filtered_feature_bc_matrix_dir, "spatial"),
    bin_size = NULL,
    assay = "Spatial", probe_set = NULL, ...) {
    requireNamespace("hdf5r")
    
    # VisiumHD?
    if (dir.exists(file.path(filtered_feature_bc_matrix_dir, "binned_outputs"))) {
        if (is.null(bin_size)) {
            bin_size <- 8
            flog.warn("bin_size not specified, defaulting to %i um.", bin_size)
        } 
        bin_size_pretty <- paste0(sprintf("%03d", bin_size), "um")
        filtered_feature_bc_matrix_dir <- file.path(filtered_feature_bc_matrix_dir, "binned_outputs", paste0("square_", 
                             bin_size_pretty))
        spatial_dir <- file.path(filtered_feature_bc_matrix_dir, "spatial")
    } else if (!dir.exists(spatial_dir)) {
        flog.warn("%s does not exist.", spatial_dir)
        if (file.exists(file.path(filtered_feature_bc_matrix_dir, 
            "scalefactors_json.json"))) {
            spatial_dir <- filtered_feature_bc_matrix_dir
        }
    }
    filename <- file.path(filtered_feature_bc_matrix_dir, 
        "filtered_feature_bc_matrix.h5")
    flog.info("Reading h5 feature matrix %s.", filename)
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

    ndata <- read_spatial(Matrix::t(raw_data), barcodes = barcodes, image = image, ...)
    ndata <- read_spaceranger_deconvolution(ndata, filtered_feature_bc_matrix_dir)
    metrics <- read_spaceranger_metrics(filtered_feature_bc_matrix_dir)
    if (!is.null(probe_set)) {
        probes <- read_spaceranger_probe_set(probe_set)        
        # if ndata contains more stable gene_id, use that one for mapping probes to features 
        if ("gene_id" %in% colnames(ndata[[assay]][[]])) {
            idx <- !rownames(ndata[[assay]]) %in% rownames(probes)
            map_missing <- ndata[[assay]][[]][idx, , drop = FALSE]
            map_missing <- map_missing[!is.na(map_missing$gene_id), , drop = FALSE]
            probes_missing <- probes[match(map_missing$gene_id, probes$gene_id), ]
            rownames(probes_missing) <- rownames(map_missing)
            probes <- rbind(probes, probes_missing)
        }
        ndata[[assay]] <- AddMetaData(ndata[[assay]], probes)
    }
    return(ndata)
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
        if (grepl("hg19|hg38|GRCh3", prefix[1])) {
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

#' read_spaceranger_metrics
#'
#' Parse 10X SpaceRanger metrics
#' @param filtered_feature_bc_matrix_dir Path to SpaceRanger filtered matrix
#' @param metrics_file Path to SpaceRanger \code{metrics_summary.csv} file
#' @export read_spaceranger_metrics
#' @examples
#' #read_spaceranger_metrics
read_spaceranger_metrics <- function(filtered_feature_bc_matrix_dir,
        metrics_file = file.path(filtered_feature_bc_matrix_dir, "metrics_summary.csv")) {

    if (!file.exists(metrics_file)) {
        flog.warn("Not finding %s. Skipping check of those numbers.", metrics_file)
        return()
    }
    metrics <- read.csv(metrics_file)
    if (!is.null(metrics[["Fraction.Reads.in.Spots.Under.Tissue"]]) &&
                 metrics[["Fraction.Reads.in.Spots.Under.Tissue"]] < 0.5) {
        flog.warn("Low 'Fraction Reads in Spots Under Tissue' %.2f. Ideal > 0.5. Application performance may be affected. Many of the reads were not assigned to tissue covered spots. This could be caused by high levels of ambient RNA resulting from inefficient permeabilization, because the incorrect image was used, or because of poor tissue detection. The latter case can be addressed by using the manual tissue selection option through Loupe.", metrics[["Fraction.Reads.in.Spots.Under.Tissue"]])
    } 
    return(metrics)
}

#' read_spaceranger_probe_set
#'
#' Parse 10X SpaceRanger probe sets
#' @param file Path to SpaceRanger probe sets file
#' @export read_spaceranger_probe_set
#' @examples
#' #read_spaceranger_probe_set
read_spaceranger_probe_set <- function(file) {

    if (!file.exists(file)) {
        flog.warn("Not finding %s. Skipping flagging of features.", file)
        return(NULL)
    }
    probes <- read.csv(file, comment.char = "#")
    expected_cols <- c("gene_id", "probe_id",  "included")
    if (!all(expected_cols %in% colnames(probes))) {
        flog.warn("Unrecognized probe set file %s.", file)
        return(NULL)
    }

    probes$symbol <- sapply(strsplit(probes$probe_id, "\\|"), function(x) x[2])

    px <- split(probes, probes$symbol)
    
    probes_by_symbol <- do.call(rbind, lapply(px, function(x) {
        data.frame(
            gene_id = x$gene_id[1],
            symbol = x$symbol[1],
            probe_seqs = paste(x$probe_seq, collapse = "|"),
            all.included = paste(x$included, collapse = "|"),
            included = any(x$included),
            regions = paste(x$region, collapse = "|"),
            row.names = x$symbol[1])}))

    return(probes_by_symbol)
}

.annotate_features_gtf <- function(ndata, assay, gtf, mane) {
    if (!is.null(gtf)) {
        if (requireNamespace("rtracklayer")) {
            flog.info("Loading gtf %s...", basename(gtf))    
            gtf_ref <- rtracklayer::mcols(rtracklayer::import(gtf))
            col_feature <- names(which.max(apply(gtf_ref,2,function(x) length(intersect(rownames(ndata), x)))))
            if (col_feature == "gene_id") {
                flog.info("Skipping gene_id annotation because features appear to be gene ids, not gene names.")
            } else {
                if ("gene_id" %in% colnames(gtf_ref)) {
                    idx <- !duplicated(paste(gtf_ref[[col_feature]], gtf_ref$gene_id))
                    idx2 <- duplicated(gtf_ref[[col_feature]][idx])
                    if (any(idx2) && !is.null(mane)) {
                        if (is(mane, "character")) {
                            mane <- fread(mane)
                        }    
                        dup_symbols <- gtf_ref$gene_name[idx][idx2]
                        dup_ids <- unique(gtf_ref$gene_id[gtf_ref$gene_name %in% dup_symbols])
                        ignore_ids <- dup_ids[!sapply(dup_ids, function(x) any(grepl(x, mane$Ensembl_Gene)))]
                        idx <- idx & !gtf_ref[["gene_id"]] %in% ignore_ids
                        idx2 <- duplicated(gtf_ref[[col_feature]][idx])
                    }
                    if (any(idx2)) {
                        flog.warn("Ambiguous 'gene_name' and 'gene_id' mapping for following features: %s",
                            paste(gtf_ref[[col_feature]][idx][idx2], collapse = ","))
                    }    
                    probes <- data.frame(
                        gene_id = gtf_ref[["gene_id"]][idx][!idx2],
                        row.names = gtf_ref[[col_feature]][idx][!idx2]
                    )
                    ndata[[assay]] <- AddMetaData(ndata[[assay]], probes)
                } else {
                    flog.info("Skipping gene_id annotation because %s does not contain required 'gene_name' and 'gene_id' fields.", gtf)
                }
            }
        } else {
            flog.warn("Install rtracklayer for parsing GTF file.")
        }    
    }    
    return(ndata)
}    
