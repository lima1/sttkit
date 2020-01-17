suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help="Infile SpatialTranscriptomics data tsv file. Rows spots, columns genes. "),
    make_option(c("--spaceranger_dir"), action = "store", type = "character", default = NULL,
        help="Path to SpaceRanger output for Visum data."),
    make_option(c("--sampleid"), action = "store", type = "character", default = NULL,
        help="Sample id."),
    make_option(c("--transpose"), action = "store_true", default = FALSE, 
        help="In case data is stored as rows genes, spots columns."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--num_features"), action="store", type = "integer", default = 2500, 
        help="Calculate that many variable features [default %default]"),
    make_option(c("--regressout"), action = "store", type = "character",
        default = NULL, 
        help="Variables to regress out [default %default]"),
    make_option(c("--normalization_method"), action = "store", default = "sctransform", 
        help="Which normalization method to use, seurat, sctransform or scran [default %default]."),
    make_option(c("--hejpeg"), action = "store", type = "character", default = NULL,
        help="Optional path to a JPEG containing cropped HE image."),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.5,
        help="Size of dots on H&E."),
    make_option(c("--min_spots"), action="store", type = "integer", default = 2, 
        help="Gene filter: Keep genes detected at that many spots or more [default %default]"),
    make_option(c("--min_features"), action="store", type = "double", default = 300, 
        help="Spot filter: Keep spots that detected that many genes or more [default %default]"),
    make_option(c("--gmt"), action = "store", type = "character", 
        default = NULL, 
        help="GMT file including genes of interest"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE, 
        help="Overwrite existing files")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$infile) && is.null(opt$spaceranger_dir)) {
    stop("Need --infile or --spaceranger_dir")
}
if (is.null(opt$outprefix)) {
    stop("Need --outprefix")
}
if (!dir.exists(dirname(opt$outprefix))) {
    dir.create(dirname(dirname(opt$outprefix)))
    dir.create(dirname(opt$outprefix))
}
    
flog.info("Loading Seurat...")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sttkit))
log_file <- paste0(opt$outprefix, "_normalize.log")
if (!is.null(log_file)) flog.appender(appender.tee(log_file))

.plot_he_ptx <- function(ndata, prefix, hejpeg) {
    flog.info("Plotting counts on H&E...")
    if (sum(grep("^hg19", rownames(ndata)))) {
        filename <- paste0(prefix, "_he_ptx.pdf")
        pdf(filename, width=8, height=3.9)
        plot_features(ndata, features = c("mm10", "hg19"), 
            hejpeg = hejpeg,  labels = scales::percent, 
            reorder_clusters = FALSE, size = opt$dot_size,
            plot_correlations = FALSE)
        if ("nFeature_RNA_mm10" %in% colnames(ndata@meta.data)) {
            plot_features(ndata, 
                features = c("nFeature_RNA_mm10", "nFeature_RNA_hg19"), 
                hejpeg = hejpeg,  labels = function(x) sprintf("%.0f", x), 
                labels_title ="", trans = TRUE, reorder_clusters = FALSE,
                plot_correlations = FALSE)
        }
        dev.off()
    }     
    filename <- paste0(prefix, "_he_counts.pdf")
    pdf(filename, width = 4, height = 3.6)
    plot_features(ndata, features = c("nCount_RNA"), 
        hejpeg = hejpeg,  labels = function(x) sprintf("%.0f", x), labels_title = "", trans = TRUE,
        reorder_clusters = FALSE)
    plot_features(ndata, features = c("nFeature_RNA"), 
        hejpeg = hejpeg,  labels = function(x) sprintf("%.0f", x), labels_title = "", trans = TRUE,
        reorder_clusters = FALSE)
    dev.off()
}
.get_serialize_path <- function(prefix, suffix) {
    s_dir <- file.path(dirname(prefix), "serialize")
    file.path(s_dir, paste0(basename(prefix), suffix))
}

.write_tsv <- function(object, prefix) {
    data <- GetAssayData(object, slot = "scale.data") 
    filename <- paste0(prefix, "_all.tsv.gz")
    flog.info("Writing normalized data to %s...", basename(filename))
    m <- t(as.matrix(data))
    data.table::fwrite(data.table::as.data.table(m), file = filename,
        sep = "\t", quote = FALSE)
    m
}

.parse_regressout <- function() { 
    if (is.null(opt$regressout)) NULL else strsplit(opt$regressout, ":")[[1]]
}

filename <- .get_serialize_path(opt$outprefix, "_scaled.rds")
if (!opt$force && file.exists(filename)) {
    flog.warn("%s exists. Skipping normalization. Use --force to overwrite.", filename)
    ndata <- readRDS(filename)
} else {
    required_features <- NULL
    if (!is.null(opt$gmt)) {
        gmt <- read_signatures(opt$gmt)
        required_features <- Reduce(union, gmt)
    }    
    if (!is.null(opt$infile)) {
        ndata <- read_spatial(opt$infile, min_spots = opt$min_spots, 
                             min_features = opt$min_features,
                             required_features = required_features, 
                             transpose = opt$transpose,
                             sampleid = opt$sampleid,
                             prefix = opt$outprefix)
    } else {
        ndata <- read_visium(
                             file.path(opt$spaceranger_dir, "filtered_feature_bc_matrix"),
                             file.path(opt$spaceranger_dir, "spatial"),
                             min_spots = opt$min_spots, 
                             min_features = opt$min_features,
                             required_features = required_features, 
                             transpose = opt$transpose,
                             sampleid = opt$sampleid,
                             prefix = opt$outprefix)
    }    

    ndata <- normalize_spatial(ndata, nfeatures = opt$num_features, 
                         scale = opt$normalization_method != "scran",
                         center = opt$normalization_method != "scran",
                         regressout = .parse_regressout(), 
                         method = opt$normalization_method,
                         prefix = opt$outprefix)
}

.plot_he_ptx(ndata, opt$outprefix, opt$hejpeg)

m <- .write_tsv(ndata, opt$outprefix)

.plot_he_scran_cluster <- function(ndata, prefix, hejpeg) {
    filename <- paste0(prefix, "_he_scran_cluster.pdf")
    flog.info("Plotting clusters on H&E...")
    pdf(filename, width = 4, height = 3.9)
    plot_features(ndata, features = "scran.cluster", hejpeg = hejpeg, 
        labels = waiver(), labels_title = "", 
        reorder_clusters = FALSE, size = opt$dot_size)
    dev.off()
}

if (opt$normalization_method == "scran") {
    .plot_he_scran_cluster(ndata, opt$outprefix, opt$hejpeg)
}
