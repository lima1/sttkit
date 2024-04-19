suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help = "Infile SpatialTranscriptomics data tsv file. Rows spots, columns genes. "),
    make_option(c("--spaceranger_dir"), action = "store", type = "character", default = NULL,
        help = "Path to SpaceRanger output for Visum data."),
    make_option(c("--spaceranger_probe_set"), action = "store", type = "character", default = NULL,
        help = "Path to SpaceRanger probe set file for feature flagging. Only useful when probe set filtering was turned off."),
    make_option(c("--sampleid"), action = "store", type = "character", default = NULL,
        help = "Sample id."),
    make_option(c("--transpose"), action = "store_true", default = FALSE,
        help = "In case data is stored as rows genes, spots columns."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help = "Outfile."),
    make_option(c("--num_features"), action = "store", type = "integer", default = 3000,
        help = "Calculate that many variable features [default %default]"),
    make_option(c("--regressout"), action = "store", type = "character",
        default = NULL,
        help = "Variables to regress out [default %default]"),
    make_option(c("--normalization_method"), action = "store", default = "sctransform",
        help = "Which normalization method to use, seurat, sctransform, sctransform2, or scran [default %default]."),
    make_option(c("--normalization_backend_method"), action = "store", default = "poisson",
        help = "Which backend normalization method to use, for example poisson or glmGamPoi [default %default]. Only relevant for sctransform and sctransform2."),
    make_option(c("--hejpeg"), action = "store", type = "character", default = NULL,
        help = "Optional path to a JPEG containing cropped HE image (Spatial Transcriptomics data)."),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.5,
        help = "Size of dots on H&E."),
    make_option(c("--min_spots"), action = "store", type = "integer", default = 2,
        help = "Gene filter: Keep genes detected at that many spots or more [default %default]"),
    make_option(c("--min_features"), action = "store", type = "double", default = 300,
        help = "Spot filter: Keep spots that detected that many genes or more [default %default]"),
    make_option(c("--downsample_prob"), action = "store", type = "double", default = NULL,
        help = "Downsample count matrix. 0.2 randomly picks 20 percent of UMIs in each spot [default %default]"),
    make_option(c("--gmt"), action = "store", type = "character",
        default = NULL,
        help = "GMT file including genes of interest"),
    make_option(c("--gtf"), action = "store", type = "character", 
        default = NULL, 
        help = "Optional GTF for including feature meta data like alternative gene ids."),
    make_option(c("--mane"), action = "store", type = "character", 
        default = NULL, 
        help = "Optional NCBI MANE file for resolving ambiguous mapping of gene name and id in --GTF."),
    make_option(c("--output_counts"), action = "store_true", default = FALSE,
        help = "Output count matrix as TSV file."),
    make_option(c("--no_crop"), action = "store_true", default = FALSE,
        help = "Do not crop H&E image."),
    make_option(c("--image_formats"), action = "store", type = "character", 
        default = "png", 
        help = "Image format(s) of output plots. 'png' and 'pdf' supported. Multiple formats are seperated by colon ('png:pdf')."),
    make_option(c("--png"), action = "store_true", default = FALSE,
        help = "Generate PNG version of output plots. DEPRECATED."),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE,
        help = "Overwrite existing files")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$infile) && is.null(opt$spaceranger_dir)) {
    stop("Need --infile or --spaceranger_dir")
}
if (is.null(opt$outprefix)) {
    stop("Need --outprefix")
}
if (opt$png) {
    flog.warn("--png is deprecated.")
    if (is.null(opt$image_formats)) {
        # old default
        opt$image_formats <- "pdf:png"
    }
}    
if (is.null(opt$image_formats)) {
    opt$image_formats <- c()
} else {
    opt$image_formats <- sapply(strsplit(opt$image_formats, ":")[[1]], tolower)
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

.plot_he_ptx <- function(ndata, prefix, assay = "Spatial") {
    flog.info("Plotting counts on H&E...")
    if (sum(grep("^hg19", rownames(ndata)))) {
        features <- c("mm10", "hg19")
        if (paste0("nFeature_", assay, "_mm10") %in% colnames(ndata@meta.data)) {
            features <- c(features,
                          paste0("nFeature_", assay, "_mm10"),
                          paste0("nFeature_", assay, "_hg19"))
        }

        plot_features(ndata, features = features,
            pt.size.factor = opt$dot_size,
            prefix = prefix, suffix = "_he_ptx.pdf",
            ggcode = theme(legend.position = "right"),
            pdf = "pdf" %in% opt$image_formats,
            png = "png" %in% opt$image_formats,
            labels = scales::percent, width = 8, height = 4,
            crop = !opt$no_crop)
    }

    features <- c(paste0("nCount_", assay),
                  paste0("nFeature_", assay),
                  "percent.mito",
                  "percent.ribo")
    plot_features(ndata, features = features,
        pt.size.factor = opt$dot_size,
        prefix = prefix, suffix = "_he_counts.pdf",
        ggcode = theme(legend.position = "right"),
        pdf = "pdf" %in% opt$image_formats,
        png = "png" %in% opt$image_formats,
        crop = !opt$no_crop,
        labels = scales::percent, width = 8, height = 4)
    if (sum(grep("S.Score|G2M.Score", colnames(ndata@meta.data)))) {
        plot_features(ndata, features = c("S.Score", "G2M.Score"),
            pt.size.factor = opt$dot_size,
            prefix = prefix, suffix = "_he_cell_cycle.pdf",
            ggcode = theme(legend.position = "right"),
            pdf = "pdf" %in% opt$image_formats,
            png = "png" %in% opt$image_formats,
            crop = !opt$no_crop,
            width = 8, height = 4)
    }
}
.get_serialize_path <- function(prefix, suffix) {
    s_dir <- file.path(dirname(prefix), "serialize")
    file.path(s_dir, paste0(basename(prefix), suffix))
}

.write_tsv <- function(object, prefix, slot = "scale.data", suffix = "_all.tsv.gz") {
    data <- GetAssayData(object, slot = "scale.data")
    filename <- paste0(prefix, suffix)
    flog.info("Writing data to %s...", basename(filename))
    m <- t(as.matrix(data))
    pos <- GetTissueCoordinates(object)[rownames(m), ]
    data.table::fwrite(data.table::as.data.table(data.frame(barcode = rownames(m), pos, m)), file = filename,
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
                             image = opt$hejpeg,
                             transpose = opt$transpose,
                             sampleid = opt$sampleid,
                             downsample_prob = opt$downsample_prob,
                             gtf = opt$gtf,
                             mane = opt$mane,
                             prefix = opt$outprefix)
    } else {
        ndata <- read_visium(
                             opt$spaceranger_dir,
                             probe_set = opt$spaceranger_probe_set,
                             min_spots = opt$min_spots,
                             min_features = opt$min_features,
                             required_features = required_features,
                             transpose = opt$transpose,
                             sampleid = opt$sampleid,
                             downsample_prob = opt$downsample_prob,
                             gtf = opt$gtf,
                             mane = opt$mane,
                             prefix = opt$outprefix)
    }
    if (opt$output_counts) {
        m <- .write_tsv(ndata, opt$outprefix, slot = "counts",
            suffix = "_all_counts.tsv.gz")
    }
    ndata <- normalize_spatial(ndata, nfeatures = opt$num_features,
                         scale = opt$normalization_method != "scran",
                         center = opt$normalization_method != "scran",
                         regressout = .parse_regressout(),
                         method = opt$normalization_method,
                         backend_method = opt$normalization_backend_method,
                         prefix = opt$outprefix)
}

plot_qc_read(ndata, opt$outprefix, assay = Assays(ndata))
.plot_he_ptx(ndata, opt$outprefix, assay = Assays(ndata))

m <- .write_tsv(ndata, opt$outprefix)

.plot_he_scran_cluster <- function(ndata, prefix) {
    flog.info("Plotting clusters on H&E...")
    plot_features(ndata, features = "scran.cluster",
        labels = waiver(), labels_title = "",
        prefix = prefix,
        suffix = "_he_scran_cluster.pdf",
        pdf = "pdf" %in% opt$image_formats,
        png = "png" %in% opt$image_formats,
        width = 4,
        crop = !opt$no_crop,
        pt.size.factor = opt$dot_size)
}

if (opt$normalization_method == "scran") {
    .plot_he_scran_cluster(ndata, opt$outprefix)
}
flog.info("Done.")
