suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help="Infile Seurat RDS of st_normalize. Should be unscaled unless normalization is sctransform"),
    make_option(c("--singlecell"), action = "store", type = "character", default = NULL,
        help="Path to a RDS file containing a (list of) Seurat single cell object(s) for spot deconvolution."),
    make_option(c("--labels_singlecell"), action = "store", type = "character", default = NULL,
        help="Optional list of labels --singlecell"),
    make_option(c("--refdata"), action = "store", type = "character", default = "type",
        help="Meta data column with prediction labels in --singlecell"),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--num_integration_features"), action="store", type = "integer", default = 3000, 
        help="Integration: Use that many features [default %default]"),
    make_option(c("--simulation"), action = "store_true", 
        default = FALSE, 
        help="Subcluster the reference cells, specified by the call attribute [default %default]"),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.6,
        help="Size of dots on H&E [default %default]"),
    make_option(c("--png"), action = "store_true", default = FALSE, 
        help="Generate PNG version of output plots."),
    make_option(c("--verbose"), action = "store_true", default = FALSE, 
        help="Verbose output"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE, 
        help="Overwrite existing files")
)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


opt <- parse_args(OptionParser(option_list=option_list))


if (is.null(opt$infile)) {
    stop("Need --infile")
}
if (is.null(opt$outprefix)) {
    stop("Need --outprefix")
}
if (!dir.exists(dirname(opt$outprefix))) {
    dir.create(dirname(dirname(opt$outprefix)))
    dir.create(dirname(opt$outprefix))
}
if (is.null(opt$singlecell)) {
    stop("Need --singlecell")
}

flog.info("Loading Seurat...")
suppressPackageStartupMessages(library(Seurat))
library(sttkit)

singlecell <- opt$singlecell
if (grepl("list$",opt$singlecell)) {
    singlecell <- cli_check_file_list(opt$singlecell)
}
if (!is.null(opt$labels_singlecell)) {
    labels <- cli_check_file_list(opt$labels_singlecell, check_exists = FALSE)
    if (length(singlecell) != length(labels)) {
        stop("--labels_singlecell requires the same length as --singlecell.")
    }
} else {
    labels <- paste0("reference_", seq(length(singlecell)))
}    

flog.info("Reading --infile (%s)...",
    basename(opt$infile))
infile <- readRDS(opt$infile)

filename_predictions <- sttkit:::.get_serialize_path(opt$outprefix, "_transfer_predictions.rds")
if (!opt$force && file.exists(filename_predictions)) {
    flog.warn("%s exists. Skipping finding transfer predictions. Use --force to overwrite.", filename_predictions)
    prediction.assay <- readRDS(filename_predictions)
} else {
    flog.info("Reading --singlecell (%s)...",
        basename(opt$singlecell))
    singlecell <- lapply(singlecell, readRDS)

    singlecell <- lapply(singlecell, function(x) {
        if ("SCT" %in% Assays(x)) return(x)
        flog.info("Running sctransform --singlecell...")
        SCTransform(x, ncells = opt$num_integration_features, verbose = FALSE)
    })

    singlecell <- lapply(singlecell, function(x) {
        if ("umap" %in% Reductions(x) &&
            "pca" %in% Reductions(x)) return(x)
        flog.info("Running PCA and UMAP on --singlecell...")
        RunPCA(x, verbose = TRUE) %>% RunUMAP(dims = 1:30)
    })
    singlecell <- lapply(singlecell, FindNeighbors, verbose = FALSE)

    flog.info("Calculating transfer anchors...")
    anchors <- lapply(singlecell, function(x)
        FindTransferAnchors(reference = x, query = infile,
            normalization.method = "SCT"))
    flog.info("Calculating transfer predictions....")
    prediction.assay <- lapply(seq_along(anchors), function(i)
        TransferData(anchorset = anchors[[i]], refdata = singlecell[[i]][[opt$refdata]][,1],
            prediction.assay = TRUE, weight.reduction = infile[["pca"]]))
    flog.info("Writing R data structure to %s...", filename_predictions)
    saveRDS(prediction.assay, filename_predictions)
}

.plot_he <- function(x, i) {
    x$predictions <- prediction.assay[[i]]
    DefaultAssay(x) <- "predictions"
    features <- rownames(x)
    ratio <- sttkit:::.get_image_ratio(length(features))
    label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
    flog.info("Generating output plots for %s ...", label)
    filename <- sttkit:::.get_sub_path(opt$outprefix, "he", 
        paste0("_he_labels", label, ".pdf"))
    sttkit:::.plot_spatial_with_image(filename, x, features, width = 10,
        ratio = ratio, png = opt$png, pt.size.factor = opt$dot_size)
}
for (i in seq_along(singlecell)) {
    .plot_he(infile, i)
}
