suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help="Infile RDS file from st_normalize.R."),
    make_option(c("--gmt"), action = "store", type = "character", default = NULL, 
        help="Input GMT file. Can be a list of files, separated by ':' (e.g. --gmt a.gmt:b.gmt)."),
    make_option(c("--name"), action = "store", type = "character", default = NULL, 
        help="Signature name used in output files"),
    make_option(c("--labels"), action = "store", type = "character", default = NULL,
        help="Optional list of labels for multi-sample analyses"),
    make_option(c("--method"), action = "store", type = "character", default = "seurat",
        help="Method to summarize gene signature scores"),
    make_option(c("--seurat_nbin"), action = "store", type = "integer", default = 24,
        help="Value of AddModuleScore nbin argument [default %default]"),
    make_option(c("--seurat_ctrl"), action = "store", type = "integer", default = 30,
        help="Value of AddModuleScore ctrl argument [default %default]"),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--scale"), action = "store_true", default = FALSE, 
        help="Scale input data"),
    make_option(c("--serialize"), action = "store_true", default = FALSE, 
        help="Serialize infile RDS objects with scores added. By default, only meta data is stored as CSV file"),
    make_option(c("--max_percent_mito"), action = "store", type = "double", default = 15,
        help="Remove low quality spots [default %default]"),
    make_option(c("--single_features"), action = "store_true", default = FALSE, 
        help="Plot all genes in the signatures individually"),
    make_option(c("--palette"), action = "store", type = "character", default = "viridis",
        help="The color palette for dots on H&E plots. Can be either viridis, brewer_single_hue_red/blue/green/grey/orange/purple or color1:color2 [default %default]"),
    make_option(c("--palette_inverse"), action = "store_true", default = FALSE, 
        help="Flip the low:high colors in --palette"),
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
if (is.null(opt$gmt)) {
    stop("Need --gmt")
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

log_file <- paste0(opt$outprefix, "_score.log")
if (!is.null(log_file)) flog.appender(appender.tee(log_file))

single_input <- TRUE
labels <- NULL
if (grepl("list$",opt$infile)) {
    infiles <- cli_check_file_list(opt$infile)
    if (length(infiles)>1) single_input <- FALSE
}
if (!is.null(opt$labels)) {
    labels <- cli_check_file_list(opt$labels, check_exists = FALSE)
    if (length(infiles) != length(labels)) {
        stop("--labels requires the same length as --infile.")
    }
}

.plot_signature <- function(ndata, prefix, gmt, name = NULL, num = "", cells = NULL) {
    name_no_dash <- sub("^_", "", name)
    # make sure he subdir exists
    filename <- sttkit:::.get_sub_path(opt$outprefix, "he", "tmp.pdf")

    filename <- sttkit:::.get_sub_path(prefix, file.path("he", name_no_dash), 
                    paste0("_signature_scores", name, num, ".pdf"))
    ndata <- sttkit:::.filter_object(ndata, opt$max_percent_mito, 0, 0)
    ndata <- plot_signatures(ndata, file = filename, gmt = gmt, 
        cells = cells, pt.size.factor = opt$dot_size,
        nbin = opt$seurat_nbin, ctrl = opt$seurat_ctrl,
        method = opt$method, png = opt$png)
    if (opt$single_features) {
        ndata_rna <- ndata
        DefaultAssay(ndata_rna) <- names(ndata@assays)[1]
        features <- unique(unlist(gmt))
        features <- features[features %in% rownames(ndata_rna)]

        flog.info("Plotting single feature counts...")
        plot_features(filename, object = ndata_rna, 
              features = features, 
              prefix = prefix, 
              suffix = paste0("_feature_counts", name, num, ".pdf"),
              subdir = file.path("he", name_no_dash),
              cells = cells, pt.size.factor = opt$dot_size,
              image.alpha = 0,
              png = opt$png, plot_correlations = TRUE)

        flog.info("Plotting scaled single feature counts...", filename)
        plot_features(filename, object = ndata, 
              features = features, 
              prefix = prefix, 
              suffix = paste0("_feature_scaled", name, num, ".pdf"),
              subdir = file.path("he", name_no_dash),
              cells = cells, pt.size.factor = opt$dot_size, 
              image.alpha = 0,
              png = opt$png, plot_correlations = TRUE)
    }    
    ndata
}
name <- if (is.null(opt$name)) "" else paste0("_", opt$name)
name_no_dash <- sub("^_", "", name)
group.term <- "library"

if (single_input) {
    ndata <- readRDS(opt$infile)
    gmt <- read_signatures(opt$gmt, ndata)
    ndata_merged <- .plot_signature(ndata, opt$outprefix, gmt, name)
    sig_names <- sttkit:::.get_signature_names(ndata_merged, gmt)
} else {
    flog.info("Loading infiles %s...", paste(sapply(infiles, basename), collapse=", "))
    ndata <- lapply(infiles, readRDS)
    if (!is.null(labels)) {
        ndata <- lapply(seq_along(ndata), function(i) {
            ndata[[i]]$label <- labels[i]
            ndata[[i]]})
        group.term <- "label"
    }
    ndata_merged <- sttkit:::.merge_safely(ndata)
    if (opt$scale) {
        ndata_merged <- ScaleData(ndata_merged)
    }
    gmt <- read_signatures(opt$gmt, ndata_merged)
    
    for (i in seq_along(infiles)) {
       if (!is.null(ndata[[i]]@meta.data$label)) {
           num <- paste0("_", ndata[[i]]$label[1], "_", ndata[[i]]$library[1])
       } else {
           num <- paste0("_", ndata[[i]]$library[1])
       }         
       ndata_merged <- .plot_signature(ndata_merged, opt$outprefix,  
            cells = colnames(ndata_merged[, ndata_merged$library==ndata[[i]]$library[1]]),
            gmt = gmt, name = name, num = num)
    }
    sig_names <- sttkit:::.get_signature_names(ndata_merged, gmt)

    flog.info("Plotting signature/library associations...")
    filename <- sttkit:::.get_advanced_path(opt$outprefix, paste0("_signature_scores_library", name, ".pdf"))
    pdf(filename, width = 12, height= 12 * sttkit:::.get_image_ratio(length(sig_names)))
    print(plot_violin(ndata_merged, features = sig_names, group.by = group.term, pt_size = 0))
    print(DotPlot(ndata_merged, features = sig_names, group.by = group.term)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
    dev.off()

    if (length(sig_names) > 2) {
        # make sure correlations subdir exists
        filename <- sttkit:::.get_sub_path(opt$outprefix, "correlations", "tmp.pdf")
         
        filename <- sttkit:::.get_sub_path(opt$outprefix, file.path("correlations", name_no_dash), 
                        paste0("_signature_scores_correlations", name, ".pdf"))
        pdf(filename, width = 10, height = 10, onefile = FALSE)
        plot_correlation_heatmap(lapply(seq_along(ndata), function(i) 
            ndata_merged[, ndata_merged$library == ndata[[i]]$library[1]]),
            sig_names, feature_labels = sig_names)
        dev.off()
        filename <- sttkit:::.get_sub_path(opt$outprefix, file.path("correlations", name_no_dash), 
                        paste0("_signature_scores_correlations_nn", name, ".pdf"))
        pdf(filename, width = 10, height = 10, onefile = FALSE)
        plot_correlation_heatmap(lapply(seq_along(ndata), function(i) 
            ndata_merged[, ndata_merged$library == ndata[[i]]$library[1]]),
            sig_names, feature_labels = sig_names, average_nn = TRUE)
        dev.off()
    }
    # make sure fake_bulk subdir exists
    filename <- sttkit:::.get_sub_path(opt$outprefix, "fake_bulk", "tmp.pdf")

    filename <- sttkit:::.get_sub_path(opt$outprefix, file.path("fake_bulk", name_no_dash),
                    paste0("_signature_normalized_counts", name, ".pdf"))
    flog.info("Plotting normalized signature counts...")
    pdf(filename, width = 10, height = 10 * sttkit:::.get_image_ratio(length(sig_names)))
    gg_data <- plot_signatures_fake_bulk(ndata, plot_pairs = FALSE, 
        plot_bar = TRUE, plot_heatmaps = FALSE, log_trans = FALSE, gmt = gmt)
    dev.off()
    filename <- sttkit:::.get_sub_path(opt$outprefix, file.path("fake_bulk", name_no_dash),
                    paste0("_signature_normalized_counts", name, ".csv"))
    write.csv(gg_data$gmt, file = filename, row.names = FALSE)
    if (opt$png) {
        filename <- sttkit:::.get_sub_path(opt$outprefix, file.path("fake_bulk", name_no_dash),
                        paste0("_signature_normalized_counts", name, ".png"))
        png(filename, width = 10, height = 10 * sttkit:::.get_image_ratio(length(sig_names)), 
            units = "in", res = 150)
        gg_data <- plot_signatures_fake_bulk(ndata, plot_pairs = FALSE, 
            plot_bar = TRUE, plot_heatmaps = FALSE, log_trans = FALSE, gmt = gmt)
        dev.off()
    }

    filename <- sttkit:::.get_sub_path(opt$outprefix, file.path("fake_bulk", name_no_dash),
                    paste0("_signature_normalized_counts_heatmaps", name, ".pdf"))
    pdf(filename, width = 10, height = 10)
    gg_data <- plot_signatures_fake_bulk(ndata, plot_pairs = FALSE, 
        plot_bar = FALSE, plot_heatmaps = TRUE, log_trans = TRUE, gmt = gmt)
    dev.off()
}

flog.info("Plotting signature/cluster associations...")
filename <- sttkit:::.get_advanced_path(opt$outprefix, paste0("_signature_scores_cluster", name, ".pdf"))
pdf(filename, width = 12, height= 12 * sttkit:::.get_image_ratio(length(sig_names)))
print(plot_violin(ndata_merged, features = sig_names, group.by = group.term, pt_size = 0))
dev.off()

filename <- sttkit:::.get_advanced_path(opt$outprefix, paste0("_score", name, ".csv"))
write.csv(ndata_merged@meta.data, file = filename)

if (opt$serialize) {
    sttkit:::.serialize(ndata_merged, opt$outprefix, 
        paste0("_score", name, ".rds"))
}
