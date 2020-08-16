suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infiles"), action = "store", type = "character", default = NULL,
        help="List of infile RDS from st_cluster.R."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--nmf_ident"), action = "store", type = "integer", 
        default = NULL, 
        help="Set Idents(infile) to NMF of specified rank after NMF [default %default]"),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.5,
        help="Size of dots on H&E"),
    make_option(c("--verbose"), action = "store_true", default = FALSE, 
        help="Verbose output"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE, 
        help="Overwrite existing files")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$infiles)) {
    stop("Need --infiles")
}
if (is.null(opt$outprefix)) {
    stop("Need --outprefix")
}
if (!dir.exists(dirname(opt$outprefix))) {
    dir.create(dirname(dirname(opt$outprefix)))
    dir.create(dirname(opt$outprefix))
}

.plot_he_cluster <- function(ndata, i, j, prefix, hejpeg) {
    num <- paste0(i, "_", j)
    filename <- paste0(prefix, "_he_cluster_benchmark", num, ".pdf")
    flog.info("Plotting clusters on H&E...")
    pdf(filename, width = 4, height = 3.9)
    ndata[[i]]$Cluster <- as.factor(paste0(
        as.character(Idents(ndata[[i]])), "/",
        as.character(Idents(ndata[[j]]))
    ))
    plot_features(ndata[[i]], features = "Cluster", hejpeg = hejpeg, 
        labels = waiver(), labels_title = "", 
        reorder_clusters = FALSE, size = opt$dot_size, 
        plot_map = TRUE)
    dev.off()
}

    
flog.info("Loading Seurat...")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sttkit))

infiles <- cli_check_file_list(opt$infiles)
flog.info("Loading infiles %s...", paste(sapply(infiles, basename), collapse=", "))
ndata <- lapply(infiles, readRDS)
ndata <- cli_check_lib_ids(ndata)
if (!is.null(opt$nmf_ident)) {
    library(NMF)
    ndata <- lapply(ndata, function(x) {
        old_idents <- Idents(x)
        x <- set_idents_nmf(x, k = opt$nmf_ident)
        if (!identical(old_idents, Idents(x))) {
            flog.info("Setting idents to NMF %i clustering. This is not serialized.", opt$nmf_ident)
        }
        x
    })
}

for (i in seq_along(ndata)) {
    for (j in seq_along(ndata)) {
        if (i <= j) next
        cluster_prediction_strength(ndata[[i]], ndata[[j]], 
            prefix = opt$outprefix)
    }    
}
