cl <- NULL
if ("--mpi" %in% commandArgs()) {
    library(doMPI)
    cl <- startMPIcluster()
    registerDoMPI(cl)
    library(NMF)
#    library(future)
#    plan(cluster, workers = cl)
}
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help="Infile RDS file from st_normalize.R."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--num_clusters"), action = "store", type = "integer", 
        default = formals(sttkit::cluster_bayesspace)$num_clusters, 
        help="Number of clusters. Default tries to find a reasonble number [default %default]"),
    make_option(c("--test_num_clusters"), action = "store", type = "character", 
        default = paste(default = range(eval(formals(sttkit::cluster_bayesspace)$test_num_clusters)), collapse=":"),
        help="List of number of clusters to check [default %default]"),
    make_option(c("--num_iter"), action = "store", type = "integer", 
        default = 100000,
        help="Number of optimization iterations [default %default]"),
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
   
tmp <- strsplit(opt$test_num_clusters, ":")[[1]]
opt$test_num_clusters <- seq(tmp[1], tmp[2])
    
flog.info("Loading Seurat...")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sttkit))
suppressPackageStartupMessages(library(BayesSpace))
log_file <- paste0(opt$outprefix, "_enhance.log")
if (!is.null(log_file)) flog.appender(appender.tee(log_file))

single_input <- TRUE
gmt <- NULL
extra_gmt <- NULL
infile <- opt$infile
labels <- NULL
if (grepl("list$",opt$infile)) {
    stop("st_enhance.R currently does not support lists of multiple input files. Provide a single RDS file")
}    

filename <- sttkit:::.get_serialize_path(opt$outprefix, "_singlecellexperiment.rds")
if (!opt$force && file.exists(filename)) {
    flog.warn("%s exists. Skipping clustering. Use --force to overwrite.", filename)
    ndata <- readRDS(filename)
} else {
    flog.info("Loading infile %s...", paste(sapply(infile, basename), collapse=", "))
    ndata <- readRDS(opt$infile)
    log.normalize <- "SCT" %in% Seurat::Assays(ndata)

    if (!log.normalize) { 
        flog.info("Skipping log normalization in BayesSpace, using provided SCT logcounts instead.")
    }    
    ndata <- as_SingleCellExperiment(ndata, log.normalize = log.normalize)
    ndata <- cluster_bayesspace(ndata,
        test_num_clusters = opt$test_num_clusters,
        num_clusters = opt$num_clusters
    )
    flog.info("Writing R data structure to %s...", filename)
    saveRDS(ndata, file = filename)
}    

filename <- sttkit:::.get_serialize_path(opt$outprefix, "_bayesspace_enhanced.rds")
if (!opt$force && file.exists(filename)) {
    flog.warn("%s exists. Skipping enhancement. Use --force to overwrite.", filename)
    ndata_enhanced <- readRDS(filename)
} else {
    if (is.null(opt$num_clusters)) {
        opt$num_clusters <- attr(ndata, "q.auto.selected")
    }    

    if (is.null(opt$num_clusters)) {
        opt$num_clusters <- length(unique(colData(ndata)$spatial.cluster))
    } 
    flog.info("Enhancing %i clusters. This will take a while...", opt$num_clusters)
    ndata_enhanced <- spatialEnhance(ndata,
        q = opt$num_clusters,
        nrep = opt$num_iter,
        verbose = TRUE)
    flog.info("Writing R data structure to %s...", filename)
    saveRDS(ndata_enhanced, file = filename)
}

