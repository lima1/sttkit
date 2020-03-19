suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(digest))

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
    make_option(c("--condition"), action = "store", type = "character", default = NULL,
        help="Optional meta data column with conditions in --singlecell. If provided, will split by condition."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--num_integration_features"), action="store", type = "integer", default = 3000, 
        help="Integration: Use that many features [default %default]"),
    make_option(c("--resolution"), action = "store", type = "double", 
        default = 0.8, 
        help="Resolution values for clustering [default %default]"),
    make_option(c("--markers"), action = "store_true", 
        default = FALSE, 
        help="Find markers for --singlecell clusters."),
    make_option(c("--simulation"), action = "store_true", 
        default = FALSE, 
        help="Subcluster the reference cells, specified by the call attribute [default %default]"),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.6,
        help="Size of dots on H&E [default %default]"),
    make_option(c("--png"), action = "store_true", default = FALSE, 
        help="Generate PNG version of output plots."),
    make_option(c("--serialize"), action = "store_true", default = FALSE, 
        help="Serialize processed --singlecell object. Can be large."),
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

filename_predictions_old <- sttkit:::.get_serialize_path(opt$outprefix, "_transfer_predictions.rds")

filename_predictions <- sttkit:::.get_serialize_path(opt$outprefix, paste0("_", digest(labels), "_transfer_predictions.rds"))
#TODO remove
if (file.exists(filename_predictions_old)) {
    file.copy(filename_predictions_old, filename_predictions)
    file.remove(filename_predictions_old)
}    
if (!opt$force && file.exists(filename_predictions)) {
    flog.warn("%s exists. Skipping finding transfer predictions. Use --force to overwrite.", filename_predictions)
    prediction.assay <- readRDS(filename_predictions)
} else {
    filename_singlecell <- sttkit:::.get_serialize_path(opt$outprefix, "_singlecell.rds")
    if (!opt$force && file.exists(filename_singlecell)) {
        flog.warn("%s exists. Skipping normalization and clustering. Use --force to overwrite.",
            filename_singlecell)
        singlecell <- readRDS(filename_singlecell)
        labels <- names(singlecell)
    } else {
        flog.info("Reading --singlecell (%s)...",
            basename(opt$singlecell))

        singlecell <- unlist(lapply(singlecell, function(x) {
            if(grepl(".rds$", tolower(x))) readRDS(x)
            else if(grepl("h5ad$", tolower(x))) ReadH5AD(x)
        }))
        if (!is.null(opt$condition)) {
            singlecell <- lapply(singlecell, SplitObject, opt$condition)
            labels_new <- lapply(seq_along(labels), function(i) 
                paste0(labels[[i]], "_", names(singlecell[[i]])))
            singlecell <- unlist(singlecell)
            labels <- unlist(labels_new)
        }
        if (!is.null(names(singlecell))) {
            flog.warn("--singlecell already contains labels.")
            labels <- names(singlecell)
        } 

        singlecell <- lapply(singlecell, function(x) {
            if ("SCT" %in% Assays(x)) return(x)
            flog.info("Running sctransform --singlecell...")
            SCTransform(x, ncells = opt$num_integration_features, verbose = FALSE)
        })

        singlecell <- lapply(singlecell, function(x) {
            DefaultAssay(x) <- "SCT"
            if ("umap" %in% Reductions(x) && "pca" %in% Reductions(x))
                if( DefaultAssay(x[["pca"]]) == "SCT" && DefaultAssay(x[["umap"]]) == "SCT")
                    return(x)
            flog.info("Running PCA and UMAP on --singlecell...")
            RunPCA(x, verbose = TRUE) %>% RunUMAP(dims = 1:30)
        })
        flog.info("Clustering --singlecell...")
        singlecell <- lapply(singlecell, FindNeighbors, verbose = FALSE)
        singlecell <- lapply(singlecell, FindClusters, resolution = opt$resolution, verbose = FALSE)
        if (opt$serialize) {
            flog.info("Writing R data structure to %s...", filename_singlecell)
            names(singlecell) <- labels
            sttkit:::.serialize(singlecell, opt$outprefix, "_singlecell.rds")
        }
    }
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

    if (opt$markers) {
        flog.info("Finding single cell cluster markers...")
        invisible(lapply(seq_along(singlecell), function(i) {
            Idents(singlecell[[i]]) <- singlecell[[i]][[opt$refdata]][,1]
            features <- VariableFeatures(singlecell[[i]])
            find_markers(singlecell[[i]][features, ],
                label = labels[[i]], prefix = opt$outprefix,
                resolution = opt$resolution, force = opt$force,
                min.pct = 0.5,  min.diff.pct = 0.25,
                logfc.threshold = 0.5, test.use = "roc"
                )
        }))
    }
}
 
.plot_he <- function(x, i) {
    x$predictions <- prediction.assay[[i]]
    DefaultAssay(x) <- "predictions"
    features <- names(Matrix::rowSums(GetAssayData(x)) > 0)
    label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
    flog.info("Generating output plots for %s ...", label)
    if (length(Images(x)) > 1 && "library" %in% colnames(x@meta.data)) {
        x_split <- SplitObject(x, split.by = "library")
        libs <- sapply(x_split, function(y) y$library[1])
        libs_label <- rep("", length(libs)) 
        if ("label" %in% colnames(x@meta.data)) {
            libs_label <- paste0("_", sapply(x_split, function(y) y$label[1]))
        }
        for (j in seq_along(libs)) {
            plot_features(object = x_split[[j]], features = features,
                prefix = opt$outprefix, subdir = "he",
                suffix = paste0("_he_labels", label, "_", libs[j], libs_label[j],".pdf"),
                png = opt$png, pt.size.factor = opt$dot_size)
        }
    } else {    
        plot_features(object = x, features = features,
            prefix = opt$outprefix, subdir = "he",
            suffix = paste0("_he_labels", label, ".pdf"),
            png = opt$png, pt.size.factor = opt$dot_size)
    }
}
for (i in seq_along(singlecell)) {
    .plot_he(infile, i)
}

if (length(prediction.assay) > 1) {
    common_labels <- Reduce(intersect, lapply(prediction.assay, function(x) rownames(GetAssayData(x))))
    if (length(common_labels)) {
        m <- Reduce("+", lapply(prediction.assay, function(x) GetAssayData(x)[common_labels,]))/length(prediction.assay)
        common_assay <- CreateAssayObject(data = m)
        prediction.assay <- c(prediction.assay, common_assay)
        labels <- c(labels, "consensus")
        .plot_he(infile, length(singlecell) + 1)
    }
}
