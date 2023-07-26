suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(data.table))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help = "Infile Seurat RDS of st_normalize. Should be unscaled unless normalization is sctransform"),
    make_option(c("--infile_raw"), action = "store", type = "character", default = NULL,
        help = "Infile Seurat RDS of st_normalize. Expects the not normalized RDS file containing raw counts."),
    make_option(c("--singlecell"), action = "store", type = "character", default = NULL,
        help = "Path to a RDS file containing a (list of) Seurat single cell object(s) for spot deconvolution."),
    make_option(c("--labels_singlecell"), action = "store", type = "character", default = NULL,
        help = "Optional list of labels --singlecell"),
    make_option(c("--refdata"), action = "store", type = "character", default = "type",
        help = "Meta data column with prediction labels in --singlecell"),
    make_option(c("--integration_method"), action = "store", type = "character", default = "seurat",
        help = "Integration: Choose between 'seurat' (default), 'celltrek', 'rctd' (alias for 'rctd_multi'), 'rctd_full', 'giotto', 'scvi_destvi', 'scvi_cell2location' as integration method"),
    make_option(c("--downsample_cells"), action = "store", type = "integer", default = 3000,
        help = "Integration: Use that many random cells per refdata [default %default]"),
    make_option(c("--condition"), action = "store", type = "character", default = NULL,
        help = "Optional meta data column with conditions in --singlecell. If provided, will split by condition."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help = "Outfile."),
    make_option(c("--num_integration_features"), action = "store", type = "integer", default = 3000,
        help = "Integration: Use that many features [default %default]"),
    make_option(c("--resolution"), action = "store", type = "double",
        default = 0.8,
        help = "Resolution values for clustering [default %default]"),
    make_option(c("--nmf_ident"), action = "store", type = "integer",
        default = NULL,
        help = "Set Idents(infile) to NMF of specified rank [default %default]"),
    make_option(c("--markers"), action = "store_true",
        default = FALSE,
        help = "Find markers for --singlecell clusters."),
    make_option(c("--keep_mito"), action = "store_true",
        default = FALSE,
        help = "By default, remove mitochondrial genes from --singlecell."),
    make_option(c("--keep_ribo"), action = "store_true",
        default = FALSE,
        help = "By default, remove ribosomal genes from --singlecell."),
    make_option(c("--num_cells_per_spot"), action = "store", type = "integer", default = 20,
        help = "Deconvolution: Average number of cells per spot [default %default]"),
    make_option(c("--simulation"), action = "store_true",
        default = FALSE,
        help = "Subcluster the reference cells, specified by the call attribute [default %default]"),
    make_option(c("--gmt"), action = "store", type = "character", default = NULL,
        help = "Input GMT file(s) for celltrek visualization. Signature names must pattern match cell types in --refdata."),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.6,
        help = "Size of dots on H&E [default %default]"),
    make_option(c("--infer_cna"), action = "store_true", default = FALSE,
        help = "Try to infer copy number alterations to label tumor clusters (works with 'seurat' integration only)."),
    make_option(c("--cna_cutoff"), action = "store", type = "double", default = 0.3,
        help = "Default infercnv cutoff, requires --infer_cna [default %default]"),
    make_option(c("--cna_hmm"), action = "store_true", default = FALSE,
        help = "Run infercnv HMM."),
    make_option(c("--cna_output_normal_counts"), action = "store_true", default = FALSE,
        help = "Output counts of normal clusters, useful for creating a normal reference database."),
    make_option(c("--cna_scale_data"), action = "store_true", default = FALSE,
        help = "Scale data in inferCNV, useful when --cna_pool_of_normals contains very different data."),
    make_option(c("--cna_pool_of_normals"), action = "store", type = "character", default = NULL,
        help = "Matrix or (matrices in .list file) used for infercnv normalization."),
    make_option(c("--cna_idents_ignore"), action = "store", type = "character", default = NULL,
        help = "(low quality) clusters ignored in infercnv [default none]."),
    make_option(c("--cna_idents_reference"), action = "store", type = "character", default = NULL,
        help = "Clusters used for reference in infercnv [default auto]."),
    make_option(c("--cna_subdir"), action = "store", type = "character", default = "",
        help = "Put inferCNV output in this sub-directory [default %default]."),
    make_option(c("--cna_assay"), action = "store", type = "character", default = "Spatial",
        help = "Extract counts from this assay [default %default]."),
    make_option(c("--run_coloc"), action = "store", default = FALSE,
        help = "Run Co-locolization analysis (--integration_method should be 'celltrek')"),
    make_option(c("--run_coexp"), action = "store", default = FALSE,
        help = "Run Co-expression analysis (--integration_method should be 'celltrek')"),
    make_option(c("--coexp_cell_types"), action = "store", type = "character", default = NULL,
        help = "Cell type(s) for co-expression analysis (use with '--run_coexp')"),
    make_option(c("--no_crop"), action = "store_true", default = FALSE,
        help = "Do not crop H&E image."),
    make_option(c("--png"), action = "store_true", default = FALSE,
        help = "Generate PNG version of output plots."),
    make_option(c("--serialize"), action = "store_true", default = FALSE,
        help = "Serialize processed --singlecell object. Can be large."),
    make_option(c("--verbose"), action = "store_true", default = FALSE,
        help = "Verbose output"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE,
        help = "Overwrite existing files")
)

opt <- parse_args(OptionParser(option_list = option_list))
saveRDS(opt, "opt.RDS")

if (is.null(opt[["infile"]])) {
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

if (!opt$integration_method %in% c("seurat", "celltrek", "rctd", "rctd_full", "rctd_multi", "giotto", "scvi", "scvi_destvi", "scvi_cell2location")) {
  stop("Integration: unknown integration method seleted.")
}

flog.info("Loading Seurat...")
suppressPackageStartupMessages(library(Seurat))
library(sttkit)
library(grid)

singlecell <- opt$singlecell
if (grepl("list$", opt$singlecell)) {
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
    basename(opt[["infile"]]))
infile <- readRDS(opt[["infile"]])

if (!is.null(opt$nmf_ident)) {
    library(NMF)
    infile <- set_idents_nmf(infile, k = opt$nmf_ident, stop_if_unavail = TRUE)
}

find_pred <- TRUE
if (opt$integration_method == "rctd") opt$integration_method <- "rctd_multi"
if (opt$integration_method == "scvi") opt$integration_method <- "scvi_destvi"
label_integration_method <- opt$integration_method
tmp <- strsplit(opt$integration_method, "_")[[1]]
opt$integration_method <- tmp[1]
if (length(tmp) > 1) opt$sub_integration_method <- tmp[2]


if (opt$integration_method == "seurat") {
    filename_predictions_old <- sttkit:::.get_serialize_path(opt$outprefix,
        paste0("_", digest(labels), "_transfer_predictions.rds"))
    filename_predictions <- sttkit:::.get_serialize_path(opt$outprefix,
        paste0("_", digest(labels), "_", label_integration_method, "_transfer_predictions.rds"))
#TODO remove
    if (file.exists(filename_predictions_old)) {
          file.copy(filename_predictions_old, filename_predictions)
          file.remove(filename_predictions_old)
    }
    if (!opt$force && file.exists(filename_predictions)) {
        flog.warn("%s exists. Skipping finding transfer predictions. Use --force to overwrite.", filename_predictions)
        prediction.assay <- readRDS(filename_predictions)
        find_pred <- FALSE
    }
} else if (opt$integration_method %in% c("rctd", "giotto", "scvi")) {
    filename_predictions <- sttkit:::.get_serialize_path(opt$outprefix,
        paste0("_", digest(labels), "_", label_integration_method, "_transfer_predictions.rds"))
    if (!opt$force && file.exists(filename_predictions)) {
        flog.warn("%s exists. Skipping %s prediction. Use --force to overwrite.",
            label_integration_method, filename_predictions)
        prediction.assay <- readRDS(filename_predictions)
        find_pred <- FALSE
    }
    if (opt$integration_method == "rctd") suppressPackageStartupMessages(library(spacexr))
    if (opt$integration_method == "giotto") suppressPackageStartupMessages(library(Giotto))
    if (opt$integration_method == "scvi") {
        suppressPackageStartupMessages(library(reticulate))
        suppressPackageStartupMessages(library(sceasy))
        suppressPackageStartupMessages(library(anndata))
        if (dir.exists(Sys.getenv("CONDA_PREFIX"))) {
            reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"))
        }
    }         
} else if (opt$integration_method == "celltrek") {
    filename_traint <- sttkit:::.get_serialize_path(opt$outprefix, "_traint.rds")
    filename_celltrek <- sttkit:::.get_serialize_path(opt$outprefix, "_celltrek.rds")
    if (!require("CellTrek", quietly = TRUE)) {
         stop("This function requires the CellTrek package.")
    }
    if (!opt$force && file.exists(filename_celltrek) && file.exists(filename_traint)) {
        flog.warn("%s and %s exist. Skipping running celltrek. Use --force to overwrite.", filename_celltrek, filename_traint)
        train <- readRDS(filename_traint)
        celltrek_predictions <- readRDS(filename_celltrek)
        find_pred <- FALSE
    }
}

filename_giotto_results <- sttkit:::.get_serialize_path(opt$outprefix,
    paste0("_", digest(labels), "_giotto_results.rds"))
filename_giotto_matrix <- sttkit:::.get_serialize_path(opt$outprefix,
    paste0("_", digest(labels), "_giotto_sign_matrix.rds"))

sign_matrix <- NULL
if (file.exists(filename_giotto_matrix)) {
    sign_matrix <- readRDS(filename_giotto_matrix)
}
    
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
        if (grepl(".rds$", tolower(x))) readRDS(x)
        else if (grepl("h5ad$", tolower(x))) ReadH5AD(x)
    }))
    if (!is.null(opt$downsample_cells)) {
        flog.info("Downsampling --singlecell to %i cells per %s annotation",
            opt$downsample_cells, opt$refdata)
        singlecell <- lapply(singlecell, function(x) {
            if (max(table(x[[opt$refdata]]), na.rm = TRUE) <= opt$downsample_cells * 1.1) return(x)
            x[, unlist(lapply(split(Cells(x), x[[opt$refdata]]), function(y)
                sample(y, min(length(y), opt$downsample_cells), replace = FALSE)))]
        })
    }
    if (!opt[["keep_mito"]]) {
        flog.info("Removing mitochondrial genes from --singlecell.")
        singlecell <- lapply(singlecell, function(x) {
            mito_feats <- grep(pattern = regex_mito(), x = rownames(x = x), value = TRUE)
            x[!rownames(x) %in% mito_feats,]
        })
    }    
    if (!opt[["keep_ribo"]]) {
        flog.info("Removing ribosomal genes from --singlecell.")
        singlecell <- lapply(singlecell, function(x) {
            ribo_feats <- grep(pattern = regex_ribo(), x = rownames(x = x), value = TRUE)
            x[!rownames(x) %in% ribo_feats,]
        })
    }
    # remove cells from too rare cell-types
    singlecell <- lapply(singlecell, function(x) {
        idx <- table(x[[opt$refdata]][,1])[x[[opt$refdata]][,1]] >= 5
        if (any(!idx)) { 
            flog.warn("Removing cells from rare cell-types (< 5 cells).")
            x <- x[, idx]
        }
        return(x)
    })
    if (!is.null(opt$condition)) {
        flog.info("Splitting --singlecell according --condition %s", opt$condition)
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
    method_wants_sct <- c("seurat", "celltrek")

    flog.info("Making sure --singlecell is compatible with current Seurat version...")
    singlecell <- lapply(singlecell, function(x) {
        x <- UpdateSeuratObject(x)
        if ("SCT" %in% Assays(x) || !opt$integration_method %in% method_wants_sct) return(x)
        flog.info("Running sctransform --singlecell...")
        vst.flavor <- !is.null(infile@commands$SCTransform.Spatial) &&
            any(grepl("vst.flavor", infile@commands$SCTransform.Spatial@call.string))
        if (vst.flavor) {
            flog.info("I think v2 of sctransform was used. Try normalizing the single cell data manually if it fails.")
            SCTransform(x, ncells = opt$num_integration_features, vst.flavor = "v2", verbose = FALSE)
        } else {
            SCTransform(x, ncells = opt$num_integration_features, verbose = FALSE)
        }    
    })
    method_wants_clustering <- c("seurat", "celltrek")
    if (opt$integration_method %in% method_wants_clustering) {
        singlecell <- lapply(singlecell, function(x) {
            DefaultAssay(x) <- "SCT"
            if ("umap" %in% Reductions(x) && "pca" %in% Reductions(x))
                if (DefaultAssay(x[["pca"]]) == "SCT" && DefaultAssay(x[["umap"]]) == "SCT")
                    return(x)
            flog.info("Running PCA and UMAP on --singlecell...")
            RunPCA(x, verbose = TRUE) %>% RunUMAP(dims = 1:30)
        })
        flog.info("Clustering --singlecell...")
        singlecell <- lapply(singlecell, FindNeighbors, verbose = FALSE)
        singlecell <- lapply(singlecell, FindClusters, resolution = opt$resolution, verbose = FALSE)
    }
    if (opt$serialize) {
        flog.info("Writing R data structure to %s...", filename_singlecell)
        names(singlecell) <- labels
        sttkit:::.serialize(singlecell, opt$outprefix, "_singlecell.rds")
    }
    flog.info("Done preprocessing --singlecell. If you plan to use this reference for multiple samples, use --serialize to skip those steps for future samples.")
}

if (find_pred == TRUE) {
    if (opt$integration_method == "seurat") {
        flog.info("Calculating transfer anchors...")
        anchors <- lapply(singlecell, function(x)
            FindTransferAnchors(reference = x, query = infile,
                normalization.method = "SCT", dims = 1:30))
        flog.info("Calculating transfer predictions...")
        prediction.assay <- lapply(seq_along(anchors), function(i)
            TransferData(anchorset = anchors[[i]], refdata = singlecell[[i]][[opt$refdata]][, 1],
                prediction.assay = TRUE, weight.reduction = infile[["pca"]], dims = 1:30))
        flog.info("Writing R data structure to %s...", filename_predictions)
        saveRDS(prediction.assay, filename_predictions)

    } else if (opt$integration_method == "rctd") {
        flog.info("Converting --singlecell to spacexr format...")
        singlecell_rctd <- lapply(singlecell, as_Reference, refdata = opt$refdata, require_int = FALSE)
        flog.info("Converting --infile to spacexr format...")
        infile_rctd <- as_SpatialRNA(infile)
        flog.info("Preparing RCTD...")
        myRCTDs <- lapply(singlecell_rctd, function(sc) spacexr::create.RCTD(infile_rctd, sc, max_cores = 1, CELL_MIN_INSTANCE = 3))
        flog.info("Running RCTD...")
        doublet_mode <- ifelse(is.null(opt$sub_integration_method), "multi", opt$sub_integration_method)
        flog.info("Using doublet mode %s.", doublet_mode)
        myRCTDs <- lapply(myRCTDs, spacexr::run.RCTD, doublet_mode = doublet_mode)
        rctd_results <- lapply(myRCTDs, function(x) x@results)
        singlecell_rctd <- NULL
        infile_rctd <- NULL
        filename_rctd_results <- sttkit:::.get_serialize_path(opt$outprefix,
            paste0("_", digest(labels), "_rctd_results.rds"))
        flog.info("Writing R data structure to %s...", filename_rctd_results)
        saveRDS(rctd_results, filename_rctd_results)
        prediction.assay <- lapply(myRCTDs, as_AssayObject)
        myRCTDs <- NULL
        flog.info("Writing R data structure to %s...", filename_predictions)
        saveRDS(prediction.assay, filename_predictions)
    } else if (opt$integration_method == "giotto") {
        instrs <- createGiottoInstructions()

        flog.info("Converting --singlecell to giotto format...")
        singlecell_giotto <- lapply(singlecell, function(x) {
           tmp <- x@meta.data
           cell_metadata <- data.table(cell_ID = rownames(tmp), tmp)

            createGiottoObject(
                expression = GetAssayData(x, assay = "RNA", slot = "counts"),
                cell_metadata = list(cell = list(rna = cell_metadata)),
                instructions = instrs
            )
        })

        flog.info("Normalizing --singlecell with giotto...")
        singlecell_giotto <- lapply(singlecell_giotto, function(gobject) {
            gobject <- normalizeGiotto(gobject = gobject)
            gobject <- calculateHVF(gobject = gobject, show_plot = FALSE, return_plot = FALSE)
            return(gobject)
        })
        flog.info("Clustering --singlecell with giotto...")
        singlecell_giotto <- lapply(singlecell_giotto, function(gobject) {
            gene_metadata <- fDataDT(gobject)
            featgenes <- gene_metadata[hvf == 'yes']$feat_ID
            gobject <- Giotto::runPCA(gobject = gobject, feats_to_use = featgenes, scale_unit = FALSE)
            signPCA(gobject, feats_to_use = featgenes, scale_unit = FALSE, show_plot = FALSE, return_plot = FALSE)
            return(gobject)
        })
        flog.info("Finding markers for --singlecell with giotto...")
        sign_matrix <- find_giotto_dwls_matrix(singlecell_giotto, opt$refdata)
        flog.info("Writing R data structure to %s...", filename_giotto_matrix)
        saveRDS(sign_matrix, filename_giotto_matrix)

        flog.info("Converting --infile to giotto format...")
        infile_giotto <- as_GiottoObject(infile, instructions = instrs)
        flog.info("Normalizing --infile with giotto...")
        infile_giotto <- normalizeGiotto(gobject = infile_giotto)
        infile_giotto <- calculateHVF(gobject = infile_giotto, show_plot = FALSE, return_plot = FALSE)
        flog.info("Clustering --infile with giotto...")
        gene_metadata <- fDataDT(infile_giotto)
        featgenes <- gene_metadata[hvf == 'yes']$feat_ID
        infile_giotto <- runPCA(gobject = infile_giotto, feats_to_use = featgenes, scale_unit = FALSE)
        signPCA(infile_giotto, feats_to_use = featgenes, scale_unit = FALSE, show_plot = FALSE, return_plot = FALSE)
        infile_giotto <- runUMAP(infile_giotto, dimensions_to_use = 1:10)
        infile_giotto <- createNearestNetwork(gobject = infile_giotto, dimensions_to_use = 1:10, k = 15)
        infile_giotto <- doLeidenCluster(gobject = infile_giotto, resolution = 0.4, n_iterations = 1000)
        flog.info("Running DWLS. Might take a while...")
        for (i in seq_along(singlecell_giotto)) {
            infile_giotto <- runDWLSDeconv(infile_giotto, sign_matrix = sign_matrix[[i]]$matrix[sign_matrix[[i]]$sig_feats, ],
                n_cell = opt$num_cells_per_spot, name = paste0("DWLS.", i))
        }

        singlecell_giotto <- NULL
        flog.info("Writing R data structure to %s...", filename_giotto_results)
        saveRDS(infile_giotto, filename_giotto_results)

        prediction.assay <- lapply(infile_giotto@spatial_enrichment$cell$rna, as_AssayObject)
        flog.info("Writing R data structure to %s...", filename_predictions)
        saveRDS(prediction.assay, filename_predictions)
        infile_giotto <- NULL
    } else if (opt$integration_method == "scvi") {
        flog.info("Loading scanpy and scvi python packages...")
        sc <- import("scanpy", convert = FALSE)
        scvi <- import("scvi", convert = FALSE)
        if (opt$sub_integration_method == "destvi") {
            prediction.assay <- lapply(singlecell, function(sc_seurat) {
                feats <- intersect(rownames(sc_seurat), rownames(infile))
                feats <- feats[feats %in% union(VariableFeatures(sc_seurat), VariableFeatures(infile))]
                flog.info("Converting --singlecell and --infile to anndata...")
                sc_adata <- convertFormat(sc_seurat[feats,], from = "seurat",
                    to = "anndata", main_layer = "counts", drop_single_values = FALSE)
                infile_adata <- convertFormat(infile[feats,], from = "seurat", to = "anndata",
                    assay = "Spatial", main_layer = "counts", drop_single_values = FALSE)
                flog.info("Running DestVI scLVM. Will take a while...")
                scvi$model$CondSCVI$setup_anndata(sc_adata, labels_key = opt$refdata)
                sclvm <- scvi$model$CondSCVI(sc_adata, weight_obs = TRUE)
                sclvm$train(max_epochs = as.integer(250))
                scvi$model$DestVI$setup_anndata(infile_adata)
                stlvm <- scvi$model$DestVI$from_rna_model(infile_adata, sclvm)
                flog.info("Running DestVI stLVM. Will take a while...")
                stlvm$train(max_epochs = as.integer(2500))
                infile_adata$obsm["proportions"] <- stlvm$get_proportions()
                as_AssayObject(infile_adata)
            })
        } else if (opt$sub_integration_method == "cell2location") {
            flog.info("Loading cell2location, pandas and numpy python packages...")
            cell2location <- import("cell2location", convert = FALSE)
            pd <- import("pandas", convert = FALSE)
            np <- import("numpy", convert = FALSE)
            prediction.assay <- lapply(singlecell, function(sc_seurat) {
                sc_seurat <- SetAssayData(sc_seurat, slot = "counts",
                    new.data = round(GetAssayData(sc_seurat, slot = "counts")))
                feats <- intersect(rownames(sc_seurat), rownames(infile))
                feats <- feats[feats %in% union(VariableFeatures(sc_seurat), VariableFeatures(infile))]
                flog.info("Converting --singlecell and --infile to anndata...")
                sc_adata <- convertFormat(sc_seurat[feats,], from = "seurat",
                    to = "anndata", main_layer = "counts", drop_single_values = FALSE)
                infile_adata <- convertFormat(infile[feats,], from = "seurat", to = "anndata",
                    assay = "Spatial", main_layer = "counts", drop_single_values = FALSE)
                sc$pp$filter_genes(sc_adata,min_cells = 1)
                sc$pp$filter_cells(sc_adata,min_genes = 1)
                selected <- cell2location$utils$filter_genes(
                    sc_adata, cell_count_cutoff = 5, cell_percentage_cutoff2 = 0.03, nonz_mean_cutoff = 1.12)

                sc_adata <- r_to_py(py_to_r(sc_adata)[, selected])
                cell2location$models$RegressionModel$setup_anndata(
                    adata = sc_adata,
                    labels_key = opt[["refdata"]]
                )
                mod <- cell2location$models$RegressionModel(sc_adata)
                flog.info("Running cell2location on --singlecell. Will take a while...")
                mod$train(max_epochs = as.integer(250), batch_size = as.integer(2500), train_size = 1, lr = 0.002)
                sc_adata <- mod$export_posterior(
                    sc_adata, sample_kwargs = list('num_samples' = as.integer(1000), 'batch_size'= as.integer(2500))
                )
                if (!is.null(py_to_r(sc_adata$varm)$means_per_cluster_mu_fg)) {
                    cols <- paste0("means_per_cluster_mu_fg_", py_to_r(sc_adata$uns['mod']['factor_names']))
                    inf_aver <- py_to_r(sc_adata$varm)[['means_per_cluster_mu_fg']][, cols]
                    colnames(inf_aver) <- py_to_r(sc_adata$uns['mod']['factor_names'])
                    inf_aver <- r_to_py(inf_aver)
                }
                intersect <- np$intersect1d(infile_adata$var_names, inf_aver$index)
                sc_adata <- r_to_py(py_to_r(sc_adata)[, intersect])
                infile_adata <- r_to_py(py_to_r(infile_adata)[, intersect])
                inf_aver <- r_to_py(py_to_r(inf_aver)[py_to_r(intersect), ])
                cell2location$models$Cell2location$setup_anndata(adata = infile_adata)
                flog.info("Running cell2location on --infile. Will take a while...")
                mod <- cell2location$models$Cell2location(
                    infile_adata, cell_state_df = inf_aver,
                    N_cells_per_location = opt[["num_cells_per_spot"]],
                    detection_alpha = 200
                )
                mod$train(
                    max_epochs = as.integer(30000),    
                    batch_size = py_none(),
                    train_size = as.integer(1)
                )    
                flog.info("Running cell2location posterior sampling...")
                infile_adata <- mod$export_posterior(
                    infile_adata, sample_kwargs = list('num_samples' = as.integer(1000), 'batch_size'= mod$adata$n_obs))
                as_AssayObject(infile_adata)
            })
        }    
        flog.info("Writing R data structure to %s...", filename_predictions)
        saveRDS(prediction.assay, filename_predictions)
    } else if (opt$integration_method == "celltrek") {
        infile <- RenameCells(infile, new.names = make.names(Cells(infile)))
        singlecell <- lapply(singlecell, function(sc) RenameCells(sc, new.names = make.names(Cells(sc))))
        train <- lapply(singlecell, function(x)
            CellTrek::traint(st_data = infile, sc_data = x,
                sc_assay = "RNA", cell_names = opt$refdata))
        flog.info("Writing R data structure to %s ...", filename_traint)
        saveRDS(train, filename_traint)

        celltrek_predictions <- lapply(seq_along(singlecell), function(i) {
            ctm <- CellTrek::celltrek(st_sc_int = train[[i]], int_assay = "traint",
                     sc_data = singlecell[[i]], sc_assay = "RNA",
                            reduction = "pca", intp = TRUE, intp_pnt = 5000,
                            intp_lin = FALSE, nPCs = 30, ntree = 1000,
                            dist_thresh = 0.55, top_spot = 5, spot_n = 5, repel_r = 20,
                            repel_iter = 20, keep_model = TRUE)
            cidx <- which(colnames(ctm$celltrek@meta.data) == opt$refdata)
            ctm$celltrek$cell_type <- factor(ctm$celltrek@meta.data[, cidx],
                levels = sort(unique(ctm$celltrek@meta.data[, cidx])))
            ctm$celltrek
        })
        flog.info("Writing R data structures to %s ...", filename_celltrek)
        saveRDS(celltrek_predictions, filename_celltrek)
    }
    if (opt$markers) {
        flog.info("Finding single cell cluster markers...")
        invisible(lapply(seq_along(singlecell), function(i) {
            Idents(singlecell[[i]]) <- singlecell[[i]][[opt$refdata]][, 1]
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

.run_infer <- function(seurat_obj, seurat_raw_obj = seurat_obj,
                       ref_group_names = NULL, out_dir = tempfile(),
                       feature = "ident", HMM = FALSE, cutoff = 0.4,
                       scale_data = opt$cna_scale_data) {

    if (!file.exists(out_dir)) {
        out_dir <- gsub(" ", "_", out_dir)
        dir.create(out_dir)
    }
    count_matrix <- GetAssayData(seurat_raw_obj, slot = "counts",
        assay = opt$cna_assay)[, colnames(seurat_obj)]
    annotations <- FetchData(seurat_obj, feature)
    types <- unique(annotations[, 1])
    if (is.null(ref_group_names)) {
        ref_group_names <- types[!types %in% c("Tumor", "?", "Unassigned")]
    }
    if (!length(ref_group_names) || is.na(ref_group_names)) ref_group_names <- NULL
    if (!is.null(opt$cna_output_normal_counts) &&
        !is.null(ref_group_names)) {
        filename <- file.path(out_dir, "normal_counts.tsv.gz")
        flog.info("Writing normal count data to %s...", basename(filename))
        idx <- annotations[, 1] %in% ref_group_names
        data.table::fwrite(data.table::as.data.table(count_matrix[, idx]), file = filename,
            sep = "\t", quote = FALSE)
    }
    pon_count_matrix <- .get_infer_cna_pon()
    if (!is.null(pon_count_matrix)) {
        isc <- intersect(rownames(count_matrix), rownames(pon_count_matrix))
        count_matrix <- cbind(count_matrix[isc, ], pon_count_matrix[isc, ])
        pon_annotations <- data.frame(
            ident = rep("PoN", ncol(pon_count_matrix)),
            row.names = colnames(pon_count_matrix)
        )
        colnames(pon_annotations) <- colnames(annotations)
        annotations[, 1] <- as.character(annotations[, 1])
        annotations <- rbind(annotations, pon_annotations)
        annotations[, 1] <- as.factor(annotations[, 1])
        ref_group_names <- c(ref_group_names, "PoN")
    }
    message("Refgroups: ", paste(ref_group_names, collapse = ","))
    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = count_matrix,
        gene_order_file = infercnv_genes_example,
        annotations_file = annotations,
        ref_group_names = ref_group_names
    )
    infercnv_obj <- infercnv::run(infercnv_obj,
                                   cutoff = cutoff,
                                   out_dir = out_dir,
                                   cluster_by_groups = TRUE,
                                   denoise = TRUE,
                                   HMM = HMM,
                                   num_threads = 2,
                                   scale_data = scale_data,
                                   sd_amplifier = 1.0
                                   )
    saveRDS(infercnv_obj, file = file.path(out_dir, "infercnv_obj.rds"))
    infercnv_obj
}

.get_infer_cna_pon <- function() {
    if (is.null(opt$cna_pool_of_normals)) return(NULL)
    .load_infer_cna_pon <- function(file, i) {
        if (tolower(tools::file_ext(file)) == "rds") {
            seurat_obj <- readRDS(file)
            count_matrix <- GetAssayData(seurat_obj, slot = "counts", assay = opt$cna_assay)
        } else {
            count_matrix <- data.table::fread(file)
        }
        colnames(count_matrix) <- gsub("-\\d$",
            paste0("-pon-", i), colnames(count_matrix))
        return(count_matrix)
    }
    if (grepl("list$", opt$cna_pool_of_normals)) {
        files <- cli_check_file_list(opt$cna_pool_of_normals)
    } else {
        files <- opt$cna_pool_of_normals
    }
    m <- lapply(seq_along(files), function(i) .load_infer_cna_pon(files[i], i))
    if (length(m) < 2) return(m[[1]])
    isc <- Reduce(intersect, lapply(m, rownames))
    return(Reduce(cbind, lapply(m, function(x) x[isc, ])))
}

.infer_cna <- function(x,
                       idents_ignore = opt$cna_idents_ignore,
                       idents_ref = opt$cna_idents_reference,
                       tumor_ps = c("Unassigned", "Tumor", "?")) {
     out_dir <- file.path(dirname(opt$outprefix), "infer_cna", opt$cna_subdir)
     if (!file.exists(out_dir)) {
         out_dir <- gsub(" ", "_", out_dir)
         dir.create(out_dir)
     }

     if (is.null(idents_ref)) {
         p <- lapply(prediction.assay, function(pa) {
            x$predictions <- pa
            GetTransferPredictions(x)
         })
         p <- do.call(cbind, p)
         tumor_fraction <- apply(p, 1, function(x) sum(x %in% tumor_ps) / length(x))
         tumor_fraction_mean <- sapply(split(tumor_fraction, Idents(x)), mean)
         tumor_fraction_o <- apply(p, 1, function(x) sum(x %in% "Tumor") / length(x))
         tumor_fraction_o_mean <- sapply(split(tumor_fraction_o, Idents(x)), mean)
         ref_group_names <- names(tumor_fraction_mean[tumor_fraction_mean < median(tumor_fraction_mean)])
         if (max(tumor_fraction_o_mean) > 0) {
            obs_group_names <- names(tumor_fraction_o_mean[tumor_fraction_o_mean > median(tumor_fraction_o_mean)])
            ref_group_names <- ref_group_names[!ref_group_names %in% obs_group_names]
        }
    } else {
        ref_group_names <- idents_ref[idents_ref %in% Idents(x)]
        if (length(ref_group_names)) {
            filename <- file.path(out_dir, "normal_counts.rds")
            xx <- x[, Idents(x) %in% ref_group_names]
            if (ncol(xx)) {
                flog.info("Writing normal data to %s...", basename(filename))
                saveRDS(xx, filename)
            }
        }
    }
    if (!is.null(idents_ignore)) {
        x <- x[, !Idents(x) %in% strsplit(idents_ignore, ":")[[1]]]
    }
    if (!is.null(opt[["infile_raw"]])) {
        flog.info("Reading --infile_raw (%s)...",
            basename(opt[["infile_raw"]]))
        seurat_raw_obj <- readRDS(opt[["infile_raw"]])
    } else {
        seurat_raw_obj <- x
    }
    .run_infer(x, seurat_raw_obj = seurat_raw_obj, ref_group_names = ref_group_names, out_dir = out_dir, cutoff = opt$cna_cutoff, HMM = opt$cna_hmm)
}

.plot_signature <- function(ndata, prefix, label, gmt, cells = NULL) {
    # make sure he subdir exists
    filename <- sttkit:::.get_sub_path(opt$outprefix, "he", "tmp.pdf")

    filename <- sttkit:::.get_sub_path(prefix, "he", paste0("_he_celltrek_signature_scores", label, ".pdf"))
    ndata <- plot_signatures(ndata, file = filename, gmt = gmt, 
        cells = cells, pt.size.factor = opt$dot_size, png = opt$png)

    ndata_rna <- ndata
    DefaultAssay(ndata_rna) <- names(ndata@assays)[1]
    features <- sort(unique(unlist(gmt)))
    features <- features[features %in% rownames(ndata_rna)]

    feature.matrix <- GetAssayData(ndata, slot = "data")
    cells <- which(colSums(feature.matrix[which(rownames(feature.matrix) %in% features),])>0)

    #feature.matrix[apply(feature.matrix, 2, function(x) {x==0})] <- NA
    #ndata.trimmed <- SetAssayData(ndata, slot = "data", new.data = feature.matrix)

    #my.pal <- brewer.pal(9, "Reds")
    flog.info("Plotting single feature counts...")
    plot_features(filename, object = ndata_rna, 
          features = features, 
          prefix = prefix, 
          suffix = paste0("_he_celltrek_features", label, ".pdf"),
          subdir = "he",
          cells = cells, pt.size.factor = opt$dot_size,
          image.alpha = 0, png = opt$png,
          crop = !opt$no_crop,
          plot_correlations = TRUE)

    #flog.info("Plotting scaled single feature counts...", filename)
    #plot_features(filename, object = ndata,
    #      features = features, 
    #      prefix = prefix, 
    #      suffix = paste0("_he_celltrek_features_scaled", label, ".pdf"),
    #      subdir = "he",
    #      cells = cells, pt.size.factor = opt$dot_size, 
    #      image.alpha = 0, cols = my.pal, png = opt$png,
    #      plot_correlations = TRUE)
}

.plot_he_ct <- function(x, y, i) {
    cidx <- which(colnames(x[[i]]@meta.data) == opt$refdata)
    features <- unique(x[[i]]@meta.data[, cidx])
    features <- features[!is.na(features)]
    label <- if (is.null(labels[i])) "" else paste0("_", labels[i])

    if(!is.null(opt$gmt)) {
      gmt <- read_signatures(opt$gmt, y[[i]])

      ctypes <- unique(gsub(" ", "", y[[i]]@meta.data$cell_type))
      ctypes <- unlist(strsplit(ctypes, split="/"))
      gmt <- gmt[names(gmt) %in% grep(paste(ctypes, collapse="|"), names(gmt), value = TRUE)]
      if(length(gmt) > 0) {
        .plot_signature(y[[i]], prefix = opt$outprefix, label = label, gmt = gmt)
      }
    }

    c2 <- SpatialDimPlot(y[[i]], group.by = "cell_type",
        pt.size.factor = opt$dot_size, cols = rep("red", length(features)))[[1]] +
              facet_wrap(cell_type ~ .) + theme(legend.position = "none")
    filename <- sttkit:::.get_sub_path(opt$outprefix, "he",
            suffix = paste0("_he_celltrek_labels", label, ".pdf"))
    #filename <- file.path(dirname(opt$outprefix), "he", paste0(basename(opt$outprefix), "_he_celltrek_labels", label, ".pdf"))
    pdf(filename, width = 10, height = 10)
    print(c2)
    dev.off()
    if (opt$png) {
        png(gsub(".pdf$", ".png", filename), width = 10, height = 10, units = "in", res = 150)
        print(patchwork::wrap_plots(c2))
        dev.off()
    #    dev.off() # TODO: device flush not happening properly earlier
    }

    library(Polychrome)
    if (length(features <= 26)) {
        umap_pal <- alphabet.colors(length(features))
    } else if (length(features <= 36)) {
        umap_pal <- palette36.colors(length(features))
    } else {
        flog.warn("Too many cell-types to plot UMAP for (%d). Skipping...", length(features))
    }

    names(umap_pal) <- features
    t1 <- DimPlot(x[[i]], label = TRUE, label.size = 4.5, group.by = "type", shuffle=TRUE)
    t2 <- DimPlot(x[[i]][, which(x[[i]]$type == "sc")], label = TRUE, label.size = 4.5,
                  group.by = opt$refdata, cols = umap_pal, shuffle = TRUE,
                  na.value = "white", repel = TRUE)
    filename <- sttkit:::.get_sub_path(opt$outprefix, "umap",
            suffix = paste0("_umap_labels", label, ".pdf"))
    pdf(filename, width = 10, height = 4.5)
    print(t1 + t2 + patchwork::plot_layout(widths = c(1.2, 1.8)))
    dev.off()
    if (opt$png) {
        png(gsub(".pdf$", ".png", filename), width = 10, height = 4.5, units = "in", res = 150)
        print(t1 + t2 + patchwork::plot_layout(widths = c(1.2, 1.8)))
        dev.off()
    }

    filename <- sttkit:::.get_sub_path(opt$outprefix, "umap",
            suffix = paste0("_single_cell_barplot", label, ".pdf"))
    pdf(filename, width = 5, height = 4.5)
    barplot(sort(table(celltrek_predictions[[i]]$cell_type)), col=umap_pal, las=2)
    dev.off()
    if (opt$png) {
        png(gsub(".pdf$", ".png", filename), width = 5, height = 4.5, units = "in", res = 150)
        barplot(sort(table(celltrek_predictions[[i]]$cell_type)), col=umap_pal, las=2)
        dev.off()
    }
}

for (i in seq_along(singlecell)) {
    if (opt$integration_method %in% c('seurat', 'rctd', 'giotto', 'scvi') ) {
        sttkit::plot_predictions(infile, prediction.assay[[i]], 
            label = labels[i], label_integration_method = label_integration_method, 
            prefix = opt[["outprefix"]], png = opt$png, pt.size.factor = opt$dot_size,
            crop = !opt$no_crop)
    } else if (opt$integration_method == 'celltrek') {
        .plot_he_ct(train, celltrek_predictions, i)
    }
}

if (opt$infer_cna) {
    if (opt$integration_method == 'celltrek') {
        flog.warn("Works with 'seurat' integration only, not 'celltrek'), skipping...")
    } else if (opt$integration_method == 'seurat') {
        if (!require("infercnv")) {
            flog.warn("--infercnv requires the infercnv package.")
        } else {
            .infer_cna(infile)
        }
    }
}

if (opt$run_coloc) {
    if (opt$integration_method == "seurat") {
        flog.warn("Co-localization analysis works with 'celltrek' integration only, not 'seurat'), skipping...")
    } else if (opt$integration_method == "celltrek") {
        opath <- file.path(dirname(opt$outprefix), "co-analysis")
        if (!dir.exists(opath)) {
            dir.create(opath)
        }
        for(i in seq_along(celltrek_predictions)) {
            ct <- celltrek_predictions[[i]]

            label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
            fileprefix <- file.path(opath, paste0(basename(opt$outprefix)))

            p <- plot_scoloc(ct)
            pdf(paste0(fileprefix, "_colocalization", label, ".pdf"))
            print(p)
            dev.off()
            if (opt$png) {
                png(paste0(fileprefix, "_colocalization", label, ".png"))
                print(p)
                dev.off()
            }
        }
    }
}

if (opt$run_coexp) {
    if (opt$integration_method == "seurat") {
        flog.warn("Co-expression analysis works with 'celltrek' integration only, not 'seurat'), skipping...")
    } else if (opt$integration_method == "celltrek") {
        for(i in seq_along(celltrek_predictions)) {
            ctp <- celltrek_predictions[[i]]
            if (!is.null(opt$coexp_cell_types)) {
                cell_types <- unlist(strsplit(opt$coexp_cell_types, ","))
            } else {
                cell_types <- levels(ctp@meta.data$cell_type)
            }
            label <- if (is.null(labels[i])) "" else paste0("_", labels[i], "_")
            for (ctype in cell_types) {
                flog.info("Co-expression analysis for %s", ctype)
                if (length(which(ctp$cell_type == ctype)) < 10)
                    next
                ct <- subset(ctp, subset = cell_type==ctype)
                ct@assays$RNA@scale.data <- matrix(NA, 1, 1)
                fileprefix <- file.path(dirname(opt$outprefix), "co-analysis",
                                        paste0(basename(opt$outprefix), label, ctype))

                ct <- FindVariableFeatures(ct)
                vst_df <- ct@assays$RNA@meta.features %>% data.frame %>% mutate(id = rownames(.))
                nz_test <- apply(as.matrix(ct[['RNA']]@data), 1, function(x) mean(x != 0) * 100)
                hz_gene <- names(nz_test)[nz_test < 20]
                mt_gene <- grep('^MT-', rownames(ct), value = TRUE)
                rp_gene <- grep('^RPL|^RPS', rownames(ct), value = TRUE)
                vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>%
                                     arrange(., -vst.variance.standardized)
                if (nrow(vst_df) > 2000) {
                    feature_temp <- vst_df$id[1:2000]
                } else {
                    feature_temp <- vst_df$id
                }

                ct_scoexp_res_cc <- scoexp(celltrek_inp=ct, assay = 'RNA', approach = 'cc',
                                           gene_select = feature_temp, sigm = 140, avg_cor_min = .4,
                                           zero_cutoff = 3, min_gen = 40, max_gen = 400)
                if (length(ct_scoexp_res_cc$gs) == 0)
                    next
                ct <- AddModuleScore(ct, features = ct_scoexp_res_cc$gs, name = 'CC_',
                                     nbin = 10, ctrl = 50, seed = 42)

                ct_k <- data.frame(gene = unlist(ct_scoexp_res_cc$gs))
                ct_k$G <- toupper(substr(rownames(ct_k), 0, 2))

                write.csv(ct_k, paste0(fileprefix, "_cluster_genes.csv"), row.names=FALSE)
                rownames(ct_k) <- ct_k$gene
                ct_k$gene <- NULL

                library(viridis)
                library(gridBase)
                ph <- pheatmap::pheatmap(ct_scoexp_res_cc$wcor[rownames(ct_k), rownames(ct_k)],
                           clustering_method = 'ward.D2', annotation_row = ct_k,
                           show_rownames = FALSE, show_colnames = FALSE,
                           treeheight_row = 10, treeheight_col = 10,
                           annotation_legend = TRUE, fontsize = 8,
                           color = viridis(10), main = paste(ctype, 'co-expression'))
                dp <- DimPlot(ct, group.by = 'seurat_clusters')
                fp <- FeaturePlot(ct, grep('CC_', colnames(ct@meta.data), value = TRUE))
                sp <- SpatialFeaturePlot(ct, grep('CC_', colnames(ct@meta.data), value = TRUE))

                vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
                pdf(paste0(fileprefix, "_coexpression_plots.pdf"), width = 10)
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(3, 4)))
                pushViewport(vplayout(1:2, 1:2)); print(dp, newpage = FALSE); popViewport()
                pushViewport(vplayout(3, 1:2)); print(fp, newpage = FALSE); popViewport()
                pushViewport(vplayout(3, 3:4)); print(sp, newpage = FALSE); popViewport()
                pushViewport(vplayout(1:2, 3:4)); par(fig = gridFIG(), new = TRUE); print(ph); popViewport()
                dev.off()

                if (opt$png) {
                    png(paste0(fileprefix, "_coexpression_plots.png"), height = 720, width=1080)
                    grid.newpage()
                    pushViewport(viewport(layout = grid.layout(3, 4)))
                    pushViewport(vplayout(1:2, 1:2)); print(dp, newpage = FALSE); popViewport()
                    pushViewport(vplayout(3, 1:2)); print(fp, newpage = FALSE); popViewport()
                    pushViewport(vplayout(3, 3:4)); print(sp, newpage = FALSE); popViewport()
                    pushViewport(vplayout(1:2, 3:4)); par(fig = gridFIG(), new = TRUE); print(ph); popViewport()
                    dev.off()
                }
            }
        }
    }
}


if (!opt$integration_method %in% c("celltrek")) {
    if (length(prediction.assay) > 1) {
        common_labels <- Reduce(intersect, lapply(prediction.assay, function(x) rownames(GetAssayData(x))))
        if (length(common_labels)) {
            common_assay <- sttkit::find_assayobject_consensus(prediction.assay)
            sttkit::plot_predictions(infile, predictions = common_assay,
                label = "consensus", label_integration_method = label_integration_method,
                prefix = opt[["outprefix"]], png = opt$png, pt.size.factor = opt$dot_size, crop = !opt$no_crop)
        }
    }

    if (require("Giotto", quietly = TRUE)) {
        filename_interactions <- sttkit:::.get_serialize_path(opt$outprefix,
            paste0("_", digest(labels), "_", label_integration_method, "_celltype_interactions.rds"))
        filename_interactions_feats <- sttkit:::.get_serialize_path(opt$outprefix,
            paste0("_", digest(labels), "_", label_integration_method, "_celltype_interactions_features.rds"))
        if (!opt$force && file.exists(filename_interactions) && file.exists(filename_interactions_feats)) {
            flog.warn("%s and %s exist. Skipping cell-type proximity analysis. Use --force to overwrite.",
                filename_interactions, filename_interactions_feats)
        } else {    
            flog.info("Running Giotto spatial cell-type proximity enrichment analysis. Might take a while...")
            giotto_object <- as_GiottoObject(infile)

            cpes <- lapply(seq_along(prediction.assay), function(i) {
                giotto_enrichment <- as_spatEnrObj(prediction.assay[[i]])
                giotto_object <- set_spatial_enrichment(gobject = giotto_object, spatenrichment = giotto_enrichment)
                giotto_object <- createSpatialNetwork(giotto_object)
                cpe <- cellProximityEnrichmentSpots(giotto_object,
                    spatial_network_name = "Delaunay_network", cells_in_spot = opt$num_cells_per_spot)
                return(cpe)
            })
            names(cpes) <- labels
            flog.info("Writing R data structure to %s...", basename(filename_interactions))
            saveRDS(cpes, filename_interactions)
            if (0) {
            if (!is.null(sign_matrix)) {
                flog.info("Running Giotto spatial cell-type proximity differential expression analysis. Might take a while...")

                icfs <- lapply(seq_along(prediction.assay), function(i) {
                    giotto_enrichment <- as_spatEnrObj(prediction.assay[[i]])
                    giotto_object <- set_spatial_enrichment(gobject = giotto_object, spatenrichment = giotto_enrichment)
                    giotto_object <- createSpatialNetwork(giotto_object)
                    icf <- findICFSpot(giotto_object, ave_celltype_exp = sign_matrix[[i]]$matrix,
                        selected_features = rownames(sign_matrix[[i]]$matrix), feat_type = "rna", spat_unit = "cell")
                    return(icf)
                })
                names(icfs) <- labels
                flog.info("Writing R data structure to %s...", basename(filename_interactions_feats))
                saveRDS(icfs, filename_interactions_feats)
            }
            }
            for (i in seq_along(cpes)) {
                giotto_enrichment <- as_spatEnrObj(prediction.assay[[i]])
                giotto_object <- set_spatial_enrichment(gobject = giotto_object, spatenrichment = giotto_enrichment)
                giotto_object <- createSpatialNetwork(giotto_object)
                filename_heatmap <- sttkit:::.get_sub_path(opt$outprefix, 
                    "advanced", suffix = paste0("_", labels[i], "_",
                        label_integration_method, "_cell_proximity_enrichment.pdf"))
                pdf(filename_heatmap)
                ret <- try(cellProximityHeatmap(giotto_object, cpes[[i]]))
                dev.off()
                if (is(ret, "try_error")) {
                    flog.warn("Could not generate cell proximity heatmap.")
                    file.remove(filename_heatmap)
                } else if (opt$png) {
                    png(gsub(".pdf$", ".png", filename_heatmap), width = 7, height = 7, units = "in", res = 150)
                    ret <- try(cellProximityHeatmap(giotto_object, cpes[[i]]))
                    dev.off()
                }
            }
        }
    } else {
        flog.info("Install the Giotto package for cell-type interaction analyses.")
    }    
}
flog.info("Done.")
