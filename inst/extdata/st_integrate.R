suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(digest))

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
    make_option(c("--downsample_cells"), action="store", type = "integer", default = 3000, 
        help = "Integration: Use that many random cells per refdata [default %default]"),
    make_option(c("--condition"), action = "store", type = "character", default = NULL,
        help = "Optional meta data column with conditions in --singlecell. If provided, will split by condition."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help = "Outfile."),
    make_option(c("--num_integration_features"), action="store", type = "integer", default = 3000, 
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
    make_option(c("--simulation"), action = "store_true", 
        default = FALSE, 
        help = "Subcluster the reference cells, specified by the call attribute [default %default]"),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.6,
        help = "Size of dots on H&E [default %default]"),
    make_option(c("--infer_cna"), action = "store_true", default = FALSE, 
        help = "Try to infer copy number alterations to label tumor clusters."),
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
    make_option(c("--png"), action = "store_true", default = FALSE, 
        help = "Generate PNG version of output plots."),
    make_option(c("--serialize"), action = "store_true", default = FALSE, 
        help = "Serialize processed --singlecell object. Can be large."),
    make_option(c("--verbose"), action = "store_true", default = FALSE, 
        help = "Verbose output"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE, 
        help = "Overwrite existing files")
)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

opt <- parse_args(OptionParser(option_list=option_list))

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
    basename(opt[["infile"]]))
infile <- readRDS(opt[["infile"]])

if (!is.null(opt$nmf_ident)) {
    library(NMF)
    infile <- set_idents_nmf(infile, k = opt$nmf_ident, stop_if_unavail = TRUE)
}

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
        if (!is.null(opt$downsample_cells)) {
            flog.info("Downsampling --singlecell to %i cells per %s annotation",
                opt$downsample_cells, opt$refdata)
            singlecell <- lapply(singlecell, function(x) 
                x[, unlist(lapply(split(Cells(x), x[[opt$refdata]]), function(y) 
                    sample(y,min(length(y), opt$downsample_cells), replace = FALSE)))])
        }
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
            x <- UpdateSeuratObject(x)
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
            normalization.method = "SCT", dims = 1:30))
    flog.info("Calculating transfer predictions....")
    prediction.assay <- lapply(seq_along(anchors), function(i)
        TransferData(anchorset = anchors[[i]], refdata = singlecell[[i]][[opt$refdata]][,1],
            prediction.assay = TRUE, weight.reduction = infile[["pca"]], dims = 1:30))
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

.run_infer <- function(seurat_obj, seurat_raw_obj = seurat_obj, 
                       ref_group_names = NULL, out_dir = tempfile(),
                       feature = "ident", HMM=FALSE, cutoff = 0.4,
                       scale_data = opt$cna_scale_data) {

    if (!file.exists(out_dir)) {
        out_dir <- gsub(" ", "_", out_dir)
        dir.create(out_dir)
    }
    count_matrix <- GetAssayData(seurat_raw_obj, slot= "counts",
        assay = opt$cna_assay)[, colnames(seurat_obj)]
    annotations <- FetchData(seurat_obj, feature)
    types <- unique(annotations[,1])
    if (is.null(ref_group_names)) {
        ref_group_names <- types[!types %in% c("Tumor", "?", "Unassigned")]
    }
    if (!length(ref_group_names) || is.na(ref_group_names)) ref_group_names <- NULL
    if (!is.null(opt$cna_output_normal_counts) &&
        !is.null(ref_group_names)) {
        filename <- file.path(out_dir, "normal_counts.tsv.gz")
        flog.info("Writing normal count data to %s...", basename(filename))
        idx <- annotations[,1] %in% ref_group_names
        data.table::fwrite(data.table::as.data.table(count_matrix[,idx]), file = filename,
            sep = "\t", quote = FALSE)
    }    
    pon_count_matrix <- .get_infer_cna_pon()
    if (!is.null(pon_count_matrix)) {
        isc <- intersect(rownames(count_matrix), rownames(pon_count_matrix))
        count_matrix <- cbind(count_matrix[isc,], pon_count_matrix[isc,])
        pon_annotations <- data.frame(
            ident = rep("PoN", ncol(pon_count_matrix)),
            row.names = colnames(pon_count_matrix)
        )
        colnames(pon_annotations) <- colnames(annotations)
        annotations[,1] <- as.character(annotations[,1])
        annotations <- rbind(annotations, pon_annotations)
        annotations[,1] <- as.factor(annotations[,1])
        ref_group_names <- c(ref_group_names, "PoN")
    }    
    message("Refgroups: ", paste(ref_group_names, collapse=","))
    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = count_matrix,
        gene_order_file=genes,
        annotations_file = annotations,
        ref_group_names = ref_group_names
    )
    infercnv_obj <- infercnv::run(infercnv_obj,
                                   cutoff=cutoff,
                                   out_dir=out_dir,
                                   cluster_by_groups=TRUE,
                                   denoise=TRUE,
                                   HMM=HMM,
                                   num_threads=2,
                                   scale_data = scale_data,
                                   sd_amplifier=1.0
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
    if (grepl("list$",opt$cna_pool_of_normals)) {
        files <- cli_check_file_list(opt$cna_pool_of_normals)
    } else {
        files <- opt$cna_pool_of_normals
    }
    m <- lapply(seq_along(files), function(i) .load_infer_cna_pon(files[i], i))
    if (length(m) < 2) return(m[[1]])
    isc <- Reduce(intersect, lapply(m, rownames))
    return(Reduce(cbind, lapply(m, function(x) x[isc,])))    
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
         p <- do.call(cbind,p)
         tumor_fraction <- apply(p,1, function(x) sum(x %in% tumor_ps)/length(x))
         tumor_fraction_mean <- sapply(split(tumor_fraction, Idents(x)), mean)
         tumor_fraction_o <- apply(p,1, function(x) sum(x %in% "Tumor")/length(x))
         tumor_fraction_o_mean <- sapply(split(tumor_fraction_o, Idents(x)), mean)
         ref_group_names <- names(tumor_fraction_mean[tumor_fraction_mean < median(tumor_fraction_mean)])
         if (max(tumor_fraction_o_mean)>0) {
            obs_group_names <- names(tumor_fraction_o_mean[tumor_fraction_o_mean > median(tumor_fraction_o_mean)])
            ref_group_names <- ref_group_names[!ref_group_names %in% obs_group_names]
        }
    } else {
        ref_group_names <- idents_ref[idents_ref %in% Idents(x)]
        if (length(ref_group_names)) {
            filename <- file.path(out_dir, "normal_counts.rds")
            xx <- x[,Idents(x) %in% ref_group_names]
            if (ncol(xx)) {
                flog.info("Writing normal data to %s...", basename(filename))
                saveRDS(xx, filename)
            }
        }
    }    
     if (!is.null(idents_ignore)) {
         x <- x[,!Idents(x) %in% strsplit(idents_ignore, ":")[[1]]]

     }     
    if (!is.null(opt[["infile_raw"]])) {
        flog.info("Reading --infile_raw (%s)...",
            basename(opt[["infile_raw"]]))
        seurat_raw_obj <- readRDS(opt[["infile_raw"]])
    } else {
        seurat_raw_obj <- x
    }        
    .run_infer(x, seurat_raw_obj = seurat_raw_obj, ref_group_names = ref_group_names, out_dir = out_dir, cutoff = opt$cna_cutoff, HMM = opt$cna_hmm )
}

.plot_he <- function(x, i) {
    x$predictions <- prediction.assay[[i]]
    DefaultAssay(x) <- "predictions"
    Idents(x) <- GetTransferPredictions(x)
    features <- names(Matrix::rowSums(GetAssayData(x)) > 0)
    label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
    flog.info("Generating output plots for %s ...", label)
    if (length(Images(x)) > 1 && "library" %in% colnames(x@meta.data)) {
        x_split <- SplitObject(x, split.by = "library")
        libs <- sapply(x_split, function(y) y$library[1])
        libs_label <- rep("", length(libs)) 
        field <- "library"
        if ("label" %in% colnames(x@meta.data)) {
            libs_label <- paste0("_", sapply(x_split, function(y) y$label[1]))
            field <- "label"
        }
        for (j in seq_along(libs)) {
            plot_features(object = x_split[[j]], features = features,
                prefix = opt$outprefix, subdir = "he",
                suffix = paste0("_he_labels", label, "_", libs[j], libs_label[j],".pdf"),
                png = opt$png, pt.size.factor = opt$dot_size)
            filename <- sttkit:::.get_sub_path(opt$outprefix, "he", 
                    suffix = paste0("_he_labels_call", label, "_", libs[j], libs_label[j],".pdf"))
            gp <- SpatialDimPlot(x_split[[j]], label = TRUE, 
                image = sttkit:::.get_image_slice(x_split[[j]]), 
                pt.size.factor = opt$dot_size, label.size = 3)
            pdf(filename, width = 4, height = 3.9)
            print(gp)
            invisible(dev.off())
            if (opt$png) {
                png(gsub("pdf$", "png", filename), width = 4,
                    height = 3.9, units = "in", res = 150)
                print(gp)
                invisible(dev.off())
            }
        }
        filename <- sttkit:::.get_sub_path(opt$outprefix, "advanced", 
                suffix = paste0("_labels", label, "_", libs[j], libs_label[j],".pdf"))
        ratio <- sttkit:::.get_image_ratio(min(6,length(features)))
        glist <- VlnPlot(x, features = features, group.by = field, pt.size = 0.25, combine = FALSE)
        glist <- lapply(glist, function(p) ggplotGrob( p + theme(legend.position = "none") ))
        if (length(features) > 6) {
            glist <- gridExtra::marrangeGrob(glist, ncol = 3, nrow = 2)
        } else {
            glist <- patchwork::wrap_plots(glist)
        }    
        ggsave(filename, glist,
               width = 10, height = 10 * ratio)
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

if (opt$infer_cna) {
    if (!require("infercnv")) {
        flog.warn("--infercnv requires the infercnv package.")
    } else {
        .infer_cna(infile)
    }
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
