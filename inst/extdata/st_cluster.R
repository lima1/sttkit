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
    make_option(c("--regressout"), action = "store", type = "character", 
        default = NULL, 
        help="Variables to regress out when --infile is a list of unscaled RDS files [default %default]"),
    make_option(c("--resolution"), action = "store", type = "character", 
        default = "0.8", 
        help="Resolution values for clustering. When multiple are provided, the last one is the main one [default %default]"),
    make_option(c("--hejpeg"), action = "store", type = "character", default = NULL,
        help="Optional path to a JPEG containing cropped HE image."),
    make_option(c("--labels"), action = "store", type = "character", default = NULL,
        help="Optional list of labels for multi-sample analyses"),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.5,
        help="Size of dots on H&E."),
    make_option(c("--markergenes"), action="store", type = "double", default = 30, 
        help="Heatmap: Use that many marker genes per cluster [default %default]"),
    make_option(c("--gmt"), action = "store", type = "character", 
        default = NULL, 
        help="GMT file including genes of interest"),
    make_option(c("--extra_gmt"), action = "store", type = "character", 
        default = NULL, 
        help="GMT file including pathways of interest. The genes are not forced to be included in any merged dataset, i.e. might be removed due to low variance."),
    make_option(c("--suffix_extra_gmt"), action = "store", type = "character", 
        default = "_pathway_", 
        help="Filename suffix for enrichment analysis of --extra_gmt."),
    make_option(c("--min_features"), action="store", type = "double", default = 100, 
        help="Integration: Keep spots that detected that many genes or more [default %default]"),
    make_option(c("--min_spots"), action="store", type = "double", default = 200, 
        help="Integration: Merge samples with fewer detected spots [default %default]"),
    make_option(c("--nmf"), action = "store_true", default = FALSE, 
        help="Do additional NMF clustering"),
    make_option(c("--nmf_ranks"), action = "store", type = "character", default = "4:9",
        help="List of ranks (clusters) to test with --nmf"),
    make_option(c("--nmf_nruns"), action = "store", type = "integer", default = 5,
        help="Number of runs with --nmf"),
    make_option(c("--nmf_max_features"), action = "store", type = "integer", default = NULL,
        help="Speedup by using only the top features"),
    make_option(c("--nmf_randomize"), action = "store_true", default = FALSE, 
        help="Do additional NMF clustering on randomized data to find better rank"),
    make_option(c("--nmf_method"), action = "store", type = "character", default = NULL,
        help="If provided, will use a different than default NMF method."),
    make_option(c("--mpi"), action = "store_true", default = FALSE, 
        help="Use doMPI package for parallel NMF."),
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
    
flog.info("Loading Seurat...")
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sttkit))
log_file <- paste0(opt$outprefix, "_cluster.log")
if (!is.null(log_file)) flog.appender(appender.tee(log_file))

.plot_he_cluster <- function(ndata, prefix, hejpeg, num = "") {
    filename <- sttkit:::.get_sub_path(prefix, "snn", paste0("_he_cluster", num, ".pdf"))
    flog.info("Plotting clusters on H&E for %s...", ndata$library[1])
    ndata$Cluster <- as.factor(as.character(Idents(ndata)))
    pdf(filename, width = 4, height = 3.9)
    plot_features(ndata, features = "Cluster", hejpeg = hejpeg, 
        labels = waiver(), labels_title = "", 
        reorder_clusters = FALSE, size = opt$dot_size, 
        plot_map = TRUE)
    .plot_clustering_overlap(ndata)
    dev.off()
    if (opt$png) {
        filename <- sttkit:::.get_sub_path(prefix, "snn", paste0("_he_cluster", num, ".png"))
        png(filename, width = 4, height = 3.9, units = "in", res = 150)
        plot_features(ndata, features = "Cluster", hejpeg = hejpeg, 
            labels = waiver(), labels_title = "", 
            reorder_clusters = FALSE, size = opt$dot_size, 
            plot_map = FALSE, plot_correlations = FALSE)
        .plot_clustering_overlap(ndata)
        dev.off()
    }
}

.plot_clustering_overlap <- function(x) {
    cluster_cols <- grep("_snn_res", colnames(x@meta.data))
    if (length(cluster_cols) < 2) return()
    cluster_cols <- tail(cluster_cols, 2)
    tbl <- table(x@meta.data[, cluster_cols])
    if (require(pheatmap)) {
        pheatmap(tbl)
    } else {
        heatmap(tbl)
    }      
}
    
.get_serialize_path <- function(prefix, suffix) {
    s_dir <- file.path(dirname(prefix), "serialize")
    file.path(s_dir, paste0(basename(prefix), suffix))
}

.write_clusters <- function(obj, prefix, single_input) {
    sids <- grep("snn_res", colnames(obj@meta.data))
    for (i in sids) {
        filename <- sttkit:::.get_sub_path(prefix, "snn", paste0("_cluster_", 
            colnames(obj@meta.data)[i], ".txt"))
        m <- obj@meta.data[, i, drop = FALSE]
        m[,1] <- as.numeric(m[,1]) + 1
        write.table(m, file = filename, sep = "\t", quote = FALSE, col.names = FALSE)
    } 
    export_snn_loupe(obj, libs = libs, labels = labels, prefix)
}
.plot_cluster_qc <- function(object, prefix) {
	m <- melt(object@meta.data)
	colnames(m)[max(grep("snn_res", colnames(m)))] <- "res"
    filename <- sttkit:::.get_sub_path(prefix, "snn", "_cluster_qc.pdf")
    pdf(filename, width=8, height=4)
	require(ggplot2)
    flog.info("Plotting cluster QC...")
	print(ggplot(m, aes(res, value))+
	      geom_boxplot()+
	      xlab("Cluster Id")+
	      facet_wrap(~variable, scales="free_y"))
	dev.off()
}
.order_features <- function(obj, features) {
    m <- GetAssayData(obj, "scale.data")
    m <- m[rownames(m) %in% features,]
    hc <- hclust(dist(m))
    hc$labels[hc$order]
}
.plot_cluster_heatmaps <- function(obj, prefix, markergenes, single_input, group.by = "ident") {
    filename <- sttkit:::.get_sub_path(prefix, "snn", "_cluster_heatmap.pdf")
    filename_markers <- .get_serialize_path(opt$outprefix, "_snn_markers.rds")
    if (!opt$force && file.exists(filename_markers)) {
        flog.warn("%s exists. Skipping differential expression analysis. Use --force to overwrite.", filename_markers)
        markers <- readRDS(filename_markers)
    } else {
        markers <- FindAllMarkers(obj)
        flog.info("Writing R data structure to %s...", filename_markers)
        sttkit:::.serialize(markers, opt$outprefix, "_snn_markers.rds")
    }
    m <- GetAssayData(obj, slot="scale.data")
    genes <- unique(unlist(lapply(split(markers, markers$cluster), function(x) 
        head(x[which(x$avg_logFC > 0 & x$gene %in% rownames(m)), "gene"], markergenes))))
    genes <- .order_features(obj, genes)
    flog.info("Plotting cluster heatmap...")
    pdf(filename, height = length(genes) / 60 * 6, width = 8)
    print(DoHeatmap(ndata, features = genes, group.by = group.by))
    if (!single_input) {
        print(DoHeatmap(ndata, features = genes, group.by = "library"))
    }
    dev.off()
    flog.info("Plotting PCA heatmap...")
    filename <- paste0(prefix, "_pca_heatmap.pdf")
    pdf(filename, height = 6, width = 8)
    print(DimHeatmap(ndata, dims=1:6, reduction="pca"))
    dev.off()
    filename <- sttkit:::.get_sub_path(prefix, "snn", "_cluster.csv") 
    write.csv(markers, file = filename, row.names = FALSE)
}

single_input <- TRUE
gmt <- NULL
extra_gmt <- NULL
infiles <- opt$infile
labels <- NULL
hejpegs <- NULL
if (grepl("list$",opt$infile)) {
    infiles <- cli_check_file_list(opt$infile)
    if (length(infiles)>1) single_input <- FALSE
    if (!is.null(opt$hejpeg)) {
        hejpegs <- cli_check_file_list(opt$hejpeg)
        if (length(infiles) != length(hejpegs)) {
            stop("--infile as list requires --hejpeg as list.")    
        }
    }
    flog.info("Loading infiles %s...", paste(sapply(infiles, basename), collapse=", "))
    reference_list <- lapply(infiles, readRDS)
    if (!is.null(opt$labels)) {
        labels <- cli_check_file_list(opt$labels, check_exists = FALSE)
        if (length(infiles) != length(labels)) {
            stop("--labels requires the same length as --infile.")    
        }
        reference_list <- lapply(seq_along(reference_list), function(i) {
            reference_list[[i]]$label <- labels[i]
            reference_list[[i]]})
    }    
    if (!is.null(opt$gmt)) {
        ndata_merged <- Reduce(merge, reference_list)
        gmt <- read_signatures(opt$gmt, ndata_merged)
        if (!length(gmt)) {
            stop("No signatures available in --gmt.")
        }    
    }    
    if (!is.null(opt$extra_gmt)) {
        ndata_merged <- Reduce(merge, reference_list)
        extra_gmt <- read_signatures(opt$extra_gmt, ndata_merged)
        if (!length(extra_gmt)) {
            stop("No signatures available in --extra_gmt.")
        }    
    }    
} 


.find_integration_features <- function(reference_list, gmt, prefix) {
    n <-  length(VariableFeatures(reference_list[[1]]))
    flog.info("Selecting %i integration features.", n)
    variable_features <- SelectIntegrationFeatures(reference_list, 
        nfeatures = n)
    if (is.null(gmt)) return(variable_features)
    gp <- plot_gmt_availability(reference_list, gmt)
    filename <- paste0(prefix, "_signatures_availability.pdf")
    pdf(filename, height = 10, width = 10)
    print(gp)
    dev.off()
    # if GMT is provided, add requested genes if present in at least one sample    
    flog.info("Adding additional features provided in --gmt.")
    all_features <- Reduce(intersect, lapply(reference_list, rownames))
    all_non0_features <- Reduce(union, lapply(reference_list, function(x) { 
        m <- GetAssayData(x, "counts")
        rownames(m)[apply(m, 1, max) > 0]
    }))
    wanted_features <- union(variable_features, unlist(gmt))
    wanted_features <- wanted_features[wanted_features %in% all_non0_features]
    features <- wanted_features[wanted_features %in% all_features]
    features
}

.write_labels_diff <- function(ndata, prefix) {
    # TODO remove
    ndata$label <- gsub(" Mouse", "_Mouse", ndata$label)
    Idents(ndata) <- ndata$label

    filename <- sttkit:::.get_sub_path(prefix, "advanced", "_labels_diff.csv") 
    markers <- FindAllMarkers(ndata)
    write.csv(markers, file = filename, row.names = FALSE)
    if (length(grep("_", levels(ndata)))) {
        Idents(ndata) <- gsub("_.*$", "", ndata$label)
        filename <- sttkit:::.get_sub_path(prefix, "advanced", "_labels_diff_2.csv") 
        markers <- FindAllMarkers(ndata)
        write.csv(markers, file = filename, row.names = FALSE)
    }
}    


.parse_regressout <- function() { 
    if (is.null(opt$regressout)) NULL else strsplit(opt$regressout, ":")[[1]]
}

filename <- .get_serialize_path(opt$outprefix, ".rds")
if (!opt$force && file.exists(filename)) {
    flog.warn("%s exists. Skipping clustering. Use --force to overwrite.", filename)
    ndata <- readRDS(filename)
} else {
    if (!single_input) {
        features <- .find_integration_features(reference_list, gmt, prefix = opt$outprefix)
        scale <- if (grepl("unscaled.rds", infiles[1])) TRUE else FALSE
        if (!scale) flog.info("Not scaling input files.")
        ndata <- integrate_spatial(reference = reference_list, features = features,
            reference_technology = "spatial", min_features = opt$min_features, 
            min_spots = opt$min_spots, 
            min_max_counts = 0, scale = scale, force = opt$force,
            prefix = opt$outprefix, verbose = opt$verbose)
        flog.info("Clustering integrated data...")
        ndata <- cluster_integrated(ndata, regressout = .parse_regressout(), 
            force = opt$force, prefix = opt$outprefix, verbose = opt$verbose)
        ndata <- cluster_reference(ndata, resolution = as.numeric(strsplit(opt$resolution, ":")[[1]]),
                             prefix = opt$outprefix, force = opt$force, verbose = opt$verbose)
    } else {    
        ndata <- readRDS(infiles[1])
        ndata <- cluster_spatial(ndata, 
                             resolution = as.numeric(strsplit(opt$resolution, ":")[[1]]))
        if (!is.null(opt$extra_gmt)) {
            extra_gmt <- read_signatures(opt$extra_gmt, ndata)
        }    
    }
    flog.info("Writing R data structure to %s...", filename)
    sttkit:::.serialize(ndata, opt$outprefix, ".rds")
}

if (single_input) {
    libs <- ndata$library[1]
    .plot_he_cluster(ndata, opt$outprefix, opt$hejpeg)
} else {
    libs <- sapply(reference_list, function(x) x$library[1])
    for (i in seq_along(libs)) {
        num <- paste0("_", libs[i])
       .plot_he_cluster(ndata[,ndata$library == libs[i]], opt$outprefix, hejpegs[i],
            num = num)
        #Idents(reference_list[[i]]) <- Idents(ndata[,ndata$library == libs[i]])
        sttkit:::.serialize(ndata[,ndata$library == libs[i]], prefix = opt$outprefix, 
            paste0(num, ".rds"))
    }
    plot_clusters(ndata, opt$outprefix)
    if (!is.null(gmt)) {
        filename <- paste0(opt$outprefix, "_signatures_normalized_counts.pdf")
        flog.info("Plotting normalized signature counts...")
        pdf(filename, width = 10, height = 10 * sttkit:::.get_image_ratio(length(gmt)))
        plot_signatures_fake_bulk(reference_list, plot_pairs = FALSE, 
            plot_bar = TRUE, plot_heatmaps = FALSE, log_trans = FALSE, gmt = gmt)
        dev.off()
        filename <- paste0(opt$outprefix, "_signatures_normalized_counts_heatmaps.pdf")
        pdf(filename, width = 10, height = 10)
        plot_signatures_fake_bulk(reference_list, plot_pairs = FALSE, 
            plot_bar = FALSE, plot_heatmaps = TRUE, log_trans = TRUE, gmt = gmt)
        dev.off()
        gp <- plot_gmt_availability(list(ndata), gmt)
        filename <- paste0(opt$outprefix, "_signatures_availability_after_integration.pdf")
        pdf(filename, height = 10, width = 10)
        print(gp)
        dev.off()
    }
}

.write_gse_results <- function(gse, gse_method, gse_perm_n, ks, 
                               suffix = "_cluster_") {
     eps <- if (gse_method == "perm") 1/gse_perm_n else eps = .Machine$double.eps
     x_df  <- lapply(gse, function(x) do.call(rbind, 
        lapply(seq_along(x), function(i) { 
            idx <- x[[i]]$pvalue < 0.05  
            x[[i]]$fdr <-  p.adjust(x[[i]]$pvalue, method = "BH")
            ret <- data.frame(
                Cluster = i, 
                Signatures = if (any(idx)) rownames(x[[i]])[idx] else "",
                Pvalue = if (any(idx)) format.pval(x[[i]]$pvalue[idx], digits = 2, eps = eps) else "",
                FDR = if (any(idx)) format.pval(x[[i]]$fdr[idx], digits = 2, eps = eps) else "",
                Genes = if (any(idx)) x[[i]]$genes[idx] else "",
                Description = ""
            )
            if (!is.null(x[[i]]$Descr) ) {
                ret$Description <- if (any(idx)) x[[i]]$Descr[idx] else ""
            }
            ret
        })))
     for (i in seq_along(x_df)){
         if (all(x_df[[i]]$Description %in% c("NULL", "", "NA", NA))) {
             x_df[[i]]$Description <- NULL
         }    
         filename <- sttkit:::.get_sub_path(opt$outprefix, "nmf/signatures", 
            paste0("_nmf_", gse_method, suffix, ks[i], ".csv"))
         write.csv(x_df[[i]], file = filename, row.names = FALSE)
     }
}    

if (!is.null(opt$labels)) .write_labels_diff(ndata, opt$outprefix)

.write_clusters(ndata, opt$outprefix, single_input = single_input)
.plot_cluster_qc(ndata, opt$outprefix)
.plot_cluster_heatmaps(ndata, opt$outprefix, opt$markergenes, single_input = single_input)

if (opt$nmf) {

    suppressPackageStartupMessages(library(NMF))
    filename <- .get_serialize_path(opt$outprefix, "_nmf.rds")
    ks <- sort(as.numeric(strsplit(opt$nmf_ranks, ":")[[1]]))
    if (length(ks) > 1) ks <- seq(ks[1], ks[2])

    .run_nmf <- function(randomize = FALSE) {
        r <- if (randomize) "_randomize" else ""
        filesuffix <- paste0("_nmf_", min(ks), "_", max(ks), r, ".rds")
        filename <- .get_serialize_path(opt$outprefix, filesuffix)
        r <- if (randomize) "randomized " else ""
        if (!opt$force && file.exists(filename)) {
            flog.warn("%s exists. Skipping %sNMF clustering. Use --force to overwrite.", 
                filename, r)
            ndata <- readRDS(filename)
        } else {
            flog.info("Running %sNMF clustering with %i runs and %i different ranks. This will probably take a while...", 
                r, opt$nmf_nruns, length(ks)) 
            nmf_method <- if (!is.null(opt$nmf_method)) opt$nmf_method else nmf.getOption('default.algorithm')
            if (is.null(cl)) {
                ndata <- cluster_nmf(ndata, ks, .options = 'v3',
                    nrun = opt$nmf_nruns, randomize = randomize, 
                    max_features = opt$nmf_max_features, method = nmf_method)
            } else {    
                flog.info("Trying to use MPI as requested.")
                ndata <- cluster_nmf(ndata, ks, .options = 'v3P', .pbackend=NULL,
                    nrun = opt$nmf_nruns, randomize = randomize,
                    max_features = opt$nmf_max_features, method = nmf_method)
            }
            flog.info("Writing R data structure to %s...", filename)
            sttkit:::.serialize(ndata, opt$outprefix, filesuffix)
        }
        ndata
    }
    ndata <- .run_nmf()
    if (opt$nmf_randomize) ndata <- .run_nmf(TRUE)
    flog.info("Done with NMF clustering!")
    plot_nmf(ndata, libs, hejpegs, labels, rank = ks, prefix = opt$outprefix, png = opt$png, 
        size = opt$dot_size)
     loupe <- lapply(ks, function(i) 
        export_nmf_loupe(obj = ndata, rank = ks, k = i, libs = libs, 
            labels = labels, prefix = opt$outprefix))
    if (!is.null(opt$gmt)) {
         for (gse_method in c("fisher")) {
             gse_perm_n <- 10000
             gse <- lapply(ks, function(i) calculate_nmf_gse(ndata, gmt,
                rank = ks, k = i, method = gse_method, perm_n = gse_perm_n, 
                alternative = "less", verbose = FALSE))
             # build spreadsheets with NMF/signature associations
             .write_gse_results(gse, gse_method, gse_perm_n, ks)
             idx <- names(gmt) %in% rownames(gse[[1]][[1]])
             for (sig in names(gmt)[idx]) {
                 sig <- sig[sig %in% rownames(gse[[1]][[1]])]
                 plot_nmf_gse(gse, sig = sig, rank = ks,
                    prefix = opt$outprefix, method = gse_method,
                    png = opt$png, subdir = "nmf/signatures/advanced")
             }   
         }
         plot_signatures_nmf(ndata, gmt, gmt_name = gsub(".gmt", "", basename(opt$gmt)),
            rank = ks,
            prefix = opt$outprefix,
            subdir = "nmf/signatures")

    }    
    if (!is.null(opt$extra_gmt)) {
         gse <- lapply(ks, function(i) calculate_nmf_gse(ndata, extra_gmt,
            rank = ks, k = i, verbose = FALSE))
         .write_gse_results(gse, gse_method = "fisher", gse_perm_n, ks,
            suffix = opt$suffix_extra_gmt)
    }    
}

if (length(Images(ndata))) {
    filename_features <- .get_serialize_path(opt$outprefix, "_variable_markvariogram.rds")
    if (!opt$force && file.exists(filename_features)) {
        flog.warn("%s exists. Skipping spatial variation analysis. Use --force to overwrite.", filename_features)
        spatial_features <- readRDS(filename_features)
    } else {
        flog.info("Finding top spatially variable features. This will probably take a while...")
        ndata_split <- SplitObject(ndata, split.by = "library")
        ndata_split <- lapply(seq_along(ndata_split), function(i) 
            FindSpatiallyVariableFeatures(ndata_split[[i]], 
            selection.method = "markvariogram", 
            image = sttkit:::.get_image_slice(ndata_split[[i]]),
            features = head(VariableFeatures(ndata), 1000)))
        spatial_features <- lapply(ndata_split, SpatiallyVariableFeatures, selection.method = "markvariogram")
        libs <- as.character(sapply(ndata_split, function(x) x$library[1]))
        names(spatial_features) <- libs

        flog.info("Writing R data structure to %s...", filename_features)
        sttkit:::.serialize(spatial_features, opt$outprefix, "_variable_markvariogram.rds")
       flog.info("Done with spatial variation analysis!") 

    }
    flog.info("Plotting spatial variation...")
    top_features <- lapply(spatial_features, function(x) head(x, length(x) / 100 * 5))
    ndata_split <- SplitObject(ndata, split.by = "library")
    libs <- as.character(sapply(ndata_split, function(x) x$library[1]))
    for (i in seq_along(ndata_split)) {
        flog.info("Generating output plots for %s ...", libs[i])
        ratio <- sttkit:::.get_image_ratio(length(top_features[i]))
        label <- if (is.null(labels[i])) "" else paste0("_",labels[i])
        libs_label <- if (length(libs) < 2) "" else paste0("_",libs[i])
        filename <- sttkit:::.get_sub_path(opt$outprefix, "spatial_variation", 
            paste0("_he_variable_markvariogram", label, libs_label, ".pdf"))
        pdf(filename, width = 10, height = 10 * ratio)
        ndata_split[[i]] <- plot_features(ndata_split[[i]], features = top_features[[i]], hejpeg = NULL,
        labels = waiver(), labels_title = "", size = opt$dot_size)
        dev.off()
    }    
    
} 

if (opt$mpi) {
    closeCluster(cl)
    mpi.quit()
}
