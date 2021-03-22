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
    make_option(c("--labels"), action = "store", type = "character", default = NULL,
        help="Optional list of labels for multi-sample analyses"),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.6,
        help="Size of dots on H&E."),
    make_option(c("--ncol"), action = "store", type = "double", default = 4,
        help="Number of columns in spatial plot for multiple-sample clustering"),
    make_option(c("--nrow"), action = "store", type = "double", default = 4,
        help="Number of rows in spatial plot for multiple-sample clustering"),
    make_option(c("--markergenes"), action="store", type = "double", default = 30, 
        help="Heatmap: Use that many marker genes per cluster [default %default]"),
    make_option(c("--gmt"), action = "store", type = "character", 
        default = NULL, 
        help="GMT file including genes of interest. Can be a list of files, separated by ':' (e.g. --gmt a.gmt:b.gmt)."),
    make_option(c("--extra_gmt"), action = "store", type = "character", 
        default = NULL, 
        help="GMT file including pathways of interest. The genes are not forced to be included in any merged dataset, i.e. might be removed due to low variance. Can be a list of files, separated by ':' (e.g. --extra_gmt c.gmt:d.gmt)."),
    make_option(c("--suffix_extra_gmt"), action = "store", type = "character", 
        default = "_pathway_", 
        help="Filename suffix for enrichment analysis of --extra_gmt."),
    make_option(c("--min_features"), action="store", type = "double", default = 100, 
        help="Integration: Keep spots that detected that many genes or more [default %default]"),
    make_option(c("--min_spots"), action="store", type = "double", default = 200, 
        help="Integration: Merge samples with fewer detected spots [default %default]"),
    make_option(c("--nmf"), action = "store_true", default = FALSE, 
        help="Do additional NMF clustering"),
    make_option(c("--nmf_ranks"), action = "store", type = "character", default = NULL,
        help="List of ranks (clusters) to test with --nmf. Default will check from min(6, #idents) to 12."),
    make_option(c("--nmf_nruns"), action = "store", type = "integer", default = 5,
        help="Number of runs with --nmf"),
    make_option(c("--nmf_max_features"), action = "store", type = "integer", default = NULL,
        help="Speedup by using only the top features"),
    make_option(c("--nmf_randomize"), action = "store_true", default = FALSE, 
        help="Do additional NMF clustering on randomized data to find better rank"),
    make_option(c("--nmf_method"), action = "store", type = "character", default = NULL,
        help="If provided, will use a different than default NMF method."),
    make_option(c("--nmf_ident"), action = "store", type = "integer", 
        default = NULL, 
        help="Set Idents(infile) to NMF of specified rank after NMF [default %default]"),
    make_option(c("--nmf_cores"), action = "store", type = "integer", 
        default = NULL, 
        help="Number of cores to use when --mpi is NOT used [default %default]"),
    make_option(c("--spatially_variable_method"), action = "store", type = "character", default = "markvariogram:moransi",
        help="Method(s) to find top spatially variable features [default %default]."),
    make_option(c("--spatially_variable_nfeatures"), action = "store", type = "integer", default = 80,
        help="Plot the specified number of top spatially variable features"),
    make_option(c("--nearest_neighbors"), action = "store_true", default = FALSE, 
        help="For multi-sample analyses, visualizes the similarity of slides."),
    make_option(c("--tissue_positions_list"), action = "store", type = "character", default = NULL,
        help = "SpaceRanger 1.1+ tissue positions list, useful in case barcode suffices differ after integretion and fixing Loupe cluster annotation files [default %default]."),
    make_option(c("--cellphonedb"), action = "store_true", default = FALSE, 
        help="Generate output files for cellphonedb."),
    make_option(c("--diff_tests"), action = "store", type = "character", default = "wilcox:negbinom",
        help="Default differential expression method(s) [default %default]."),
    make_option(c("--species"), action = "store", type = "character", default = "Hs",
        help="Currently only used for --cellphonedb to convert symbols to EMSEMBL ids [default %default]."),
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

.plot_he_cluster <- function(ndata, prefix, num = "") {
    #makes sure that snn is available
    filename <- sttkit:::.get_sub_path(prefix, "snn", "tmp.pdf")

    filename <- sttkit:::.get_sub_path(prefix, "snn/he", paste0("_he_cluster", num, ".pdf"))
    flog.info("Plotting clusters on H&E for %s...", ndata$library[1])
    pdf(filename, width = 4, height = 3.9)
    gp <- SpatialDimPlot(ndata, label = TRUE, image = sttkit:::.get_image_slice(ndata), 
        pt.size.factor = opt$dot_size, label.size = 3)
    if (requireNamespace("ggthemes", quietly = TRUE) &&
        length(levels(Idents(ndata))) <= 8) {
        gp <- gp + ggthemes::scale_fill_colorblind()
    }
    print(gp)
    if (length(levels(Idents(ndata))) > 1) .plot_clustering_overlap(ndata)
    invisible(dev.off())
    if (opt$png) {
        filename <- sttkit:::.get_sub_path(prefix, "snn/he", paste0("_he_cluster", num, ".png"))
        png(filename, width = 4, height = 3.9, units = "in", res = 150)
        print(gp)
        invisible(dev.off())
    }
}

.plot_clustering_overlap <- function(x) {
    cluster_cols <- grep("_snn_res", colnames(x@meta.data))
    if (length(cluster_cols) < 2) return()
    cluster_cols <- tail(cluster_cols, 2)
    tbl <- table(x@meta.data[, cluster_cols])
    if (nrow(tbl) < 2) return()
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

.write_clusters <- function(obj, prefix) {
    flog.info("Exporting Loupe files...")
    sids <- grep("snn_res", colnames(obj@meta.data))
    export_snn_loupe(obj, libs = libs, labels = labels, prefix)
}
.plot_cluster_qc <- function(object, prefix) {
	m <- melt(object@meta.data)
	colnames(m)[max(grep("snn_res", colnames(m)))] <- "res"
    filename <- sttkit:::.get_sub_path(prefix, "snn/qc", "_cluster_qc.pdf")
    pdf(filename, width=8, height=4)
	require(ggplot2)
    flog.info("Plotting cluster QC...")
	print(ggplot(m, aes(res, value))+
	      geom_boxplot()+
	      xlab("Cluster Id")+
	      facet_wrap(~variable, scales="free_y"))
    invisible(dev.off())
}
.order_features <- function(obj, features) {
    m <- GetAssayData(obj, "scale.data")
    m <- m[rownames(m) %in% features,]
    if (nrow(m) < 4) return(features)
    hc <- hclust(dist(m))
    hc$labels[hc$order]
}
.plot_cluster_heatmaps <- function(obj, prefix, markergenes, single_input, group.by = "ident") {
    filename <- sttkit:::.get_sub_path(prefix, "snn", "")
    filename <- sttkit:::.get_sub_path(prefix, "snn/heatmap", "_cluster_heatmap.pdf")
    markers <- sttkit:::.find_all_markers(obj, prefix, "_snn_markers.rds")
    m <- GetAssayData(obj, slot="scale.data")
    key <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
    
    genes <- unique(unlist(lapply(split(markers, markers$cluster), function(x) 
        head(x[which(x[[key]] > 0), "gene"], markergenes))))
    genes <- .order_features(obj, genes)
    flog.info("Plotting cluster heatmap...")
    pdf(filename, height = length(genes) / 60 * 6, width = 8)
    print(DoHeatmap(ndata, features = genes, group.by = group.by))
    if (!single_input) {
        print(DoHeatmap(ndata, features = genes, group.by = "library"))
    }
    invisible(dev.off())
    flog.info("Plotting PCA heatmap...")
    filename <- sttkit:::.get_sub_path(prefix, "pca", "_pca_heatmap.pdf") 
    pdf(filename, height = 6, width = 8)
    print(DimHeatmap(ndata, dims=1:6, reduction="pca"))
    invisible(dev.off())
    filename <- sttkit:::.get_sub_path(prefix, "snn/advanced", "_cluster_markers.csv") 
    write.csv(markers, file = filename, row.names = FALSE)
}

single_input <- TRUE
gmt <- NULL
extra_gmt <- NULL
infiles <- opt$infile
labels <- NULL
if (grepl("list$",opt$infile)) {
    infiles <- cli_check_file_list(opt$infile)
    if (length(infiles)>1) single_input <- FALSE
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
    reference_list <- cli_check_lib_ids(reference_list)
    if (!is.null(opt$gmt)) {
        ndata_merged <- Reduce(merge, reference_list)
        gmt <- read_signatures(opt$gmt, ndata_merged)
        if (!length(gmt)) {
            stop("No signatures available in --gmt.")
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
    filename <- sttkit:::.get_sub_path(prefix, "signatures", "_signatures_availability.pdf") 
    pdf(filename, height = 10, width = 10)
    print(gp)
    invisible(dev.off())
    # if GMT is provided, add requested genes if present in at least one sample    
    flog.info("Adding additional features provided in --gmt.")
    all_features <- Reduce(intersect, lapply(reference_list, rownames))
    all_non0_features <- Reduce(union, lapply(reference_list, function(x) { 
        m <- GetAssayData(x, "counts", assay = "Spatial")
        rownames(m)[apply(m, 1, max) > 0]
    }))
    wanted_features <- union(variable_features, unlist(gmt))
    wanted_features <- wanted_features[wanted_features %in% all_non0_features]
    features <- wanted_features[wanted_features %in% all_features]
    features
}

.write_labels_diff <- function(ndata, test.use = "wilcox", prefix) {
    Idents(ndata) <- ndata$label

    filename <- sttkit:::.get_sub_path(prefix, "advanced",
        paste0("_", test.use, "_labels_diff.csv"))
    markers <- sttkit:::.find_all_markers(ndata, prefix,
        paste0("_", test.use, "_snn_labels_markers.rds"),
        test.use = test.use)
    write.csv(markers, file = filename, row.names = FALSE)
    .volcano <- function(markers, filename) {
        if (require(EnhancedVolcano)) {
            if (is.null(markers$avg_log2FC)) {
                markers$avg_log2FC <- log2(exp(markers$avg_logFC))
            }
            pdf(gsub(".csv$", ".pdf", filename)) 
            for (group in levels(markers$cluster)) {
                markers_f <- markers[markers$cluster %in% group,]
                print(EnhancedVolcano(markers_f, lab = markers_f[["gene"]],
                    x = "avg_log2FC", y = "p_val", 
                    title = paste(group, "vs others"), subtitle = ""))
            }
            dev.off()
        }
    }
    .volcano(markers, filename)
    # in case labels contain a suffix, e.g. trt_1, trt_2, ctl_1, ctl_2 
    if (length(grep("_\\d+$", levels(ndata)))) {
        Idents(ndata) <- gsub("_\\d+$", "", ndata$label)
        filename <- sttkit:::.get_sub_path(prefix, "advanced",
            paste0("_", test.use, "_labels_diff_2.csv"))
        markers <- sttkit:::.find_all_markers(ndata, prefix,
            paste0("_", test.use, "_snn_labels_2_markers.rds"),
            test.use = test.use)
        write.csv(markers, file = filename, row.names = FALSE)
        .volcano(markers, filename)
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
    }
    flog.info("Writing R data structure to %s...", filename)
    sttkit:::.serialize(ndata, opt$outprefix, ".rds")
}

# reloading after merging because some genes might drop out
if (!is.null(opt$gmt)) {
    gmt <- read_signatures(opt$gmt, ndata)
}    
if (!is.null(opt$extra_gmt)) {
    extra_gmt <- read_signatures(opt$extra_gmt, ndata)
}    

if (single_input) {
    libs <- ndata$library[1]
    .plot_he_cluster(ndata, opt$outprefix)
    plot_clusters(ndata, opt$outprefix)
} else {
    libs <- sapply(reference_list, function(x) x$library[1])
    for (i in seq_along(libs)) {
        num <- paste0("_", libs[i])
       .plot_he_cluster(ndata[,ndata$library == libs[i]], opt$outprefix, 
            num = num)
        #Idents(reference_list[[i]]) <- Idents(ndata[,ndata$library == libs[i]])
        sttkit:::.serialize(ndata[,ndata$library == libs[i]], prefix = opt$outprefix, 
            paste0(num, ".rds"))
    }
    plot_clusters(ndata, opt$outprefix)
    if (!is.null(labels)) {
        sttkit:::.plot_correlation_labels(ndata, cluster_labels = Idents(ndata), 
                prefix = opt$outprefix, file.path("snn", "advanced"), 
                suffix = paste0("_snn_cluster_label_correlations.pdf")) 
    }
    if (!is.null(gmt)) {
        filename <- sttkit:::.get_sub_path(opt$outprefix, "signatures", "_signatures_normalized_counts.pdf") 
        flog.info("Plotting normalized signature counts...")
        pdf(filename, width = 10, height = 10 * sttkit:::.get_image_ratio(length(gmt)))
        plot_signatures_fake_bulk(reference_list, plot_pairs = FALSE, 
            plot_bar = TRUE, plot_heatmaps = FALSE, log_trans = FALSE, gmt = gmt)
        invisible(dev.off())
        filename <- sttkit:::.get_sub_path(opt$outprefix, "signatures", "_signatures_normalized_counts_heatmaps.pdf") 
        pdf(filename, width = 10, height = 10)
        plot_signatures_fake_bulk(reference_list, plot_pairs = FALSE, 
            plot_bar = FALSE, plot_heatmaps = TRUE, log_trans = TRUE, gmt = gmt)
        invisible(dev.off())
        gp <- plot_gmt_availability(list(ndata), gmt)
        filename <- sttkit:::.get_sub_path(opt$outprefix, "signatures", "_signatures_availability_after_integration.pdf") 
        pdf(filename, height = 10, width = 10)
        print(gp)
        invisible(dev.off())
    }
    if (opt$nearest_neighbors) {
        flog.info("Finding nearest neighbors. This will probably take a while...")
       ndata <- find_nearest_neighbors(ndata)
       gp <-SpatialFeaturePlot(ndata, features = "int.nearest.neighbor",
            combine = FALSE)
       filename <- sttkit:::.get_sub_path(opt$outprefix, "advanced", suffix = paste0("_he_nearest_neighbors.pdf")) 
       ratio <- sttkit:::.get_image_ratio(length(libs))
       pdf(filename, height = 10 * ratio, width = 10)
       print(patchwork::wrap_plots(gp))
       invisible(dev.off())
       if (opt$png) {
           png(gsub(".pdf$", ".png", filename), width = 10, height = 10 * ratio, units = "in", res = 150)
           print(patchwork::wrap_plots(gp))
           invisible(dev.off())
       }    
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
     tmp <- sttkit:::.get_sub_path(opt$outprefix, file.path("nmf", "signatures"),"")
     for (i in seq_along(x_df)){
         if (all(x_df[[i]]$Description %in% c("NULL", "", "NA", NA))) {
             x_df[[i]]$Description <- NULL
         }    
         filename <- sttkit:::.get_sub_path(opt$outprefix, 
            file.path("nmf", "signatures", ks[i]),
            paste0("_nmf_", gse_method, suffix, ks[i], ".csv"))
         write.csv(x_df[[i]], file = filename, row.names = FALSE)
     }
}    

if (!is.null(opt$labels)) {
    for (test.use in strsplit(opt$diff_tests, ":")[[1]]) {
        .write_labels_diff(ndata, prefix = opt$outprefix, test.use = test.use)
    }
}    

.write_clusters(ndata, opt$outprefix)
.plot_cluster_qc(ndata, opt$outprefix)
.plot_cluster_heatmaps(ndata, opt$outprefix, opt$markergenes, single_input = single_input)

if (opt$nmf) {

    suppressPackageStartupMessages(library(NMF))
    filename <- .get_serialize_path(opt$outprefix, "_nmf.rds")
    if (is.null(opt$nmf_ranks)) {
        ks <- seq(min(6, length(levels(Idents(ndata)))), 12)
    } else {    
        ks <- sort(as.numeric(strsplit(opt$nmf_ranks, ":")[[1]]))
    }    
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
                ndata <- cluster_nmf(ndata, ks, .options = paste0('v3', ifelse(is.null(opt$nmf_cores),"", paste0("p", opt$nmf_cores))),
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
    plot_nmf(ndata, libs, labels = labels, rank = ks, prefix = opt$outprefix, png = opt$png, 
        pt.size.factor = opt$dot_size)
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
         plot_signatures_nmf(ndata, gmt, gmt_name = "",
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

if (!is.null(opt$nmf_ident)) {
    if (is.null(opt$nmf)) {
        flog.warn("--nmf_ident requires --nmf")
    }  
    library(NMF)
    old_idents <- Idents(ndata)
    ndata <- set_idents_nmf(ndata, k = opt$nmf_ident, stop_if_unavail = TRUE)
    if (!identical(old_idents, Idents(ndata))) {
        flog.info("Setting idents to NMF %i clustering. This is not serialized.", opt$nmf_ident)
    }
}


if (length(Images(ndata))) {
    methods <- sapply(strsplit(opt$spatially_variable_method, ":")[[1]], trimws)
    for (method in methods) {
        filename_features <- .get_serialize_path(opt$outprefix, paste0("_variable_", method, ".rds"))
        if (!opt$force && file.exists(filename_features)) {
            flog.warn("%s exists. Skipping spatial variation analysis. Use --force to overwrite.", filename_features)
            spatial_features <- readRDS(filename_features)
        } else {
            flog.info("Finding top spatially variable features. This will probably take a while...")
            ndata_split <- SplitObject(ndata, split.by = "library")
            ndata_split <- lapply(seq_along(ndata_split), function(i) { 
                flog.info("Working on %s...", ndata_split[[i]]$library[1])
                FindSpatiallyVariableFeatures(ndata_split[[i]], 
                selection.method = method, 
                image = sttkit:::.get_image_slice(ndata_split[[i]]),
                features = head(VariableFeatures(ndata), 1000))})
            spatial_features <- lapply(ndata_split, SpatiallyVariableFeatures, selection.method = method)
            libs <- as.character(sapply(ndata_split, function(x) x$library[1]))
            names(spatial_features) <- libs

            flog.info("Writing R data structure to %s...", filename_features)
            sttkit:::.serialize(spatial_features, opt$outprefix, paste0("_variable_", method, ".rds"))
            flog.info("Done with spatial variation analysis!") 
        }
        flog.info("Plotting spatial variation...")
        plot_spatially_variable(ndata, labels = labels, method = method, 
            spatial_features = spatial_features, prefix = opt$outprefix,
            number_features = opt$spatially_variable_nfeatures,
            pt.size.factor = opt$dot_size, ncol = opt$ncol, nrow = opt$nrow)
    }
}

if (opt$cellphonedb) {
    orgdb <- paste0("org.", opt$species, ".eg.db")
    assay <- NULL
    slot <- data
    if ("SCT" %in% Assays(ndata)) {
        assay <- "SCT"
        slot <- "counts"
    } else if ("Spatial" %in% Assays(ndata)) {
        assay <- "Spatial"
    }
    if (!require(orgdb, character.only = TRUE)) {
        flog.warn("Install %s to cellphonedb output", orgdb)
    } else {
        cellphone_for_seurat(ndata, get(orgdb), prefix = opt$outprefix,
            slot = slot, assay = assay)
    }
}            
if (opt$mpi) {
    closeCluster(cl)
    mpi.quit()
}
