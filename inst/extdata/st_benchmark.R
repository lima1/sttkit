suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help="Infile RDS file from st_normalize.R."),
    make_option(c("--htseq"), action = "store", type = "character", default = NULL,
        help="Infile htseq bulk from the RnaSeq pipeline."),
    make_option(c("--gmt"), action = "store", type = "character", 
        default = NULL, 
        help="Input GMT file"),
    make_option(c("--gtf"), action = "store", type = "character", 
        default = NULL, 
        help="GTF including gene symbols to include."),
    make_option(c("--name"), action = "store", type = "character", 
        default = NULL, 
        help="Signature name used in output files"),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--dot_size"), action = "store", type = "double", default = 1.5,
        help="Size of dots on H&E"),
    make_option(c("--hide_r2"), action = "store_true", default = FALSE, 
        help="Do not add R2 to plot"),
    make_option(c("--verbose"), action = "store_true", default = FALSE, 
        help="Verbose output"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE, 
        help="Overwrite existing files")
)

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

.check_file_list <- function(file) {
    flog.info("Parsing input list %s...", file)
    files <- read.delim(file, as.is = TRUE, header = FALSE)[,1]
    files <- trimws(files, which = "both")
    num_exists <- sum(file.exists(files), na.rm = TRUE)
    if (num_exists < length(files)) { 
        stop("File not exists in file ", file)
    }
    files
}

single_input <- TRUE
infiles <- opt$infile
htseqs <- opt$htseq
if (grepl("list$",opt$infile)) {
    infiles <- .check_file_list(opt$infile)
    if (length(infiles)>1) single_input <- FALSE
}

if (!is.null(opt$htseq)) {
    if (grepl("list$",opt$htseq)) {
        htseqs <- .check_file_list(opt$htseq)
    } else {
        htseqs <- opt$htseq
    }
} else {
    stop("Need --htseq.")
}

gtf <- NULL
if (!is.null(opt$gtf)) {
    if (!require("rtracklayer", quietly = TRUE)) {
        stop("Install rtracklayer for --gtf.")
    }
    flog.info("Importing %s...", basename(opt$gtf))
    gtf <- rtracklayer::import(opt$gtf)
}
    
.plot_signature <- function(ndata, prefix, gmt, name = NULL, num = "") {
    name <- if (is.null(name)) "" else paste0("_", name)
    filename <- paste0(prefix, "_signature_scores", name, num, ".pdf")
    if (!opt$force && file.exists(filename)) {
        flog.warn("%s exists. Skipping clustering. Use --force to overwrite.", filename)
    } else {
        ndata <- plot_signatures(ndata, file = filename, gmt = gmt, 
            labels = waiver(), labels_title = "", reorder_clusters = FALSE, 
            plot_correlations = FALSE)
    }
    ndata
}

.load_htseq <- function(htseq, gtf = NULL) {
    suppressPackageStartupMessages(library(edgeR))
    bulk_counts <- edgeR::readDGE(opt$htseq, 
        labels = gsub(".gene_counts.cts", "", basename(htseq)))
    if (!is.null(gtf)) {
        bulk_counts <- bulk_counts[make.names(rownames(bulk_counts)) %in% 
            make.names(as.character(gtf$gene_name)),]
    }
    bulk_norm <- cpm(bulk_counts, log = TRUE)
    # change to log from log2 because Seurat uses log
    list(normalized = log(2^bulk_norm), counts = bulk_counts)
}

.plot_spots_bulk <- function(obj, bulk_norm) {
    ncounts <- as.matrix(GetAssayData(obj))
    # make sure we don't miss genes because of make.names
    names(ncounts) <- make.names(names(ncounts))
    rownames(bulk_norm) <- make.names(rownames(bulk_norm))
    isc <- intersect(rownames(ncounts), 
                     rownames(bulk_norm))
    cor_bulk <- cor(ncounts[isc,], bulk_norm[isc,1])
    obj <- AddMetaData(obj, cor_bulk, "Correlation Bulk")
    plot_features(obj, features = "Correlation Bulk",  
        size = opt$dot_size, labels = waiver(), labels_title = "R")
    obj
}

.plot_venn_bulk <- function(obj, bulk_counts, prefix_csv) {
    counts <- as.matrix(GetAssayData(obj, slot.name = "counts"))
    # make sure we don't miss genes because of make.names
    names(counts) <- make.names(names(counts))
    rownames(bulk_counts) <- make.names(rownames(bulk_counts))
    # only check genes that were annotated in bulk
    counts <- counts[rownames(counts) %in% rownames(bulk_counts),]

    if (require(gplots, quietly = TRUE)) {
        ids <- list(Spatial = rownames(counts), 
                    Bulk = names(bulk_counts[bulk_counts[,1]>0,]))
        gplots::venn(ids)
        s1 <- sort(union(ids$Spatial, ids$Bulk))
        if (length(s1)) { 
            write.csv(matrix(s1, ncol = 1), 
                paste0(prefix_csv, "_union.csv"), 
                row.names = FALSE, quote = FALSE)
        }    
        s1 <- sort(setdiff(ids$Spatial, ids$Bulk))
        if (length(s1)) { 
            write.csv(matrix(s1, ncol = 1), 
                paste0(prefix_csv, "_spatial_only.csv"), 
                row.names = FALSE, quote = FALSE)
        }    
        s1 <- sort(setdiff(ids$Bulk, ids$Spatial))
        if (length(s1)) { 
            write.csv(matrix(s1, ncol = 1), 
                paste0(prefix_csv, "_bulk_only.csv"), 
                row.names = FALSE, quote = FALSE)
        }    
    }
}

.plot_fake_bulk <- function(objs, bulk_norm, labels = NULL, hide_r2 = FALSE) {
    idx <- sapply(objs, is, "Seurat")
    objs <- objs[idx]
    if (!is.null(labels)) labels <- labels[idx]

    gg_data <- do.call(rbind, lapply(seq_along(objs), function(i) {
        obj <- objs[[i]]
        counts <- Matrix::rowSums(GetAssayData(obj, slot.name = "counts"))
        # make sure we don't miss genes because of make.names
        names(counts) <- make.names(names(counts))
        rownames(bulk_norm) <- make.names(rownames(bulk_norm))
        counts <- log(2^cpm(counts, log = TRUE))

        isc <- intersect(rownames(counts), 
                         rownames(bulk_norm))
        
        d_f <- data.frame(Symbol = isc,
                   Id = if (is.null(labels)) as.character(obj$library[1]) else labels[[i]], 
                   Spatial = counts[isc, 1],
                   Bulk = bulk_norm[isc, 1])
        flog.info("Sample %s: Pearson %.2f, Spearman %.2f.", 
            paste(unique(as.character(obj$library)), collapse=", "), 
            cor(d_f$Spatial, d_f$Bulk, method = "pearson", use = "complete.obs"),
            cor(d_f$Spatial, d_f$Bulk, method = "spearman", use = "complete.obs"))
        d_f
    }))
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    gp <- ggplot(gg_data, aes(Spatial, Bulk))+
        geom_point(shape = ".")+
        scale_colour_manual(values=cbPalette)+
        xlab("Spatial")+
        ylab("Bulk")+
        facet_wrap(~Id)
    if (!opt$hide_r2) {    
        if (require(ggpmisc)) {
            my_formula <- y ~ x
            gp <- gp + stat_poly_eq(formula = my_formula, aes(label = ..rr.label..), parse = TRUE)
        } else {
            flog.warn("Install ggpmisc package for R2 annotation.")
        }    
    }
    gp
}

ndata <- lapply(infiles, readRDS)
libs <- lapply(ndata, function(x) x$library[1])
bdata <- lapply(htseqs, .load_htseq, gtf)

if (!is.null(opt$gmt)) {
    gmt <- read_signatures(opt$gmt, ndata[[1]])
}

.plot_benchmark <- function(ndata, j, num1, num2) {
    filename <- paste0(opt$outprefix, "_he_benchmark", num1, num2, ".pdf")
    pdf(filename, width = 4, height = 3.6)
    .plot_spots_bulk(ndata, bdata[[j]]$normalized)
    dev.off()
    filename <- paste0(opt$outprefix, "_benchmark", num1, num2, ".pdf")
    pdf(filename, width = 4, height = 3.6)
    print(.plot_fake_bulk(list(ndata), bdata[[j]]$normalized, hide_r2 = opt$hide_r2))
    dev.off()
    filename <- paste0(opt$outprefix, "_benchmark_clusters", num1, num2, ".pdf")
    pdf(filename, width = 7, 
        height = 7 * sttkit:::.get_image_ratio(length(levels(ndata))))
    print(.plot_fake_bulk(SplitObject(ndata, "ident")[levels(ndata)], 
        bdata[[j]]$normalized, labels = levels(ndata), hide_r2 = opt$hide_r2))
    dev.off()
    filename <- paste0(opt$outprefix, "_venn", num1, num2, ".pdf")
    prefix_csv <- paste0(opt$outprefix, "_venn", num1, num2)
    pdf(filename, width = 4, height = 4)
    .plot_venn_bulk(ndata, bdata[[j]]$counts$counts, prefix_csv)
    dev.off()
}
   
for (i in seq_along(ndata)) {
    num1 <- if (single_input) "" else paste0("_", libs[i])
    for (j in seq_along(bdata)) {
        num2 <- if (length(bdata) == 1) "" else paste0("_bulk", j)
        .plot_benchmark(ndata[[i]], j, num1, num2)
        if (!single_input && i == 1) {
            filename <- paste0(opt$outprefix, "_benchmark_merged", num2, ".pdf")
            pdf(filename, width = 4, height = 3.6)
            ndata_merged <- sttkit:::.merge_safely(ndata)
            print(.plot_fake_bulk(list(ndata_merged), bdata[[j]]$normalized, hide_r2 = opt$hide_r2))
            dev.off()
            filename <- paste0(opt$outprefix, "_venn_merged", num2, ".pdf")
            prefix_csv <- paste0(opt$outprefix, "_venn_merged", num2)
            pdf(filename, width = 4, height = 4)
            .plot_venn_bulk(ndata_merged, bdata[[j]]$counts$counts, prefix_csv)
            dev.off()
        }
    }
}
