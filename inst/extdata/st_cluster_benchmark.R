suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infiles"), action = "store", type = "character", default = NULL,
        help="List of infile RDS from st_cluster.R."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--hejpeg"), action = "store", type = "character", default = NULL,
        help="Optional path to a JPEG containing cropped HE image."),
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
ndata <- lapply(infiles, readRDS)
libs <- lapply(ndata, function(x) x$library[1])

isc <- Reduce(intersect, lapply(ndata, colnames))
ndata <- lapply(ndata, function(x) x[,isc])

filename <- paste0(opt$outprefix, "_cluster_benchmark_heatmap.pdf")
pdf(filename)
for (i in seq_along(ndata[-1])) {
    for (j in seq(i+1, length(ndata))) {
        ndata[[i]]$idents1 <- as.character(Idents(ndata[[i]]))
        ndata[[i]]$idents2 <- as.character(Idents(ndata[[j]]))

        tbl <- table(ndata[[i]]$idents1, ndata[[i]]$idents2) 
        flog.info("chisq statistic %.3f.", chisq.test(tbl)$statistic)
        if (require(pheatmap)) {
            pheatmap(tbl)
        } else {
            heatmap(tbl)
        }        
        if (!is.null(opt$hejpeg)) {
            .plot_he_cluster(ndata, i, j, opt$outprefix, opt$hejpeg)
        }
    }    
}
dev.off()

if (!is.null(opt$hejpeg)) {
    for (i in seq_along(ndata[-1])) {
        for (j in seq(i+1, length(ndata))) {
            m <- sttkit:::.get_sample_consistency_matrix(ndata[[i]], "idents1", "idents2")
            sttkit:::.plot_consistency_matrix(m, ndata[[i]], ndata[[i]], ndata[[j]], 
                hejpeg = opt$hejpeg, prefix = opt$outprefix)
        }
    }    

    for (i in seq_along(ndata)) {
        if ("scran.cluster" %in% colnames(ndata[[i]]@meta.data)) {
            m <- sttkit:::.get_sample_consistency_matrix(ndata[[i]], "idents1", "scran.cluster")
            sttkit:::.plot_consistency_matrix(m, ndata[[i]], ndata[[i]], ndata[[j]], 
                hejpeg = opt$hejpeg, prefix = opt$outprefix, suffix = "_scran")
        }
    }    
}

