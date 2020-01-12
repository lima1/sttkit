suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help="Infile Seurat RDS of st_normalize. Should be unscaled unless normalization is sctransform"),
    make_option(c("--singlecell"), action = "store", type = "character", default = NULL,
        help="Path to a RDS file containing a (list of) Seurat single cell object(s) for spot deconvolution."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--regressout"), action = "store", type = "character", 
        default = NULL, 
        help="Variables to regress out [default %default]"),
    make_option(c("--min_features"), action="store", type = "integer", default = 750, 
        help="Integration: Keep cells that detected that many genes or more [default %default]"),
    make_option(c("--num_integration_features"), action="store", type = "integer", default = 3000, 
        help="Integration: Use that many features [default %default]"),
    make_option(c("--resolution"), action = "store", type = "character", 
        default = "0.4:0.8:1.2", 
        help="Custering: Resolution values. When multiple are provided, the first one is the main one [default %default]"),
    make_option(c("--max_cor"), action="store", type = "double", default = 0.9, 
        help="Clustering: Maximum correlation of gene signature scores [default %default]"),
    make_option(c("--min_markers"), action="store", type = "integer", default = 7, 
        help="Clustering: Minimum number of markers per cluster [default %default]"),
    make_option(c("--max_markers"), action="store", type = "integer", default = NULL, 
        help="Clustering: Maximum number of markers per cluster [default %default]"),
    make_option(c("--idents"), action = "store", type = "character", 
        default = NULL, 
        help="Provide a Seurat object as clustering reference [default %default]"),
    make_option(c("--subclustering"), action = "store", type = "character", 
        default = NULL, 
        help="Subcluster the reference cells, specified by the call attribute [default %default]"),
    make_option(c("--simulation"), action = "store_true", 
        default = FALSE, 
        help="Subcluster the reference cells, specified by the call attribute [default %default]"),
    make_option(c("--hejpeg"), action = "store", type = "character", default = NULL,
        help="Optional path to a JPEG containing cropped HE image."),
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
suppressPackageStartupMessages(library(CellMix))

.get_pair_matrix <- function(z) {
    m <- matrix(0, nrow = nrow(z), ncol = nrow(z))
    for (i in seq(ncol(z))) {
        idx <- which(z[,i] > 0.05)
        if (length(idx)>1) {
            m[idx[1], idx[-1]] <- m[idx[1], idx[-1]] + 1
        } else {
            m[idx[2], idx[1]] <- m[idx[1], idx[1]] + 1
        }
    }
    m[lower.tri(m)] <- NA
    m
}
.get_serialize_path <- function(prefix, suffix) {
    s_dir <- file.path(dirname(prefix), "serialize")
    file.path(s_dir, paste0(basename(prefix), suffix))
}


.filterMarkersIntCor <- function(integrated, sig) {
    markers <- attr(sig, "markers")
    idx <- integrated$technology=="spatial"
    e1 <- ExpressionSet(as.matrix(GetAssayData(integrated[,idx])))
    e2 <- ExpressionSet(as.matrix(GetAssayData(integrated[,!idx])))
    e1 <- e1[rowMax(exprs(e1)) != 0,]
    e2 <- e2[rowMax(exprs(e2)) != 0,]
    merged <- mergeExprs(e1,e2)
    ic  <-metaArray::intcor(merged)
    markers$intcor <- ic$pair.cor[markers$gene,1]
    new.markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(y) y[y$intcor>0.2,]))
    rownames(new.markers) <- new.markers$gene
    new.sig <- sig[unique(new.markers$gene),]
    attr(new.sig, "markers") <- new.markers
    new.sig
}
        
merge_clusters <- function(obj_sc, sc_sig, max_correlation, min_markers, fun_find_markers) {
    old_ncol <- ncol(sc_sig)
    markers_orig <- attr(sc_sig, "markers")
    markers <- split(markers_orig$gene, markers_orig$cluster) 
    while (ncol(sc_sig) > 1) {
        markers <- markers[colnames(sc_sig)]
        lmarkers <- sapply(markers, length)
        # corner case when one cluster has no markers
        if (min(lmarkers) < min_markers) {
            p1 <- which.min(lmarkers)
            cc <- sapply(names(markers), function(i) cor(sc_sig[,names(markers)[p1]], sc_sig[,i]))
            cc[p1] <- -Inf
            p2 <- which.max(cc)
            nn <- max(cc)
        } else {
            flog.info("Calculating module scores for %i markers...", length(markers))

            x <- AddModuleScore(obj_sc, features = markers, ctrl = 30, name = names(markers))
            ss <- x[[paste0(names(markers), seq(length(markers)))]]
            cc <- cor(ss)
            nn <- sapply(seq(nrow(cc)), function(i) max(cc[-i, i]))
            p1 <- which.max(nn)
            p2 <- which(cc[, p1] == nn[p1])
        }
        if (max(nn) > max_correlation || 
                min(length(markers[[p1]]), length(markers[[p2]])) < min_markers ) {
           flog.info("Merging clusters %s and %s (%.2f correlation, sizes %i and %i).", 
                names(markers)[p1], names(markers)[p2], max(nn), 
                length(markers[[p1]]), length(markers[[p2]])) 
           new_ids <- colnames(sc_sig)
           # merge markers
           markers[[p1]] <- unlist(markers[c(p1,p2)])
           markers[[p2]] <- markers[[p1]]

           new_pick <- ifelse( 
            table(Idents(obj_sc))[names(markers)[p1]] >  
            table(Idents(obj_sc))[names(markers)[p2]], 1, 2) 
           idx <- new_ids %in% names(markers)[c(p1,p2)]
           new_ids[idx] <- new_ids[idx][new_pick] 
           names(x = new_ids) <- levels(x = obj_sc)
           obj_sc <- RenameIdents(object = obj_sc, new_ids)
           flog.info("%i idents.", length(levels(obj_sc)))
           #obj_sc_average <- AverageExpression(obj_sc, features=rownames(sc_sig), use.scale = TRUE)
           flog.info("Need to find new signature after merging")
           sc_sig <- do.call(fun_find_markers, list(obj = obj_sc, 
            merged = TRUE, serialize = FALSE, force = TRUE, write = FALSE))
        } else {
           flog.info("NOT merging clusters %s and %s (%.2f correlation).", 
                names(markers)[p1], names(markers)[p2], max(nn)) 
            break
        }    
    }
    obj_sc[["new.idents"]] <-  Idents(obj_sc)
    obj_sc 
}

.order_clusters <- function(s) {
    s <- unique(as.character(s))
    s[order(gsub("_\\d+$","", s), as.numeric(gsub("^.*_","", s)))]
}

    
infile <- readRDS(opt$infile)
references <- readRDS(opt$singlecell)

# First integrate the single cell data with the unscaled Spatial data
scale <- if (grepl("unscaled.rds", opt$infile)) TRUE else FALSE
if (!scale) flog.info("Not scaling input file.")

integrated <- integrate_spatial(obj_spatial = infile, 
                     features = opt$num_integration_features,
                     references = references, prefix = opt$outprefix, 
                     scale = scale, min_features = opt$min_features, 
                     force = opt$force, verbose = opt$verbose)

.parse_regressout <- function() { 
    if (is.null(opt$regressout)) NULL else strsplit(opt$regressout, ":")[[1]]
}

# Now scale the integrated data and calculate UMAP
integrated <- cluster_integrated(integrated, prefix = opt$outprefix, 
                     force = opt$force, regressout = .parse_regressout(),
                     verbose = opt$verbose)

flog.info("Finished integrating data. Using %i features.", nrow(integrated))

# Cluster the single cell data
obj_sc <- cluster_reference(integrated, prefix = opt$outprefix, force = opt$force, 
                     resolution = as.numeric(strsplit(opt$resolution, ":")[[1]])[1],
                     idents = opt$idents, verbose = opt$verbose)

if (!is.null(opt$subclustering) && is.null(opt$idents)) {
    flog.info("Sub-clustering %s...", opt$subclustering)
    x <- obj_sc[, grep(opt$subclustering, Idents(obj_sc))]
    x <- cluster_reference(x, prefix = opt$outprefix, force = opt$force, 
                         resolution = as.numeric(strsplit(opt$resolution, ":")[[1]])[1],
                         sub = TRUE, verbose = opt$verbose)
    Idents(obj_sc,  cells = colnames(x)) <- as.character(Idents(x))
    obj_sc[["new.idents"]] <-  Idents(obj_sc)
}

my_find_markers_fun <- function(obj, merged = FALSE, serialize = TRUE, force = TRUE, write = TRUE) {
     find_markers(obj, references = references, 
         prefix = opt$outprefix, force = force,
         max_markers = opt$max_markers,
         resolution = as.numeric(strsplit(opt$resolution, ":")[[1]])[1],
         merged = FALSE, verbose = opt$verbose, 
         logfc.threshold = log(2), min.pct = 0.2,
         serialize = serialize, write = write)
}                     

sc_sig <- my_find_markers_fun(obj_sc, force = opt$force)

# Unless idents were provided, optimize them further 
if (is.null(opt$idents)) {

    old_idents <- Idents(obj_sc)

# Merge clusters if it's unlikely we can deconvolute them due to insufficient 
# markers
    obj_sc <- merge_clusters(obj_sc, sc_sig, max_correlation = opt$max_cor,
    min_markers = opt$min_markers, fun_find_markers = my_find_markers_fun)

# If we had to merge in the previous step, generate a new signature
    if (!identical(old_idents, Idents(obj_sc))) {
        sc_sig <- my_find_markers_fun(obj_sc, merged = TRUE)
    }
}
filename <- .get_serialize_path(opt$outprefix, "_singlecell_cluster_final.rds")
if (!opt$force && file.exists(filename)) {
    flog.warn("%s exists. Use --force to overwrite.", filename)
} else {
    flog.info("Writing R data structure to %s...", filename)
    saveRDS(obj_sc, filename)
}    
    
# Provide comprehensive cluster plots
plot_clusters(obj_sc, prefix = opt$outprefix)    

# Finally deconvolute the Spatial data
sc_deconv <- deconvolute_spatial(integrated[, integrated$technology == "spatial"], 
    prefix = opt$outprefix, 
    sig = sc_sig, hejpeg = opt$hejpeg)

find_nn_spatial(integrated[, integrated$technology == "spatial"],
    obj_sc, prefix = opt$outprefix, hejpeg = opt$hejpeg)

if (!is.null(opt$simulation)) {
    simulation <- simulate_spatial(obj_sc, infile, clusters_per_spot = 3, 
        cells_per_spot = 25, write = FALSE)
    filename <- paste0(opt$outprefix, "_simulation_qc.pdf")
    pdf(filename)
    simulation <- deconvolute_simulation(simulation, references = references, 
                       max_markers = opt$max_markers, plot_qc = TRUE,
                       min_features = opt$min_features, regressout = opt$regressout, 
                       resolution = as.numeric(strsplit(opt$resolution, ":")[[1]])[1],
                       idents = Idents(obj_sc), verbose = opt$verbose)
    dev.off()
    filename <- paste0(opt$outprefix, "_simulation_accuracy.csv")
    write.csv(data.frame(Cluster = rownames(simulation$proportions), 
                         Accuracy = simulation$evaluation), file = filename,
              row.names = FALSE)           
}
