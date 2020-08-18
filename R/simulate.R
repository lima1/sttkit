#' simulate_spatial
#'
#' Simulate SpatialTranscriptomics data
#' @param obj_single_cell Seurat object containing scRNA-seq reference
#' @param obj_spatial Seurat object containing SpatialTranscriptomics 
#' template data
#' @param clusters_per_spot Assume that many clusters contribute to each spot.
#' @param cells_per_spot Assume that many cells contribute to each spot.
#' @param write Write simulated data to text files
#' @param prefix Prefix of output files
#' @export simulate_spatial
#' @examples
#' simulate_spatial()
simulate_spatial <- function(obj_single_cell, obj_spatial, clusters_per_spot, 
                             cells_per_spot, write = TRUE, prefix) {
    
    spatial_counts <- as.matrix(GetAssayData(obj_spatial, slot="counts"))
    idx_row <- intersect(rownames(spatial_counts), 
        rownames(GetAssayData(obj_single_cell, slot = "counts", assay = "RNA")))
    spatial_counts <- spatial_counts[idx_row,]
    m_ids <- apply(spatial_counts, 2, function(x) { 
        ids <- sample(clusters_per_spot, cells_per_spot, replace = TRUE)
        ids <- sample(levels(obj_single_cell))[ids]
    })
    m <- apply(m_ids, 2, function(ids) {
        cell_ids <- sapply(ids, function(x) sample(which(Idents(obj_single_cell) == x))[1])
        rowSums(as.matrix(GetAssayData(obj_single_cell, slot = "counts", assay = "RNA")[rownames(spatial_counts), cell_ids]))
    })
    
    m_prop <- apply(m_ids, 2, function(x) sapply(levels(obj_single_cell), 
        function(i) sum(i == x)))/cells_per_spot

    mm <- round( t(m)/colSums(m)*colSums(spatial_counts) ) 
    colnames(mm) <- gsub("^hg19-", "hs_", colnames(mm))
    colnames(mm) <- gsub("^mm10-", "mm_", colnames(mm))
    if (write) {
        filename <- paste0(prefix, "_singlecell_sim_prop", clusters_per_spot,
                        "_", cells_per_spot, ".csv")
        write.csv(m_prop, filename)
        filename <- paste0(prefix, "_singlecell_sim_", clusters_per_spot,
                        "_", cells_per_spot, ".tsv")
        write.table(mm, filename, sep= "\t", quote = FALSE, col.names = NA)
    }
    list(data = mm, proportions = m_prop)
}

#' evaluate_simulation
#'
#' Evalue deconvolution of simulated SpatialTranscriptomics data
#' @param simulation True cluster proportions as obtained by \code{\link{simulate_spatial}}
#' @param deconvolution Deconvoluted cluster proportions from \code{\link{deconvolute_spatial}}
#' @param plot_qc Generate concordance heatmap and scatter plots
#' @param min_fraction Minimum deconvoluted fraction to call cluster detected
#' @export evaluate_simulation
#' @examples
#' evaluate_simulation()
#' @importFrom stats cor
#' @importFrom data.table melt
evaluate_simulation <- function(simulation, deconvolution, plot_qc = TRUE, min_fraction = 0.01) {
    features <- apply(simulation$data, 1, function(x) sum(x > 0))
     
    simulation <- simulation$proportions
    colnames(deconvolution) <- make.names(colnames(deconvolution))
    colnames(deconvolution) <- gsub("_1", "", colnames(deconvolution))
    colnames(simulation) <- make.names(colnames(simulation))
    
    simulation <- simulation[, colnames(deconvolution)]
    cr <- do.call(rbind, lapply(1:nrow(simulation), function(i) 
        sapply(1:nrow(simulation), function(j) 
        cor(unlist(deconvolution[i,]), unlist(simulation[j,])))))
    colnames(cr) <- rownames(simulation)
    rownames(cr) <- rownames(simulation)
    if (plot_qc && requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(cr)
    }
    if (plot_qc) { 
        m <- merge(melt(simulation), melt(deconvolution), by=c("Var1", "Var2"))
        m$nFeature_RNA <- features[match(m$Var2, make.names(names(features)))]
        print(ggplot(m, aes_string("value.x", "value.y", color = "nFeature_RNA"))+
            geom_point()+
            viridis::scale_color_viridis(labels = function(x) sprintf("%.0f", x), 
                 trans = scales::log2_trans(),
                 breaks = scales::trans_breaks("log2", function(x) 2^x))+
            facet_wrap(~Var1)+xlab("Simulation")+ylab("Deconvolution"))
    }

    res <- sapply(seq(nrow(simulation)), function(i) 
        sum((unlist(simulation[i,]) > min_fraction & deconvolution[i,] > min_fraction) | 
            (unlist(simulation[i,]) < min_fraction & deconvolution[i,] < min_fraction))) / 
            ncol(deconvolution)
}

#' deconvolute_simulation
#'
#' Deconvolute simulated SpatialTranscriptomics data
#' @param simulation True cluster proportions as obtained by \code{\link{simulate_spatial}}
#' @param references List of Seurat scRNA-seq reference datasets 
#' @param min_features Remove cells or spots with small number of detected genes
#' @param regressout Regressout these features.
#' @param resolution Cluster resolution
#' @param idents Use idents from a serialized Seurat object or as provided
#' @param max_markers Use only top markers
#' @param plot_qc Plot concordance heatmap and scatter plot
#' @param min_fraction Minimum deconvoluted fraction to call cluster detected
#' @param verbose Verbose Seurat output
#' @export deconvolute_simulation
#' @examples
#' deconvolute_simulation()
deconvolute_simulation <- function(simulation, references, min_features, regressout, 
                                   resolution, idents, max_markers = NULL, 
                                   plot_qc = TRUE, 
                                   min_fraction = 0.01, verbose = FALSE) {
    obj_sim <- read_spatial(simulation$data, sampleid = "sim",
        serialize = FALSE)
    obj_sim <- normalize_spatial(obj_sim, scale = FALSE, serialize = FALSE)
    obj_sim_int <- integrate_spatial(obj_spatial = obj_sim, 
                         references = references, min_features = min_features,
                         force = TRUE, 
                         serialize = FALSE, verbose = verbose)
    obj_sim_int <- cluster_integrated(obj_sim_int, 
                     force = TRUE, regressout = regressout,
                     serialize = FALSE, plot_umap = FALSE, verbose = verbose)
    obj_sim_sc <- cluster_reference(obj_sim_int, force = TRUE, plot_umap = FALSE, 
                     serialize = FALSE, 
                     resolution = resolution,
                     idents = idents, verbose = verbose)
    sim_sc_sig <- find_markers(obj_sim_sc, references = references, 
        max_markers = max_markers, 
        serialize = FALSE, force = TRUE, write = FALSE)
    sim_sc_deconv <- deconvolute_spatial(obj_sim_int[, obj_sim_int$technology == "spatial"], 
        sig = sim_sc_sig, write = FALSE)
    sim_eval <- evaluate_simulation(simulation, sim_sc_deconv, min_fraction = min_fraction, plot_qc = plot_qc)
    c(simulation[c("data", "proportions")], 
      list(deconvolution = sim_sc_deconv, evaluation = sim_eval, markers = sim_sc_sig))
}
