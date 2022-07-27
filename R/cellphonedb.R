#' cell-cell interaction analysis using method CellPhoneDB
#'
#' This function allows you to use an CellPhoneDB on Seurat Object
#'
#' @param obj Seurat Object
#' @param orgdb orgdb Object
#' @param prefix Prefix of output files
#' @param slot Specific information to pull
#' (i.e. counts, data, scale.data,...)
#' @param assay Name of assay to pull data from
#'
#' @return output tables for cell-cell interction inference
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, cell-cell interactions
#'
#' @examples
#'
#' @export cellphone_for_seurat
#' @importFrom AnnotationDbi mapIds

cellphone_for_seurat <- function(obj, orgdb, prefix, slot = "data",
    assay = NULL){
    flog.info("Using counts from slot %s and assay %s.", slot,
        ifelse(is.null(assay), DefaultAssay(obj),  assay))
    counts <- as.data.frame(
        GetAssayData(obj, slot = slot, assay = assay)
    )
    counts <- .fixSymbols(counts)
    ids <- mapIds(orgdb, rownames(counts), 'ENSEMBL', 'SYMBOL')
    counts$ensembl_gene_id <- ids
    counts <- counts[!is.na(counts$ensembl_gene_id),]
    counts <- cbind(counts[,which(colnames(counts) == 'ensembl_gene_id')], counts)

    colnames(counts)[1] <- 'Gene'
    counts$ensembl_gene_id <- NULL
    colnames(counts) <- make.names(colnames(counts))

    metadata <- data.frame(Cell = make.names(Cells(obj)),
                           cell_type = make.names(Idents(obj)))


    filename <- .get_sub_path(prefix, "cellphonedb", "_counts.txt")
    flog.info("Writing count matrix to %s...", filename)
    write.table(counts,
            file = filename,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t')
    filename <- .get_sub_path(prefix, "cellphonedb", "_metadata.txt")
    write.table(metadata,
            file = filename,
            quote = FALSE,
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t')
}

.fixSymbols <- function(counts) {
    data(symbol.mapping, envir = environment())
    species_id <- which.max(apply(symbol.mapping, 2, function(x)
        length(intersect(rownames(counts),x))))
    if (species_id==1) return(counts)
    # First remove rows we cannot map to known symbols    
    row_ids <- match(rownames(counts), symbol.mapping[, species_id])

    counts <- counts[!is.na(row_ids),]
    # Now translate them; we might have cases with 1:n mappings
    # lets pick the symbol with highest mean counts
    counts <- counts[order(apply(counts,1, mean, na.rm = TRUE),decreasing = TRUE),]
    row_ids <- match(rownames(counts), symbol.mapping[, species_id])
    counts <- counts[!duplicated(symbol.mapping[row_ids,1]),]
    row_ids <- match(rownames(counts), symbol.mapping[, species_id])
    rownames(counts) <- symbol.mapping[row_ids,1]
    return(counts)
}

#' Import cell-cell interaction analysis using CellPhoneDB
#'
#' This function allows you to read CellPhoneDB output and rerank it based on local interactions.
#'
#' @param obj Seurat Object
#' @param cellphone_outpath CellPhoneDB output path
#' @param orgdb orgdb Object
#'
#' @return interaction pairs in GMT format
#'
#' @examples
#'
#' @export import_cellphone

import_cellphone <- function(obj, cellphone_outpath, orgdb) {
  if(file.exists(file.path(cellphone_outpath, "significant_means.txt"))) {
    sig_means <- read.delim(file.path(cellphone_outpath, "significant_means.txt"),
        sep="\t", as.is = TRUE, check.names = FALSE)
  } else if(file.exists(file.path(cellphone_outpath, "significant_means.csv"))) {
    sig_means <- read.delim(file.path(cellphone_outpath, "significant_means.csv"),
        sep=",", as.is = TRUE, check.names = FALSE)
  } else {
    stop("ERROR: CellPhoneDB - significant_means.(txt|csv) file not found at %s", cellphone_outpath)
  }
    col_ids <- grep("\\|", colnames(sig_means))
    sig_means <- sig_means[ apply(sig_means[, col_ids],1,function(x) sum(!is.na(x))) > 0, ]

    ids_a <- mapIds(orgdb, sig_means$gene_a, 'SYMBOL', 'ENSEMBL')
    ids_b <- mapIds(orgdb, sig_means$gene_b, 'SYMBOL', 'ENSEMBL')
    sig_means$symbol_a <- sapply(ids_a[sig_means$gene_a], function(x) ifelse(length(x), x, NA))
    sig_means$symbol_b <- sapply(ids_b[sig_means$gene_b], function(x) ifelse(length(x), x, NA))

    pairs <- apply(sig_means[, col_ids], 1, function(x) which(!is.na(x)))
    pairs <- lapply(pairs, function(x) lapply(strsplit(names(x), "\\|"), make.names))
    hood <- data.frame(find_cluster_neighborhoods(obj))
    rownames(hood) <- make.names(rownames(hood))

    if (!identical(sort(unique(colnames(hood))), sort(unique(unlist(pairs))))) {
        flog.warn("Could not map cluster ids.")
        return(sig_means)
    }

    sig_means$hood_score <- sapply(lapply(pairs, sapply, function(x)
        hood[x[1], x[2]]), function(y) 1 / (max(y) / length(y)))
    sig_means <- sig_means[order(sig_means$rank, sig_means$hood_score),]
    write.csv(sig_means, file.path(cellphone_outpath, "significant_means_ranked_spatial.csv"), row.names = FALSE)
    return(sig_means)
}

#' Plot cell-cell interaction analysis using CellPhoneDB
#'
#' This function allows you to read reranked CellPhoneDB output and plot it.
#'
#' @param obj Seurat Object
#' @param id Line of significant pair in significant_means_ranked_spatial.csv file
#' @param cellphone_outpath Path with CellPhoneDB output
#' @param features alternative genes to plot (makes sense to set \code{id} to
#' \code{NULL} in that case). 
#' @return ggplot object
#'
#' @examples
#' plot_cellphone
#' @export plot_cellphone
plot_cellphone <- function(obj, id = 1, cellphone_outpath, features = NULL) {
    if (!is.null(id)) {
        ranked <- read.csv(file.path(cellphone_outpath, "significant_means_ranked_spatial.csv"),
            as.is = TRUE)
        col_ids <- colnames(ranked)[grep("^X\\d+\\.X\\d+$", colnames(ranked))]
        pairs <- lapply(strsplit(gsub("X", "", col_ids), "\\."), as.numeric)
        cluster_pairs <- apply(t(ranked[,col_ids]),2, function(x) pairs[which(!is.na(x))])
        new_idents <- as.character(Idents(obj))
        new_idents[!Idents(obj) %in% as.character(unlist(cluster_pairs[[id]]))] <- "NS"
        Idents(obj) <- new_idents
        features <- unlist(ranked[id,c("symbol_a", "symbol_b")])
    }    
    g1 <- SpatialFeaturePlot(obj, features = features)
    g2 <- SpatialDimPlot(obj) + theme(legend.position = "top")
    if (requireNamespace("ggthemes", quietly = TRUE) &&
        length(levels(Idents(obj))) <= 8) {
        g2 <- g2 + ggthemes::scale_fill_colorblind()
    }
    return(g1+g2)
}    
   
