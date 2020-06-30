#' cell-cell interaction analysis using method CellPhoneDB
#'
#' This function allows you to use an CellPhoneDB on Seurat Object
#'
#' @param obj Seurat Object
#' @param orgdb orgdb Object
#' @param prefix Prefix of output files
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
                           cell_type = Idents(obj))


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
