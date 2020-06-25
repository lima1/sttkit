#' cell-cell interaction analysis using method CellPhoneDB
#'
#' This function allows you to use an CellPhoneDB on Seurat Object
#'
#' @param Seurat_obj Seurat Object
#'
#' @return output tables for cell-cell interction inference
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, cell-cell interactions
#'
#' @examples
#'
#'
#'
#'
#' @export cellphone_for_seurat
#' @importFrom AnnotationDbi mapIds

cellphone_for_seurat <- function(obj, orgdb, prefix){
    counts <- as.data.frame(
        GetAssayData(obj, slot = "counts")
    )
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

