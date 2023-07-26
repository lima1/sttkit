.as_AssayObject_scvi <- function(object) {
    fields <- names(py_to_r(object$obsm))
    # dest_vi
    if ("proportions" %in% fields) {
        m <- as.matrix(py_to_r(object$obsm$get("proportions")))
    } else if ("q05_cell_abundance_w_sf" %in% fields) {
        m <- as.matrix(py_to_r(object$obsm$get("q05_cell_abundance_w_sf")))
        colnames(m) <- gsub("q05cell_abundance_w_sf_", "", colnames(m))
        m <- m / 100
    } else {
        stop("Not sure where to get the cell type proportions.")
    }
    m <- t(m)
    m <- rbind(m, max = apply(m, 2, max))
    return(CreateAssayObject(data = m))
}

