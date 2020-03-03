
#' regex_mito
#'
#' Regular expression to robustly find mitochondria genes
#' @param hg Gene symbol prefix for human genes in hybrid references
#' @param mm Gene symbol prefix for mouse genes in hybrid references
#' @param rn Gene symbol prefix for rat genes in hybrid references
#' @export regex_mito
#' @examples
#' regex_mito()
regex_mito <- function(hg = "hg19.", mm = "mm10.", rn = "rn6.") {
    paste0("^MT\\.|^MT-|^", hg, "MT\\.|^", hg, "MT-|^mt\\.|^mt-|^", mm, "mt\\.|^", mm, 
        "mt-|^Mt-|^Mt\\.|^", rn, "Mt-|^", rn, "Mt\\.")
}

#' regex_ribo
#'
#' Regular expression to robustly find ribosomal genes
#' @param hg Gene symbol prefix for human genes in hybrid references
#' @param mm Gene symbol prefix for mouse genes in hybrid references
#' @export regex_ribo
#' @examples
#' regex_ribo()
regex_ribo <- function(hg = "hg19.", mm = "mm10.") {
    paste0("^RPS|^RPL|^", hg, "RPS|^", hg, "RPL|^Rps|^Rpl|^", mm, "Rps|^", mm, "Rpl")
}
