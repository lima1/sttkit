
#' cli_check_file_list
#'
#' Internal function for the command line interface
#' @param file Input .list file, alternatively a string
#' with filenames separated by :
#' @param check_exists Check input files
#' @export cli_check_file_list
#' @examples
#' cli_check_file_list()
cli_check_file_list <- function(file, check_exists = TRUE) {
    if (file.exists(file)) {
        files <- read.delim(file, as.is = TRUE, header = FALSE)[,1]
        files <- trimws(files, which = "both")
    } else {
        files <- strsplit(file, ":")[[1]]
    }    
    if (check_exists) {
        num_exists <- sum(file.exists(files), na.rm = TRUE)
        if (num_exists < length(files)) { 
            stop("File not exists in file ", file)
        }
    }
    files
}

#' cli_check_lib_ids
#'
#' Internal function for the command line interface
#' @param reference_list List of \code{Seurat} objects
#' @export cli_check_lib_ids
#' @examples
#' cli_check_lib_ids
cli_check_lib_ids <- function(reference_list) {
    libs <- sapply(reference_list, function(x) x$library[1])
    if (any(duplicated(libs))) {
        flog.warn("Found duplicated library ids, appending indices to ids.")
        reference_list <- lapply(seq_along(reference_list), function(i) {
            reference_list[[i]]$library <- paste(reference_list[[i]]$library, i, sep = "_")
            reference_list[[i]]
        })
    }
    reference_list
}    
