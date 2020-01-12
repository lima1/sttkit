
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
