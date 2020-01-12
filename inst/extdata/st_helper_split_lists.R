suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(reshape2))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help="Infile RDS file from st_normalize.R."),
    make_option(c("--outprefix"), action = "store", type = "character", default = NULL,
        help="Outfile."),
    make_option(c("--outsuffix"), action = "store", type = "character",
        default = "_spatial.list", help="Suffix for outfile list"),
    make_option(c("--hejpeg"), action = "store", type = "character", default = NULL,
        help="Optional path to a JPEG containing cropped HE image."),
    make_option(c("--labels"), action = "store", type = "character", default = NULL,
        help="Optional list of labels for multi-sample analyses"),
    make_option(c("--verbose"), action = "store_true", default = FALSE, 
        help="Verbose output"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE, 
        help="Overwrite existing files")
)


opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$infile)) {
    stop("Need --infile")
}
if (is.null(opt$outprefix)) {
    stop("Need --outprefix")
}
if (is.null(opt$labels)) {
    stop("Need --labels")
}
if (!dir.exists(dirname(opt$outprefix))) {
    dir.create(dirname(dirname(opt$outprefix)))
    dir.create(dirname(opt$outprefix))
}
if (!grepl("list$",opt$infile)) {
    stop("This tool requires input lists.")
}
suppressPackageStartupMessages(library(sttkit))

infiles <- cli_check_file_list(opt$infile)
if (length(infiles)>1) single_input <- FALSE
hejpegs <- cli_check_file_list(opt$hejpeg)
if (length(infiles) != length(hejpegs)) {
    stop("--infile as list requires --hejpeg as list.")    
}
labels <- cli_check_file_list(opt$labels, check_exists = FALSE)
if (length(infiles) != length(labels)) {
    stop("--labels requires the same length as --infile.")    
}

replicates <- sttkit:::.find_technical_replicates(labels)
labels2 <- gsub(" ", "_", gsub("_\\d+$", "", labels))
for (r in replicates) {
    cat(infiles[r], file = paste0(opt$outprefix, labels2[r[1]], opt$outsuffix), sep="\n")
    cat(gsub("_scaled.rds", ".rds", gsub("normalize", "cluster", infiles[r])), 
        file = paste0(opt$outprefix, labels2[r[1]], paste0("_cluster", opt$outsuffix)), sep="\n")
    cat(hejpegs[r], file = paste0(opt$outprefix, labels2[r[1]], "_he.list"), sep="\n")
    cat(labels[r], file = paste0(opt$outprefix, labels2[r[1]], "_labels.list"), sep="\n")
}

cat(sapply(replicates, function(r) labels2[r[1]]))
