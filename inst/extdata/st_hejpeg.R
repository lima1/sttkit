suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action = "store", type = "character", default = NULL,
        help="Infile JPEG with HE."),
    make_option(c("--outfile"), action = "store", type = "character", default = NULL,
        help="Outfile JPEG optimized for STT."),
    make_option(c("--scale"), action = "store", type = "character", default = "3000x",
        help="Scale JPEG [default %default]."),
    make_option(c("--brightness"), action = "store", type = "double", default = 200,
        help="Adjust brightness [default %default]."),
    make_option(c("--saturation"), action = "store", type = "double", default = 100,
        help="Adjust saturation [default %default]."),
    make_option(c("--hue"), action = "store", type = "double", default = 100,
        help="Adjust hue [default %default]."),
    make_option(c("--dither"), action = "store_true", default = NULL, 
        help="Apply Floyd/Steinberg error diffusion to the image."),
    make_option(c("--skip_grayscale"), action = "store_true", default = FALSE, 
        help="Skip grayscale step"),
    make_option(c("--skip_rescale"), action = "store_true", default = FALSE, 
        help="Skip rescale step"),
    make_option(c("--skip_brightness"), action = "store_true", default = FALSE, 
        help="Skip rescale step"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE, 
        help="Overwrite existing files")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (!opt$force && file.exists(opt$outfile)) {
    stop(opt$outfile, " exists. Use --force to overwrite.")
}

if (is.null(opt$infile)) {
    stop("Need --infile")
}

suppressPackageStartupMessages(library(magick))
flog.info("Reading %s...", opt$infile)
image <- image_read(opt$infile)

if (!opt$skip_rescale) { 
    flog.info("Rescaling to %s...", opt$scale)
    image <- image_scale(image, opt$scale)
}    
if (!opt$skip_grayscale) {
    flog.info("Converting to grayscale...")
    image <- image_quantize(image, max = 256, colorspace = "gray", dither = opt$dither)
}    
if (!opt$skip_brightness) {
    flog.info("Brightness to %.0f, Saturation to %.0f and Hue to %.0f...", 
        opt$brightness, opt$saturation, opt$hue)
    image <- image_modulate(image, brightness = opt$brightness, 
        saturation = opt$saturation, hue = opt$hue)
}    
flog.info("Writing JPEG image to %s...", opt$outfile)
image_write(image, path = opt$outfile, format = "jpeg")
