# Name: get_distances_between_overlaps.R
# Description: This script is designed to obtain the distance between the 5' end of the reduced GFFs and
# the 5' end of the manual anotation, for every element of the genome analysis. The same for the 3' end.
# Author: David Felipe Rendón Luna 
# Date: August-2023


# Check for optparse library to load arguments from command line ----------
if(suppressMessages(!require("optparse"))) {
  stop("optparse was not found. Exiting.")
}

# Check for aditional libraries
if(nzchar(system.file(package = "rtracklayer"))) {
  cat("-rtracklayer library found.\n")
} else {
  stop("rtracklayer was not found. Exiting.\n")
}

if(nzchar(system.file(package = "GenomicRanges"))) {
  cat("-GenomicRanges library found.\n")
} else {
  stop("GenomicRanges was not found. Exiting.\n")
}

# Load parser -------------------------------------------------------------
library("optparse")
# Create list of arguments and help asociadted to receive.
opt_list = list(make_option(opt_str = c("-q", "--query"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "Manual anotation with the genetic feature to be analyzed, as a BED formated file. \n
                            For cleaner results, filter your BED file for the desired genetic feature.", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-s", "--subject"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "Reduced GFF file, produced by gff_disambiguation.R script.", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-f", "--feature"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "The name of the genetic feature to be analized.", 
                            metavar = "[STRING_FEATURE]"),
                make_option(opt_str = c("-o", "--outfile"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "Filename for the plot to be produced.", 
                            metavar = "[FILENAME]"))

# Make the parser
opt_parser = OptionParser(option_list = opt_list)

# Load the arguments in the parser into an object
opt = parse_args(opt_parser)

# Check for arguments -----------------------------------------------------
if (length(opt) == 1) {
  print_help(opt_parser)
  stop("No arguments submitted")
}

if (is.null(opt$query)) {
  print_help(opt_parser)
  stop("A query file must be submitted", call.=FALSE)
}

if (is.null(opt$subject)) {
  print_help(opt_parser)
  stop("A subject file must be submitted", call.=FALSE)
}

if (is.null(opt$feature)) {
  print_help(opt_parser)
  stop("A genetic feature must be submitted", call.=FALSE)
}

if (is.null(opt$outfile)) {
  print_help(opt_parser)
  stop("A outfile name must be submitted", call.=FALSE)
}

# SCRIPT
library(utils)

# PREPARE QUERY AND SUBJECT

# Import BED file of manual anotation
print("Reading manual anotation.")
manual_anotation <- rtracklayer::import(opt$query, format = "BED")

# Edit names and get the element code
gene_names <- sapply(manual_anotation@elementMetadata@listData$name, function(x) gsub('\\|.*', '', x ), USE.NAMES = FALSE)

# Substitute each element with only the element code
manual_anotation@elementMetadata@listData$name <- gene_names 

# Import reduced GFF file
print("Reading reduced GFF file.")
reduced_gff <- rtracklayer::import(opt$subject, format = "GFF")

# Make overlaps
overlaps <- GenomicRanges::findOverlaps(query = manual_anotation, subject = reduced_gff)

# Initialice empty vector for results
all_dif_max <- vector(mode = "integer")
all_dif_min <- vector(mode = "integer")

# ANALICE DATA
for (GENE_NAME in unique(gene_names)) {
  # Get the position in the manual anotation that matches each element
  positions_element <- which(manual_anotation@elementMetadata@listData$name == GENE_NAME)
  
  # Filter the manual anotation to get every range of one element
  temp_gobj <- manual_anotation[positions_element]
  
  # Using the start and the width of the elements, get the higher and lower value
  manual_max <- max(c(temp_gobj@ranges@start, temp_gobj@ranges@start + temp_gobj@ranges@width -1)) 
  manual_min <- min(c(temp_gobj@ranges@start, temp_gobj@ranges@start + temp_gobj@ranges@width -1)) 
  
  # Use position element to extract the positions of only the querys that have matched
  positions_element_overlaps_sub <- which(overlaps@from %in% positions_element)
  
  #Obtain the element that makes the overlap
  temp_gobj_reduced <- reduced_gff[unique(overlaps[positions_element_overlaps_sub]@to)]
  
  # Using the start and the width of the elements, get the higher and lower value
  reduced_max <- suppressWarnings(max(c(temp_gobj_reduced@ranges@start, temp_gobj_reduced@ranges@start + temp_gobj_reduced@ranges@width -1))) 
  reduced_min <- suppressWarnings(min(c(temp_gobj_reduced@ranges@start, temp_gobj_reduced@ranges@start + temp_gobj_reduced@ranges@width -1)))
  
  # Get diference of maximums and minimums
  dif_max <- reduced_max - manual_max
  dif_min <- reduced_min - manual_min
  
  if(dif_min == Inf | dif_max == -Inf) {
    #print(GENE_NAME)
    next()
  }
  
  # Store obtained values
  all_dif_max <- c(all_dif_max, dif_max)
  all_dif_min <- c(all_dif_min, dif_min)
}
                     
print("testing")
print(opt$outfile)
print(typeof(opt$outfile))
                     
outname <- paste0(opt$outfile, "_", opt$feature, ".png")

png(outname, width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
data_hist1 <- hist(all_dif_min, plot = FALSE)
hist(all_dif_min, 
     breaks = seq(10 * floor(min(all_dif_min) / 10), 10 * ceiling(max(all_dif_min) / 10), 10), 
     xlim = c(min(all_dif_min), max(all_dif_min)), 
     xaxt = "n",
     yaxt = "n",
     main = 'Min', 
     xlab = 'Distance in bp', 
     sub = paste0("N = ", length(all_dif_min)),
     col = "blue")
axis(side = 1, at = seq(10 * floor(min(all_dif_min) / 10), 10 * ceiling(max(all_dif_min) / 10), 20))
axis(side = 2, las = 2,
     labels = seq(0, 10 * ceiling(max(data_hist1$counts) / 10), 10), 
     at = seq(0, 10 * ceiling(max(data_hist1$counts) / 10), 10))
abline(v = 0, col = 'red', lwd = 2, lty = 'dashed')

data_hist2 <- hist(all_dif_max, plot = FALSE)
hist(all_dif_max, breaks = seq(10 * floor(min(all_dif_max) / 10), 10 * ceiling(max(all_dif_max) / 10), 10), 
     xlim = c(min(all_dif_max), max(all_dif_max)), 
     xaxt = "n",
     yaxt = "n",
     main = 'Max', 
     xlab = 'Distance in bp', 
     sub = paste0("N = ", length(all_dif_max)),
     col = "orange")
axis(side = 1, at = seq(10 * floor(min(all_dif_max) / 10), 10 * ceiling(max(all_dif_max) / 10), 20))
axis(side = 2,  las = 2,
     labels = seq(0, 10 * ceiling(max(data_hist2$counts) / 10), 10), 
     at = seq(0, 10 * ceiling(max(data_hist2$counts) / 10), 10))
abline(v = 0, col = 'red', lwd = 2, lty = 'dashed')

dev.off()
