# J. Taroni 2018
# Given a list of prepared data saved as an RDS, this script can be run two 
# ways: 
# - Randomly selecting the specified number of samples and train a PLIER model. 
#   This process is repeated as many times is as specified using different 
#   random seeds. The output is a list of each of the repeats (named by the seed 
#   used) saved as an RDS file.
# - Given a list (TSV) of sample names/accession codes for the prepared data,
#   subset the expression data to these samples and train a model. This is not
#   repeated as we have found the models to be pretty stable.
# 
# USAGE:
#   
#   For repeated random subsampling
#   
#    Rscript scripts/subsampling_PLIER.R \
#      --input <PREPARED_DATA_RDS> \
#      --output <PATH_TO_OUTPUT_RDS> \
#      --num_samples <NUMBER_OF_SAMPLES> \
#      --repeats <NUMBER_OF_REPEATS> \
#      --seed <INITIAL_SEED>
#   
#   For subsetting to a provided list of samples
#  
#    Rscript scripts/subsampling_PLIER.R \
#      --input <PREPARED_DATA_RDS> \
#      --output <PATH_TO_OUTPUT_RDS> \
#      --use_sample_list TRUE \
#      --sample_list <PATH_TO_TSV_FILE>
#  
#  Arguments:
#    input: an RDS file that contains a list with the following elements 
#           (in order): an expression matrix, a prior information matrix -- 
#           the expression matrix and prior information should have been 
#           prepared via PLIER::commonRows
#    output: full path to the output RDS file that will contain a list of all
#            repeats -- each element of the list contains the sampled expression
#            matrix and the model
#    num_samples: the number of samples to be randomly selected (default is 500)
#    repeats: number of repeats (default is 5)
#    seed: initial seed (default is 123)
#    use_sample_list: whether or not a file with a list of sample names should
#                     be used for the subsetting (default is FALSE); if TRUE
#                     no repeats will be performed
#    sample_list: TSV file that contains a single column, without a header, of
#                 the samples to be included in the training set

#### custom functions ----------------------------------------------------------

PLIERWrapper <- function(exprs, pathway.mat, seed = 12345) {
  # "Pared down" PLIER wrapper function. The one we use elsewhere includes
  # collapsing to common genes between the expression matrix and the gene sets
  # that are included with the PLIER package -- we already have this information
  # for this experiment
  #
  # Args:
  #   exprs: a gene expression matrix, rows are genes, columns are samples
  #   pathway.mat: combined pathway matrix for use with PLIER, same
  #                genes and ordering as exprs
  #   seed: an integer to be supplied to set.seed() for reproducibility 
  #         purposes, default is 12345
  # 
  # Returns
  #   a list that contains the row-normalized expression data and the PLIER
  #   output
  
  # trouble with glmnet otherwise
  library(PLIER)
  
  set.seed(seed)
  
  # check row order of expression and pathway matrices
  ord.row.check <- all.equal(rownames(exprs), rownames(pathway.mat))
  if (ord.row.check != TRUE) {
    stop("rownames for exprs and pathway.mat should be the same.")
  }
  
  # rescale
  exprs <- PLIER::rowNorm(exprs)
  
  # what should we set the minimum k parameter to in PLIER? estimate the number
  # of PC for the SVD decomposition
  set.k <- PLIER::num.pc(exprs)
  
  # PLIER main function + return results
  plier.res <- PLIER::PLIER(as.matrix(exprs), as.matrix(pathway.mat),
                            k = round((set.k + set.k * 0.3), 0), trace = TRUE)
  return(list("exprs" = exprs,
              "PLIER" = plier.res))
  
}

#### parse command line arguments ----------------------------------------------

# loading the library seems a bit cleaner here
library("optparse")

option.list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Prepped expression and pathway matrices (.RDS)"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file (with full path) for models (.RDS)"),
  make_option(c("-n", "--num_samples"), type = "integer", default = 500,
              help = "Sample size for subsampling"),
  make_option(c("-r", "--repeats"), type = "integer", default = 5,
              help = "Number of repeats to perform"),
  make_option(c("-s", "--seed"), type = "integer", default = 123,
              help = "Number of repeats to perform"),
  make_option(c("-u", "--use_sample_list"), type = "logical", default = FALSE,
              help = "Number of repeats to perform"),
  make_option(c("-l", "--sample_list"), type = "character", default = NULL,
              help = "TSV file that contains samples to subset to")
)

opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)

input.file <- opt$input
output.file <- opt$output
sample.size <- opt$num_samples
num.repeats <- opt$repeats
initial.seed <- opt$seed
use.sample.list <- opt$use_sample_list
sample.file <- opt$sample_list

# if the arguments specify that we should be using a sample list and there's no
# sample list provided, throw an error
if (use.sample.list & is.null(sample.file)) {
  stop("use_sample_list == TRUE and no sample list provided!")
}

# if the output directory doesn't exist, create it
output.directory <- stringr::word(string = output.file, start = 1, 
                                  end = -2, sep = .Platform$file.sep)
dir.create(output.directory, showWarnings = FALSE, recursive = TRUE)

#### main ----------------------------------------------------------------------

# read in the prepped data list
# the expression data is the first element of prepped.data
# the pathways prior information matrix is the second element of prepped.data
prepped.data <- readRDS(input.file)

# subsetting to using a list of samples
if (use.sample.list) {
  
  # read in the TSV file that contains the names or accession codes of the 
  # samples that should be included in training
  sample.df <- readr::read_tsv(sample.file, col_names = FALSE)
  
  # find the indices
  sample.index <- which(colnames(prepped.data[[1]]) %in% sample.df[[1]])
  
  # if not all the samples are found, throw an error
  if (length(sample.index) != nrow(sample.df)) {
    stop("Sample index length != number of listed samples")
  }
  
  # expression matrix that includes the specified samples
  smpl.exprs <- prepped.data[[1]][, sample.index]
  
  plier.results <- PLIERWrapper(exprs = smpl.exprs,
                              pathway.mat = prepped.data[[2]],
                              seed = initial.seed)
  
  # save to file
  saveRDS(object = plier.results, file = output.file)
  
} else {  # Random subsampling
  
  # set seed for reproducibility
  set.seed(initial.seed)
  
  # we'll need seeds for each of the repeats
  seeds <- sample(1:10000, num.repeats)
  
  # how many total samples are in the recount2 compendium we're using?
  num.cols <- ncol(prepped.data[[1]])
  
  # initialize list to hold models
  model.list <- list()
  # repeat this subsampling experiment 10x
  for (seed in seeds) {
    
    # message for a little transparency!
    cat(paste("\n\nPerforming subsampling for seed:", seed, "\n\n\n"))
    
    # set seed -- we'll use this one for the sampling of column indices 
    # and as an argument to PLIERWrapper()
    set.seed(seed)
    
    # index of samples to include
    sample.index <- sample(1:num.cols, sample.size)
    
    # expression matrix n = samples.size
    smpl.exprs <- prepped.data[[1]][, sample.index]
    
    # rescaled exprs data & PLIER model
    model.list[[as.character(seed)]] <- 
      PLIERWrapper(exprs = smpl.exprs,
                   pathway.mat = prepped.data[[2]],
                   seed = seed)
  }
  
  # save to file
  saveRDS(object = model.list, file = output.file)
  
}