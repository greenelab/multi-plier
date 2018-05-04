# J. Taroni 2018
# Subsample recount2 to the same size as the SLE WB compendium & train PLIER
#
# Is the magic (i.e., performance gains over the SLE WB compendium) the *size* 
# of the recount2 compendium or the *heterogeneity* of the underlying samples
# and conditions?
#
# We can train PLIER models on a subset of the recount2 compendium that is the
# same size as the SLE WB compendium (n = 1640) to get at this question. We'll
# do that 10x here and we'll run this from the command line rather than using a 
# notebook due to its computationally intensive nature.
#
# USAGE: Rscript 11-subsample_recount_PLIER.R
#

# set seed for reproducibility
set.seed(123)

# results directory
results.dir <- file.path("results", "11")
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)

# #### functions -----------------------------------------------------------------

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

#### main ----------------------------------------------------------------------

# get a vector of 10 random seeds
seeds <- sample(1:10000, 10)

# read in recount2 data -- this will be a list that contains the
# recount2 expression matrix (RPKM), the gene set matrix, and the number
# of significant PCs for the entire gene expression matrix (k)
# k is no longer relevant here and will need to be recalculated
recount.file <- file.path("data", "recount2_PLIER_data", 
                          "recount_data_prep_PLIER.RDS")
recount.data.list <- readRDS(recount.file)

# how many total samples are in the recount2 compendium we're using?
num.cols <- ncol(recount.data.list$rpkm.cm)

# repeat this subsampling experiment 10x
for (seed in seeds) {
  
  # message for a little transparency!
  cat(paste("\n\nPerforming subsampling for seed:", seed, "\n\n\n"))
  
  # set seed -- we'll use this one for the sampling of column indices 
  # and as an argument to PLIERWrapper()
  set.seed(seed)

  # index of samples to include
  sample.index <- sample(1:num.cols, 1640)  # the number of SLE samples
  
  # expression matrix n = 1640
  smpl.exprs <- recount.data.list$rpkm.cm[, sample.index]
  
  # rescaled exprs data & PLIER model
  smpl.res <- PLIERWrapper(exprs = smpl.exprs,
                           pathway.mat = recount.data.list$all.paths.cm,
                           seed = seed)
  
  # save to file in the results directory
  smpl.file <- file.path(results.dir, 
                         paste0("recount2_subsampled_", seed, ".RDS"))
  saveRDS(object = smpl.res, file = smpl.file)
}
