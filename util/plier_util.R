# J. Taroni 2018
# For import only: source(file.path("util", "plier_util.R"))

PLIERNewData <- function(exprs.mat, seed = 12345) {
  # A wrapper function for applying PLIER to a data set. We use the following
  # genesets that come with PLIER: bloodCellMarkersIRISDMAP, svmMarkers, 
  # and canonicalPathways. We set the k parameter for the PLIER model by
  # identifying the number of "significant PCs" with PLIER::num.pc and then 
  # using sig PCs * 0.3. This is consistent with recommendations from the 
  # PLIER authors.
  # 
  # Args:
  #   exprs.mat: a gene expression matrix, rows are genes, columns are samples
  #   seed: an integer to be supplied to set.seed() for reproducibility 
  #         purposes, default is 12345
  #         
  # Returns:
  #   plier.res: output from PLIER::PLIER()
  #
  require(PLIER)
  
  set.seed(seed)
  
  # load PLIER pathway and cell type data
  data(bloodCellMarkersIRISDMAP)
  data(svmMarkers)
  data(canonicalPathways)
  
  # combine the pathway data from PLIER
  all.paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, 
                                   canonicalPathways)
  
  # what genes are common to the pathway data and the expression matrix
  cm.genes <- PLIER::commonRows(all.paths, exprs.mat)
  
  # row normalize
  exprs.norm <- PLIER::rowNorm(exprs.mat)
  
  # what should we set the minimum k parameter to in PLIER? estimate the number 
  # of PC for the SVD decomposition 
  set.k <- PLIER::num.pc(exprs.norm[cm.genes, ])
  
  # PLIER main function + return results
  plier.res <- PLIER::PLIER(exprs.norm[cm.genes, ], all.paths[cm.genes, ], 
                            k = round((set.k + set.k * 0.3), 0), trace = TRUE)
  
  return(plier.res)
  
}

GetOrderedRowNorm <- function(exprs.mat, plier.model) {
  # Given an expression matrix of data that was not used to train a plier model
  # and the output of PLIER::PLIER, row-normalize, remove genes not in the 
  # training data, set missing genes to zero (the mean), and reorder to match
  # plier.model$Z
  # 
  # This makes the input gene expression data suitable for projection into
  # the training latent variable space (see GetNewDataB) and for evaluating 
  # reconstruction (MASE, correlation)
  # 
  # Args:
  #   exprs.mat: a gene expression matrix, rows are genes, columns are samples
  #   plier.model: PLIER results, output from PLIER::PLIER()
  # 
  # Returns:
  #   ord.rownorm: a row normalized (z-scored) gene expression matrix that
  #                matches the Z matrix of the PLIER model -- ordered to match,
  #                contains the same genes
  
  # first, row normalize the new expression data
  exprs.norm <- PLIER::rowNorm(exprs.mat)
  
  # get Z matrix from PLIER model
  z.mat <- plier.model$Z
  # what genes were used in the model?
  genes.in.model <- rownames(z.mat)
  
  # get the genes that are common to the PLIER model and the new expression
  # data
  exprs.cg <- exprs.norm[which(rownames(exprs.norm) %in% genes.in.model), ]
  
  # add in genes that are missing in the new exprs data 
  genes.not.exprs <- 
    genes.in.model[which(!(genes.in.model %in% rownames(exprs.norm)))]
  # set all to zero -- this is the mean, due to z-scoring
  miss.mat <- matrix(0, ncol = ncol(exprs.cg), nrow = length(genes.not.exprs))
  # set gene names (rownames) to missing gene names
  rownames(miss.mat) <- genes.not.exprs
  # set colnames to the same as the expression matrix for genes present in
  # exprs.mat
  colnames(miss.mat) <- colnames(exprs.cg)
  # add into common gene expression matrix
  exprs.cg <- rbind(exprs.cg, miss.mat)
  # reorder rows
  ord.rownorm <- exprs.cg[genes.in.model, ]
  
  # check reordering
  gene.ord.chk <- all.equal(rownames(ord.rownorm), rownames(z.mat))
  
  if (gene.ord.chk) {
    return(ord.rownorm)
  } else {
    stop("Something went wrong -- Z matrix gene order doesn't match with
         the ordered gene expression matrix")
  }
  
}

GetNewDataB <- function(exprs.mat, plier.model) {
  # Apply a previously computed PLIER to a new dataset to get the LV x sample
  # matrix (B)
  # see main PLIER function: 
  # https://github.com/wgmao/PLIER/blob/a2d4a2aa343f9ed4b9b945c04326bebd31533d4d/R/Allfuncs.R#L227
  # 
  # Args:
  #   exprs.mat: a gene expression matrix, rows are genes, columns are samples
  #   plier.model: PLIER results, output from PLIER::PLIER()
  # 
  # Returns:
  #   exprs.new.b: a matrix that contains the values of each latent variable
  #                for each sample from the new dataset (exprs.mat), 
  # 
  require(PLIER)
  
  # the error handling for reordering failing is now in GetOrderedRowNorm 
  ord.rownorm <- GetOrderedRowNorm(exprs.mat, plier.model) 
  
  # get Z matrix from PLIER model
  z.mat <- plier.model$Z
  
  # get LV by sample (B) matrix from the new exprs using PLIER model
  # https://github.com/wgmao/PLIER/blob/a2d4a2aa343f9ed4b9b945c04326bebd31533d4d/R/Allfuncs.R#L465
  exprs.new.b <-
    solve(t(z.mat) %*% z.mat + plier.model$L2 * diag(ncol(z.mat))) %*%
    t(z.mat) %*% ord.rownorm
  
  # add in the rownames from the PLIER model B matrix
  rownames(exprs.new.b) <- rownames(plier.model$B)
  # return B matrix
  return(exprs.new.b)
  
}

GetReconstructedExprs <- function(z.matrix, b.matrix) {
  # Get "reconstructed" expression data from PLIER loadings and latent variables
  # 
  # Args:
  #   z.matrix: a matrix of gene loadings from PLIER
  #   b.matrix: a matrix of latent variables from PLIER
  #   
  # Returns:
  #   recon.exprs: the reconstructed expression data
  recon.exprs <- z.matrix %*% b.matrix
  return(recon.exprs)
}

#### compare reconstructed to input --------------------------------------------

GetReconstructionMASE <- function(true.mat, recon.mat){
  # This function takes two gene expression matrices (before and after 
  # reconstruction) and compares them.
  # The reconstruction error between the two (MASE: mean absolute scaled error) 
  # is calculated on a per sample basis.
  # 
  # Args:
  #   true.mat: gene expression matrix before reconstruction, samples are rows
  #             genes are columns
  #   recon.mat: gene expression matrix after reconstruction, samples are rows
  #              genes are columns
  # 
  # Returns:
  #   mase: MASE before and after reconstruction (per sample basis)
  # 
  
  # Error-handling  
  if (all.equal(dim(true.mat), dim(recon.mat)) != TRUE){
    stop("true.mat and recon.mat must be of equal dimensions")
  }
  
  if (!all.equal(colnames(true.mat), colnames(recon.mat))) {
    stop("colnames(true.mat) and colnames(recon.mat) must be equal")
  }
  
  CalculateMASE <- function(y, y.pred){
    # For a pair of vectors, a sample's measured expression value (y) and 
    # a sample's reconstructed expression value (y.pred), calculate the sample
    # mean absolute scaled error.
    # number of observations
    n <- length(y)
    # difference in gene expression and reconstructed gene expression
    abs.errors <- abs(y - y.pred)
    # mean of true expression values
    y.bar <- mean(y)
    # calculate absolute scaled error
    scaled.errors <- abs.errors/(sum(abs(y - y.bar)) / n)
    # calculate mean absolute scaled error
    mase <- mean(scaled.errors)
    return(mase)
  }
  
  # for each gene (column), calculate the MASE between the true expression 
  # values and the expression values after reconstruction
  mase <- sapply(1:ncol(true.mat), 
                 function(x) CalculateMASE(y = true.mat[, x], 
                                           y.pred = recon.mat[, x]))
  
  # return a vector of mean absolute scaled errors
  return(mase)
  
}

GetReconstructionCorrelation <- function(true.mat, recon.mat,
                                         cor.method = "spearman"){
  # This function takes two gene expression matrices (before and after 
  # reconstruction) and calculates the correlation between the input and
  # reconstructed values (by sample)
  # 
  # Args:
  #   true.mat: gene expression matrix before reconstruction, samples are rows
  #             genes are columns
  #   recon.mat: gene expression matrix after reconstruction, samples are rows
  #              genes are columns
  #   cor.method: the method to be used with cor(), default is "spearman"
  # 
  # Returns:
  #   a vector of correlation values between the input and reconstructed 
  #   samples, using cor.method
  # 
  
  # Error-handling  
  if (all.equal(dim(true.mat), dim(recon.mat)) != TRUE){
    stop("true.mat and recon.mat must be of equal dimensions")
  }
  
  if (!all.equal(colnames(true.mat), colnames(recon.mat))) {
    stop("colnames(true.mat) and colnames(recon.mat) must be equal")
  }
  
  cor.vect <- vector()
  for (col.iter in 1:ncol(true.mat)){
    # for each gene (column), calculate the MASE between the true expression 
    # values and the expression values after reconstruction
    # due to size of matrices, this is more efficient than calculating 
    # correlation between all values
    cor.vect[col.iter] <- cor(true.mat[, col.iter],
                              recon.mat[, col.iter],
                              method = cor.method)
  }
  
  # return a vector of mean absolute scaled errors
  return(cor.vect)
  
}

#### sparsity and pathway coverage ---------------------------------------------

GetPathwayCoverage <- function(plier.results, fdr.cutoff = 0.05) {
  # This function calculates the proportion of pathways (input into PLIER
  # as prior information) that are significantly associated with a latent
  # variable (to.return = "pathway"), the proportion of LVs that have a 
  # pathway significantly associated with them (to.return = "LV") and 
  # the number of significant pathways / number of LV
  # 
  # Args:
  #   plier.results: output of main PLIER function (PLIER::PLIER) -- this will
  #                  contain the information needed to calculate pathway 
  #                  coverage -- the significance of association between
  #                  LV and pathways (summary), the prior info coef matrix (U),
  #                  and the prior information supplied to PLIER (C)
  #   fdr.cutoff: what is the FDR cutoff for significance? 0.05 by default
  #                       
  # Returns:
  #   A list with the following elements:
  #     pathway: proportion of "covered" pathways
  #     lv: proportion of LVs that have pathways associated with them
  #     sig.pathway.by.lv: number of pathways "covered" / number of LVs
  #     
  
  # extract relevant information from the PLIER results
  summary.df <- plier.results$summary
  input.pathways <- colnames(plier.results$C)
  num.lvs <- ncol(plier.results$U)

  # identify the pathways that have significant associations with latent
  # variables, using the value from the fdr.cutoff argument
  sig.pathways <- 
    unique(summary.df$pathway[which(summary.df$FDR < fdr.cutoff)])
  
  # identify the LV that are (significantly) associated with pathways
  sig.lvs <- unique(summary.df$`LV index`[which(summary.df$FDR < fdr.cutoff)])
  
  # initialize a list to hold the results
  return.list <- list()
  
  # proportion of "covered" pathways
  return.list$pathway <- 
    length(sig.pathways) / length(input.pathways)
  
  # proportion of LVs that have pathways associated with them
  return.list$lv <- length(sig.lvs) / num.lvs
  
  # number of pathways "covered" / number of latent variables
  return.list$sig.pathway.by.lv <- length(sig.pathways) / num.lvs
  
  return(return.list)
  
}

CalculateUSparsity <- function(plier.results, significant.only = FALSE,
                               fdr.cutoff = 0.05) {
  # This function takes the output of PLIER::PLIER (plier.results) and 
  # calculates what proportion of each of the columns of plier.results$U 
  # (prior information coefficient matrix) have positive entries. 
  # This can be performed using all positive entries (significant.only = FALSE,
  # the default) or using only pathways/priori that are significantly
  # associated with an LV(signifcant.only = TRUE). We use plier.results$summary 
  # and the specified FDR cutoff (0.05 by default) to
  # determine with LVs are significantly associated with a pathway or gene set.
  # 
  # Args:
  #   plier.results: output of main PLIER function (PLIER::PLIER), contains the
  #                  U (prior information coefficient) matrix and summary 
  #                  data.frame required for the calculation
  #   significant.only: logical - should only pathways that are significantly 
  #                     associated with an LV information be taken 
  #                     into account? default is FALSE -> all positive entries
  #                     will be used in the calculation
  #   fdr.cutoff: the FDR cutoff used to determine which LVs should be taken
  #               into account if significant.only = TRUE
  # 
  # Returns:
  #   prop.vector: a vector of length ncol(plier.results$U) that contains the
  #                proportion of positive (significant.only = FALSE) 
  #                or significant (significant.only = TRUE) entries in each 
  #                column
  #   
  
  # extract relevant elements from plier.results
  # summary
  summary.df <- plier.results$summary
  # U matrix
  data.mat <- plier.results$U
  no.entries <- nrow(data.mat)  # how many pathways/gene sets were included?
  
  # MAIN: calculate proportion of positive or "significant" entries in each
  # column of U (data.mat)
  if (!significant.only) {  # if using all entries
    prop.vector <- apply(data.mat, 2, function(x) sum(x > 0) / no.entries )
  } else if (significant.only) {  # if using only significant LVs
    # initialize vector to hold values
    prop.vector <- vector()
    # filter to only significant entries of summary.df using the fdr.cutoff 
    # argument
    summary.df <- dplyr::filter(summary.df, FDR < fdr.cutoff)
    # for each LV (column of U)
    for (col.iter in 1:ncol(data.mat)) {
      # first, identify the significant pathways for each LV
      if (col.iter %in% summary.df$`LV index`) {  # LV index in filt.summary.df
        # number of pathways
        prop.vector[col.iter] <- 
          sum(summary.df$`LV index` == col.iter) / no.entries
      } else {  # no significant pathways associated with this LV
        prop.vector[col.iter] <- 0
      }
    }
  } else {
    stop("significant.only argument must be TRUE or FALSE")
  }
  
  # return proportion information
  return(prop.vector)
}

EvalWrapper <- function(plier.model){
  # Wrapper function for evaluation of PLIER models, takes a PLIER model 
  # (output of PLIER::PLIER, typically via the custom function 
  # PLIERNewData) and outputs the following metrics: 1) U sparsity (all LVs) 
  # 2) U sparsity (only significant LVs taken into account) 3) pathway coverage
  # metrics 4) number of LVs
  #
  # Args:
  #   plier.model: output of PLIER::PLIER, typically via the custom function 
  #                PLIERNewData
  #
  # Returns:
  #   A list with the following elements 1) all.sparsity 2) sig.sparsity 
  #   3) pathway.coverage 4) num.lvs -- see above'
  
  # Initialize list
  return.list <- list()
  # U sparsity
  return.list$all.sparsity <- CalculateUSparsity(plier.model)
  return.list$sig.sparsity <- CalculateUSparsity(plier.model, 
                                                 significant.only = TRUE,
                                                 fdr.cutoff = 0.05)
  # pathway coverage metrics (3 of them)
  return.list$pathway.coverage <- GetPathwayCoverage(plier.model, 
                                                     fdr.cutoff = 0.05)
  # number of latent variables
  return.list$num.lvs <- nrow(plier.model$B)
  
  return(return.list)
  
}

#### pathway holdout -----------------------------------------------------------

CalculateHoldoutAUC <- function(plier.result, holdout.mat, ncores = 6) {
  # This function is adapted from PLIER:::crossVal. It will calculate the
  # AUC and p-value for each pathway in a heldout prior information matrix for
  # each latent variable in the supplied PLIER model.
  # 
  # Args:
  #   plier.result: PLIER model, output of PLIER::PLIER
  #   holdout.mat: a prior information matrix, formatted as is passed to
  #                PLIER::PLIER -- genes are rows, pathways are columns, values
  #                are binary (0/1) where 1 indicates that a gene is in a 
  #                pathway
  #   ncores: number of cores to use for parallel backend, 6 by default
  #
  # Returns:
  #  auc.df: a data.frame that contains the pathway, LV index, AUC, p-value,
  #          and FDR for each heldout pathway-LV index pair
  
  require(foreach)
  `%>%` <- dplyr::`%>%`
  
  # what genes were used in the model?
  genes.in.model <- rownames(plier.result$Z)
  
  # get the genes that are common to the PLIER model and the oncogenic pathways
  holdout.cg <- holdout.mat[which(rownames(holdout.mat) %in% genes.in.model), ]
  
  # add in genes that are missing in the new oncogenic pathway data 
  genes.not.path <- 
    genes.in.model[which(!(genes.in.model %in% rownames(holdout.mat)))]
  # set all to zero -- that means they're not in the pathway
  miss.mat <- matrix(0, ncol = ncol(holdout.cg), nrow = length(genes.not.path))
  # set gene names (rownames) to missing gene names
  rownames(miss.mat) <- genes.not.path
  # add the "missing" genes
  colnames(miss.mat) <- colnames(holdout.mat)
  holdout.cg <- rbind(holdout.cg, miss.mat)
  # reorder rows
  ord.holdout <- holdout.cg[genes.in.model, ]
  
  # parallel backend 
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  # adapted from PLIER:::crossVal -- in that function, genes that are either
  # not in the pathway or are in the pathway but are heldout are supplied as 
  # labels to PLIER:::AUC 
  # here, all genes are held out, so we pass the entire set and we do this
  # for each combination of pathway (column in ord.holdout) and latent variable 
  # (column in Z from the PLIER model)
  auc.list <- foreach(path.iter = 1:ncol(ord.holdout)) %do% {
    foreach(z.iter = 1:ncol(plier.result$Z)) %dopar% {
      PLIER:::AUC(ord.holdout[, path.iter], plier.result$Z[, z.iter])
    }
  }
  
  # stop parallel backend
  parallel::stopCluster(cl)
  
  # the outermost element of this list is the heldout pathways
  names(auc.list) <- colnames(ord.holdout)
  
  # melting and wrangle to match the data.frames returned from PLIER::PLIER / 
  # PLIER:::crossVal
  auc.df <- reshape2::melt(auc.list) %>%
    dplyr::filter(L3 %in% c("pval", "auc")) %>%
    tidyr::spread(L3, value)
  colnames(auc.df) <- c("LV index", "pathway", "AUC", "p-value")
  
  # reorder to match output and calculate FDR
  auc.df <- auc.df %>%
    dplyr::select(c("pathway", "LV index", "AUC", "p-value")) %>%
    dplyr::mutate(FDR = PLIER:::BH(`p-value`))
  
  # return the full data.frame (no filtering for "significance")
  return(auc.df)
}

EvalWrapperWithHoldout <- function(plier.model, holdout.matrix, 
                                   ncores = 6, fdr.cutoff = 0.05) {
  # Given a PLIER model, calculate pathway coverage (w/ GetPathwayCoverage), 
  # find the number of latent variables in the model, and calculate the AUC 
  # for the heldout pathways in holdout.matrix (w/ CalculateHoldoutAUC); 
  # return a list of these results
  #
  # Args:
  #  plier.model: output of main PLIER function (PLIER::PLIER)
  #  holdout.matrix: a prior information matrix, formatted as is passed to
  #                  PLIER::PLIER -- genes are rows, pathways are columns, 
  #                  values are binary (0/1) where 1 indicates that a gene 
  #                  is in a pathway
  #   ncores: number of cores to use for parallel backend in holdout AUC 
  #           calculations, 6 by default
  #   fdr.cutoff: what is the FDR cutoff for significance (for pathway
  #               coverage calculations)? 0.05 by default
  #   
  # Returns a list with the following elements:
  #    pathway.coverage: output of GetPathwayCoverage; see that function's
  #                      documentation for more information
  #    num.lvs: number of latent variables in the model
  #    heldout.results: output of CalculatedHoldoutAUC; see that function's
  #                     documentation for more information
  #
  
  # pathway coverage
  pathway.coverage <- GetPathwayCoverage(plier.results = plier.model, 
                                         fdr.cutoff = fdr.cutoff)
  
  # number of latent variables
  num.lvs <- nrow(plier.model$B)
  
  # AUC, FDR, etc. for pathways in holdout.matrix
  heldout.results <- CalculateHoldoutAUC(plier.result = plier.model,
                                         holdout.mat = holdout.matrix,
                                         ncores = ncores)
  
  # return a list of all results
  return(list(pathway.coverage = pathway.coverage,
              num.lvs = num.lvs,
              heldout.results = heldout.results))
  
}