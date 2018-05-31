# J. Taroni 2018
# Use with source() only
# Functions for testing differential expression of PLIER latent variables

# adapted from ADAGEpath::build_limma, Jie Tan (@tj8901nm)
TestLVDifferences <- function(b.matrix, phenotype,
                              use.bonferroni = FALSE) {
  # This function tests for differential expression of PLIER-derived latent
  # variables between groups specified by the phenotype argument
  # using limma. Examples include different disease groups or disease v. 
  # control. It takes the B matrix from a PLIER model and a factor vector with 
  # the group labels. It will reorder the phenotype vector if control.phenotype
  # is specified. The resulting p-values are Benjamini-Hochberg corrected by 
  # default or Bonferroni corrected if use.bonferroni = TRUE
  # 
  # Args:
  #   b.matrix: a B matrix (latent variables are rows, samples are columns) 
  #             from a PLIER model. Can be the B element of the list 
  #             returned by PLIER::PLIER or the output of GetNewDataB
  #   phenotype: a named factor vector that contains group labels to be used
  #              for contrasts
  #   use.bonferroni: logical - should bonferroni correction be used? if FALSE
  #                   (default), will use "BH"
  #   
  # Returns:
  #   A limma::topTable, where the first column is the latent variable name and
  #   all pathways are returned (without filtering or sorting)
  
  ## error-handling ##
  
  if (is.null(names(phenotype))) {
    stop("phenotype should be a named factor vector -- the names ensure that
         the vector is correctly ordered")
  }
  
  # no names should be "missing"
  check.names <- all(colnames(b.matrix) %in% names(phenotype)) & 
    all(names(phenotype) %in% colnames(b.matrix)) 
  
  if (!check.names) {
    stop("Some sample(s) is missing from colnames(b.matrix) or 
         names(phenotype)")
  }
  
  # the phenotype labels should be in the same order as the b.matrix samples
  ordered.phenotype <- as.factor(phenotype[colnames(b.matrix)])
  
  # get contrasts (all pairwise)  
  num.levels <- length(levels(ordered.phenotype))
  contrast.vector <- c()
  for (lvl in levels(ordered.phenotype)[1:(num.levels - 1)]) {
    for (lvl.2 in levels(ordered.phenotype)[2:num.levels]) {
      if (lvl != lvl.2) {
        contrast.vector <- append(contrast.vector, paste(lvl, lvl.2, sep = "-"))
      }
    }
  }
  contrast.vector <- paste(contrast.vector, collapse = ", ")
  
  # prep design matrix
  design <- model.matrix(~ 0 + ordered.phenotype)
  colnames(design) <- levels(ordered.phenotype)
  
  # fit linear model
  fit <- limma::lmFit(b.matrix, design)
  contrast.matrix <-
    eval(parse(
      text = paste0(
        'limma::makeContrasts(', contrast.vector, ', levels = design)')))
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  
  # extract results as a data.frame
  if (use.bonferroni) {  # if specified, use Bonferroni correction
    limma.result <- limma::topTable(fit2, number = nrow(b.matrix),
                                    adjust.method = "bonferroni", 
                                    sort.by = "none")
  } else {  # calculate FDR
    limma.result <- limma::topTable(fit2, number = nrow(b.matrix),
                                    adjust.method = "BH", sort.by = "none")
  }
  
  # want feature (latent variable) names as column
  limma.result <- tibble::rownames_to_column(limma.result, var = "LV")
  
  return(limma.result)

}

BoxplotDiffLV <- function(tidy.b.df, phenotype.column, pdf.path) {
  # This function takes a LV x sample info (B matrix) in tidy form and makes
  # boxplots. Samples are grouped by the labels (levels) in the column specified
  # by the phenotype column argument. It outputs a PDF, with one page per
  # LV.
  # 
  # Args:
  #   tidy.b.df: tidy (long) data.frame that contains B matrix information and
  #              also phenotype information (in the phenotype.column)
  #   phenotype.column: column name that contains grouping for the boxplot
  #   pdf.path: full path for output PDF
  # 
  # Returns:
  #   outputs PDF specified by pdf.path
  
  ## error-handling ##
  
  pheno.check <- phenotype.column %in% colnames(tidy.b.df)
  if (!pheno.check) {
    stop("phenotype.column does not match a column name in tidy.b.df")
  }
  
  other.col.check <- all(c("LV", "Value") %in% colnames(tidy.b.df))
  if(!other.col.check) {
    stop("tidy.b.df must contain LV and Value columns")
  }
  
  pdf(pdf.path, onefile = TRUE)
  for (current.lv in unique(tidy.b.df$LV)) {
    print(dplyr::filter(tidy.b.df, LV == current.lv) %>%
            ggplot2::ggplot(ggplot2::aes_string(x = phenotype.column, 
                                                y = "Value", 
                                                group = phenotype.column)) +
            ggplot2::geom_boxplot(notch = TRUE) +
            ggplot2::geom_jitter(alpha = 0.4) +
            ggplot2::facet_wrap(~ LV) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                               hjust = 1))
    )
  }
  dev.off()
}

PrepExpressionDF <- function(exprs){
  # Takes a data.frame where the first columns contains gene identifiers.
  # Returns matrix of expression data, collapsed to gene level (mean), where
  # the gene identifiers are rownames.
  # 
  # Args:
  #   exprs: A data.frame of (normalized) gene expression data, where the
  #          rows are genes and the samples are columns. The first column
  #          should contain gene identifiers and have the column name "Gene".
  #          
  # Returns: 
  #   exprs.agg: a data.frame of aggregated expression data
  
  # error handling
  if (colnames(exprs)[1] != "Gene") {
    stop("The first column name of exprs must be named 'Gene'.")
  }
  
  `%>%` <- dplyr::`%>%` 
  
  # for duplicate gene identifiers, take the average
  exprs.agg <- exprs %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise_all(dplyr::funs(mean(., na.rm = TRUE)))
  
  # return aggregated expression data.frame
  return(exprs.agg)
  
}

LVTestWrapper <- function(b.matrix,
                          sample.info.df,
                          phenotype.col,
                          file.lead,
                          plot.dir = "plots",
                          results.dir = "results",
                          use.bonferroni = FALSE,
                          significant.only = FALSE,
                          sig.threshold = 0.05) {
  # A wrapper function for TestLVDifferences and BoxplotDiffLV; does the 
  # reshaping required for plotting. Produces the following files: 1) tsv of the
  # differential expression results 2) a long form of the B matrix joined with
  # the sample information (sample.info.df) and 3) a PDF of boxplots
  # 
  # Args:
  #   b.matrix: a B matrix (latent variables are rows, samples are columns) 
  #             from a PLIER model. Can be the B element of the list 
  #             returned by PLIER::PLIER or the output of GetNewDataB
  #   sample.info.df: a long form data.frame that contains sample information,
  #                   sample names that match the B matrix sample identifiers
  #                   must be in a "Sample" column & it also must contain
  #                   the factor to group by for testing differential expression
  #                   (DE) and plotting
  #   phenotype.col: the column name of the column in sample.info.df to be used
  #                  for DE and plotting; character
  #   file.lead: string that designates the "beginning" of filenames
  #   plot.dir: plot directory where the boxplots PDF should be saved
  #   results.dir: results directory for DE results and reshaped B data.frame
  #   use.bonferroni: logical - should bonferroni correction be used for DE? 
  #                   if FALSE (default), will use "BH"
  #   significant.only: logical - should only differentially expressed LVs be
  #                     plotted? if FALSE (default), all will be plotted
  #   sig.threshold: the adj. P cutoff to be used if only plotting significant 
  #                  results; default is 0.05
  #                  
  # Returns:
  #   A list with the following elements
  #       limma: the results from TestLVDifferences 
  #       b.df: sample.b.df, prior to any filtering for plotting (if applicable)
  # 
  #   the following files are written by this function (see above):
  #     1) <results.dir>/<file.lead>_LV_limma_results.tsv
  #     2) <results.dir>/<file.lead>_B_long_sample_info.tsv
  #     3) <plot.dir>/<file.lead>_LV_boxplots.pdf
  
  `%>%` <- dplyr::`%>%`
  
  # error-handling
  # we need to join by "Sample" column for boxplots
  if (!("Sample" %in% colnames(sample.info.df))) {
    stop("'Sample' must be a column in sample.info.df")
  }
  # phenotype.col needs to be in column names
  if (!(phenotype.col %in% colnames(sample.info.df))) {
    stop("phenotype.col must be a column name in sample.info.df")
  }
  
  # initialize list to hold results to be returned
  return.list <- list()
  
  #### Differential Expression ####
  # get the named vector to use as the phenotype for testing differential 
  # expression
  phenotype.vector <- as.factor(make.names(sample.info.df[[phenotype.col]]))
  names(phenotype.vector) <- sample.info.df$Sample
  # test itself
  limma.df <- TestLVDifferences(b.matrix = b.matrix,
                                phenotype = phenotype.vector,
                                use.bonferroni = use.bonferroni)
  # write to file
  dlve.file <- file.path(results.dir, 
                         paste0(file.lead, "_LV_limma_results.tsv"))
  readr::write_tsv(limma.df, path = dlve.file)
  # add to list to be returned
  return.list[["limma"]] <- limma.df
  
  #### Reshape & join with sample information ####
  b.df <- reshape2::melt(b.matrix)
  colnames(b.df) <- c("LV", "Sample", "Value")
  sample.b.df <- dplyr::inner_join(b.df, sample.info.df, by = "Sample")
  long.file <- file.path(results.dir, 
                         paste0(file.lead, "_B_long_sample_info.tsv"))
  readr::write_tsv(sample.b.df, long.file)
  # add to list to be returned
  return.list[["b.df"]] <- sample.b.df
  
  #### Plotting ####
  plot.file <- file.path(plot.dir,
                         paste0(file.lead, "_LV_boxplots.pdf"))
  
  # if we only want significant LVs plotted, filter sample.b.df using the
  # adj.P.Val cutoff sig.threshold
  if (significant.only) {
    sig.lvs <- limma.df$LV[which(limma.df$adj.P.Val < sig.threshold)]
    sample.b.df <- sample.b.df %>%
      dplyr::filter(LV %in% sig.lvs)
  }
  
  BoxplotDiffLV(tidy.b.df = sample.b.df, 
                phenotype.column = phenotype.col,
                pdf.path = plot.file)
  
  
  return(return.list)
}
