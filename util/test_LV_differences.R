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
