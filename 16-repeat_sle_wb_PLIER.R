# J. Taroni 2018
# Train PLIER models on SLE WB compendium with 10 random seeds. This will
# facilitate comparison with the subsampled recount2 data.
#
# USAGE: Rscript 16-repeat_sle_wb_PLIER.R

# functions
`%>%` <- dplyr::`%>%`
source(file.path("util", "plier_util.R"))

# set seed for reproducibility
set.seed(123)
seeds <- sample(1:10000, 10)

# create results directory
results.dir <- file.path("results", "16")
dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)

# SLE WB expression data
exprs.file <- file.path("data", "expression_data", 
                        "SLE_WB_all_microarray_QN_zto_before_with_GeneSymbol.pcl")
sle.exprs.df <- readr::read_tsv(exprs.file, progress = FALSE)
# matrix with gene symbol as rownames
exprs.mat <- dplyr::select(sle.exprs.df, -EntrezID)
rownames(exprs.mat) <- exprs.mat$GeneSymbol
exprs.mat <- as.matrix(dplyr::select(exprs.mat, -GeneSymbol))

# initialize list to hold models
model.list <- list()

# 10 repeats
for (seed in seeds) {
  
  # identifier for repeat
  base.file <- paste0("SLE-WB_PLIER_", seed)
  
  # run PLIER
  plier.model <- PLIERNewData(exprs.mat = exprs.mat,
                              seed = seed)
  
  # save model to file
  model.file <- file.path(results.dir, 
                          paste0(base.file, ".RDS"))
  saveRDS(object = plier.model, file = model.file)
  
  # add to model list
  model.list[[base.file]] <- plier.model
  
}

# evaluate models with wrapper function 
eval.list <- lapply(model.list, EvalWrapper)

# reshape list to data.frame for wrangling
eval.df <- reshape2::melt(eval.list)
colnames(eval.df) <- c("value", "pathway_coverage_type", "metric", "model")

# U sparsity -- we'll keep all and significant only in the same data.frame
sparsity.df <- eval.df %>%
  dplyr::filter(metric %in% c("all.sparsity", "sig.sparsity")) %>%
  dplyr::mutate(sparsity_type = metric) %>%
  dplyr::select(c(model, sparsity_type, value))

# number of lvs
num.lvs.df <- eval.df %>%
  dplyr::filter(metric == "num.lvs") %>%
  dplyr::mutate(num_lvs = value) %>%
  dplyr::select(c(model, num_lvs))

# pathway coverage
pathway.df <- eval.df %>%
  dplyr::filter(metric == "pathway.coverage") %>% 
  dplyr::select(c(model, pathway_coverage_type, value))

# write to file
sparsity.file <- file.path(results.dir, "sle-wb_repeated_sparsity.tsv")
readr::write_tsv(sparsity.df, sparsity.file)
num.file <- file.path(results.dir, "sle-wb_repeated_num_lvs.tsv")
readr::write_tsv(num.lvs.df, num.file)
pathway.file <- file.path(results.dir, "sle-wb_repeated_pathway.tsv")
readr::write_tsv(pathway.df, pathway.file)
