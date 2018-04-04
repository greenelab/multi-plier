# J. Taroni
# Wrapper script for knitr::purl
# .Rmd -> .R for easy code review
# USAGE: Rscript util/purl_wrapper.R <R MARKDOWN FILE>

args <- commandArgs(trailingOnly = TRUE)
input.rmd <- args[1]

# create directory for notebook-associated scripts, don't warn if it 
# already exists
dir.create("Rnotebook_scripts", showWarnings = FALSE)

# retain the same name for the file
output.file <- sub(".Rmd", ".R", input.rmd)

# purl
knitr::purl(input.rmd, 
            output = file.path("Rnotebook_scripts", output.file),
            documentation = 2)
