#!/bin/bash

INPUT_DATA=data/recount2_PLIER_data/recount_data_prep_PLIER.RDS
MODEL_DIR=models
SAMPLE_LIST_DIR=data/sample_info

# all files that end in accessions.tsv -- these are the lists of samples that
# are from particular biological contexts (e.g., cancer)
arr=($SAMPLE_LIST_DIR/*accessions.tsv)

# input data will always be the sample in this case

for f in "${arr[@]}"; 
do
  
  echo "Training models for ${f}"
  
  # the output file should be the same file name as the sample list file, but
  # we're removing .tsv extension and replacing with "_PLIER_model.RDS" and 
  # it needs to be saved in the model directory
  output_file=${f/\.tsv/_PLIER_model\.RDS}
  output_file=${MODEL_DIR}/${output_file/$SAMPLE_LIST_DIR/}
  
  # run the Rscript
  Rscript scripts/subsampling_PLIER.R \
    --input $INPUT_DATA \
    --output $output_file \
    --use_sample_list \
    --sample_list $f

done
