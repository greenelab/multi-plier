#!/bin/bash

INPUT_DATA=data/recount2_PLIER_data/recount_data_prep_PLIER.RDS
MODEL_DIR=models

# sample size "grid"
declare -a sample_size_arr=("500" "1000" "2000" "4000" "8000" "16000" "32000")

for size in "${sample_size_arr[@]}"; 
do
  printf "\n\nTraining for sample size ${size}\n\n"
  output_file=${MODEL_DIR}/subsampled_recount2_PLIER_model_${size}.RDS  
  Rscript scripts/subsampling_PLIER.R \
   --input $INPUT_DATA \
   --output $output_file \
   --num_samples $size \
   --repeats 5
done