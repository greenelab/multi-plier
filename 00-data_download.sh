#!/bin/bash

# set up directories
mkdir data && mkdir plots && mkdir results && mkdir util

# data directory subdirectories
cd data && mkdir expression_data

# get recount2 data & model from figshare, source code in 
# greenelab/rheum-data-plier
wget https://ndownloader.figshare.com/files/10881866 \
  -O recount2.zip
unzip recount2.zip && rm recount2.zip

## microarray data from greenelab/rheum-plier-data
cd expression_data

# sle-wb data
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/sle-wb/processed/aggregated_data/SLE_WB_all_microarray_QN_zto_before.pcl

# NARES
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/NARES/processed/NARES_SCANfast_ComBat.pcl

# GPA blood dataset (GSE18885)
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/gpa-blood/GSE18885_series_matrix.txt

# isolated blood cell populations from autoimmune conditions
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/isolated-cell-pop/processed/E-MTAB-2452_hugene11st_SCANfast.pcl

# glomeruli data
wget https://github.com/greenelab/rheum-plier-data/raw/55d86bb537f9e38c83fc3cca993cde48dc984411/glomeruli/ERCB_Glom_CustCDF19_forVCRC.txt

# get sample (e.g., phenotype) data
cd .. && mkdir sample_info && cd sample_info
# sle-wb sample to dataset of origin data
wget https://github.com/jaclyn-taroni/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/sle-wb/processed/sle-wb_sample_dataset_mapping.tsv
# other/single dataset sample information
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/sle-wb/arrayexpress/E-GEOD-65391/E-GEOD-65391.sdrf.txt
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/isolated-cell-pop/E-MTAB-2452.sdrf.txt
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/sle-wb/arrayexpress/E-GEOD-39088/E-GEOD-39088.sdrf.txt
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/sle-wb/arrayexpress/E-GEOD-78193/E-GEOD-78193.sdrf.txt
wget https://github.com/greenelab/rheum-plier-data/raw/4be547553f24fecac9e2f5c2b469a17f9df253f0/NARES/NARES_demographic_data.tsv
wget https://github.com/greenelab/rheum-plier-data/raw/55d86bb537f9e38c83fc3cca993cde48dc984411/glomeruli/ERCB_glom_diagnosis.tsv

