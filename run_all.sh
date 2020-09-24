#!/usr/bin/env bash

: '
this shell script was designed to automatically execute the python modules
sequentially instead of running them separately. Please modify where it is
says so in this script or leave it as is.
'

# You can modify the output directory name within each module but not here
file_path="output"
data_dir="data"
source_dir="source"

# [MODIFY] in case you want to change the default prefix
prefix="ReQTL_test"

echo "start executing python scripts sequentially ..."

# running the build_gen_exp_matrix.py with its command arguments
python $source_dir/build_gen_exp_matrix.py \
      -i $data_dir \
      -o $prefix
echo "Done with build_gen_exp_matrix.py"

# second python script build_VAF_matrix
python $source_dir/build_VAF_matrix.py \
      -r $data_dir \
      -o $prefix
echo "Done with build_VAF_matrix.py"

# NOTE: this script is optional but may help you avoid possible error or bug
python $source_dir/harmonize_matrices.py \
      -v $file_path/ReQTL_test_VAF_matrix.txt \
      -g $file_path/ReQTL_test_gene-exp_matrix.txt \
      -c $data_dir/covariates_matrix.txt
echo "Done with harmonize_matrices.py"

# for script 4, the main program

# 1: execute it by splitting cis and trans
python $source_dir/run_matrix_ReQTL.py \
      -s $file_path/ReQTL_test_VAF_matrix_harmonized.txt \
      -sl $file_path/ReQTL_test_VAF_loc_matrix.txt \
      -ge $file_path/ReQTL_test_gene-exp_matrix_harmonized.txt \
      -gl $file_path/ReQTL_test_gene-exp-loc_matrix.txt \
      -c $file_path/covariates_matrix_harmonized.txt \
      -ct T \
      -o $prefix \
      -pcis 0.001 \
      -ptra 0.00001
echo "Done with run_matrix_ReQTL with the split of cis and trans"

# 2: execute by unifying cis and trans
python $source_dir/run_matrix_ReQTL.py \
      -s $file_path/ReQTL_test_VAF_matrix_harmonized.txt \
      -sl $file_path/ReQTL_test_VAF_loc_matrix.txt \
      -ge $file_path/ReQTL_test_gene-exp_matrix_harmonized.txt \
      -gl $file_path/ReQTL_test_gene-exp-loc_matrix.txt \
      -c $file_path/covariates_matrix_harmonized.txt \
      -ct F \
      -o $prefix \
      -p 0.001
echo "Done with the unified cis and trans"

# run the annotation script
python $source_dir/annotate_cis_trans.py \
      -r $file_path/ReQTL_test_all_ReQTLs.txt \
      -ga $data_dir/gene_locations_hg38.txt \
      -o ReQTL_test
echo "Done with annotate_cis_trans"