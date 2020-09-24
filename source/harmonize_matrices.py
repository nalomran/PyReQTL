#!/usr/bin/env python3

"""Builds harmonizes matrices

Created on Aug, 29 2020

@author: Nawaf Alomran

this module allows to harmonizes matrices so that the inputs needed for
running the module "run_matrix_ReQT.py" will contain the same samples.
Note: This is an optional module but it is recommended to run it in order to
avoid possible errors.

Input + Options
----------------

    + -v: the path to the variant matrix generated from the module
    build_VAF_matrix.py

    + -g: the path to the gene expression matrix generated from the module
    build_gen_exp_matrix.py

    + -c: the path to the covariate matrix located in the "data" folder


Return
-------
    + three output files encompasses the three input files which have
    mutual samples


How to Run
----------
      python harmonize_matrices.py -v output/ReQTL_test_VAF_matrix.txt \
      -g output/ReQTL_test_gene-exp_matrix.txt -c ../data/covariates_matrix.txt

* Python runtime with time command 1.19s user 0.24s system 305% cpu 0.466 total
* R time command line 1.33s user 0.13s system 96% cpu 1.504 total


"""

import argparse
import sys
from datetime import datetime

import pandas as pd  # type: ignore

from common import (create_output_dir, extract_filename_tag,
                    output_filename_generator)


def main(args) -> None:
    """Harmonizes matrices so the matrices can have the same samples.

    Parameter
    ---------
    args: arguments from the command line

        -v: the variant matrix generated from the script
        build_VAF_matrix.py

        -g: the gene expression matrix generated from the script
        build_gen_exp_matrix.py

        -c: the covariate matrix in the data folder

    Returns
    -------
    None

    Outputs
    -------
    - three files consist of the input files which have mutual samples

    """

    start_time = datetime.now()

    # load in matrices into pandas dataframe for each input file
    vaf_df = pd.read_table(args.var_dir, sep="\t")

    gene_express_df = pd.read_table(args.gex_dir, sep="\t")

    covar_df = pd.read_table(args.cov_dir, sep="\t")

    # extract samples columns for each dataframe
    vaf_cols = vaf_df.iloc[:, 1:].columns

    gene_express_cols = gene_express_df.iloc[:, 1:].columns

    covar_cols = covar_df.iloc[:, 1:].columns

    # get the samples in all 3 matrices via pandas intersection
    cols_intersect_tmp = vaf_cols.intersection(gene_express_cols)
    cols_intersect = list(cols_intersect_tmp.intersection(covar_cols))

    # select only the matched samples from all 3 samples
    vaf_df = vaf_df[["SNV"] + cols_intersect]

    gene_express_df = gene_express_df[["gene_id"] + cols_intersect]

    covar_df = covar_df[["id"] + cols_intersect]

    # create a directory "output" if it doesn't exists
    output = create_output_dir("output")

    # writes result to files
    filename_var = output_filename_generator(output,
                                             extract_filename_tag(
                                                 args.var_dir),
                                             '_harmonized.txt')

    vaf_df.to_csv(filename_var, sep="\t", index=False, na_rep='NA')

    filename_genexp = output_filename_generator(output,
                                                extract_filename_tag(
                                                    args.gex_dir),
                                                '_harmonized.txt')

    gene_express_df.to_csv(filename_genexp, sep="\t", index=False, na_rep='NA')

    filename_cov = output_filename_generator(output,
                                             extract_filename_tag(
                                                 args.cov_dir),
                                             '_harmonized.txt')

    covar_df.to_csv(filename_cov, sep="\t", index=False, na_rep='NA')

    print(f"\nfile for variant matrix is saved to {filename_var}")

    print(f"\nfile for gene expression matrix is saved to {filename_genexp}")

    print(f"\nfile is for gene covariate matrix is saved to {filename_cov}\n")

    print(f"Analysis took "
          f"{(datetime.now() - start_time).total_seconds()} sec")


def run_command_lines() -> None:
    """Parses the command line arguments entered by the user

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    parser = argparse.ArgumentParser(description="build harmonizes matrices."
                                                 "[OPTIONAL SCRIPT]!")
    parser.add_argument('-v',
                        dest='var_dir',
                        required=True,
                        help="""the path to the variant matrix generated from 
                        the module 'build_VAF_matrix.py'""")

    parser.add_argument('-g',
                        dest='gex_dir',
                        required=True,
                        help="""the path to the gene expression matrix 
                        generated from the module 'build_gen_exp_matrix.py'""")

    parser.add_argument('-c',
                        dest='cov_dir',
                        help="""the path to the covariate. In case this is not 
                        supplemented you have to remove the references to the 
                        covariate matrix in this very module""")

    args = parser.parse_args()

    try:
        main(args)
    except KeyboardInterrupt:
        sys.exit('\nthe user ends the program')


if __name__ == '__main__':
    run_command_lines()