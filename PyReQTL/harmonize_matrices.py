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
    vaf_matrix.py

    + -g: the path to the gene expression matrix generated from the module
    gene_matrix.py

    + -c: the path to the covariate matrix located in the "data" folder


Return
-------
    + three output files encompasses the three input files which have
    mutual samples


How to Run
----------
      python -m PyReQTL.harmonize_matrices \
      -v output/ReQTL_test_VAF_matrix.txt \
      -g output/ReQTL_test_gene-exp_matrix.txt \
      -cov data/covariates_matrix.txt \
      -c True

* Python runtime with time command 1.19s user 0.24s system 305% cpu 0.466 total
* R time command line 1.33s user 0.13s system 96% cpu 1.504 total


"""

import argparse
import sys
from datetime import datetime

import pandas as pd  # type: ignore

try:
    from common import (create_output_dir, output_filename_generator,
                        bool_conv_args, extract_filename_tag)
except ModuleNotFoundError:
    from PyReQTL.common import (create_output_dir, output_filename_generator,
                                bool_conv_args, extract_filename_tag)


def harmonize_matrices(var_dir: str,
                       gex_dir: str,
                       cov_dir: str,
                       cli: bool = False):
    """Harmonizes matrices so the matrices can have the same samples.

    Parameter
    ---------
    var_dir: the variant matrix generated from the script vaf_matrix.py

    gex_dir: the gene expression matrix generated from the script
    gene_matrix.py

    cov_dir: the covariate matrix in the data folder

    cli: Whether the function is been executed with the command line.
    Default is False.

    Return
    ------
    In case cli argument is left unchanged "False" or set explicitly to False,
    the following will be returned:

    - vaf_df: a harmonized variant dataframe

    - gene_express_df: a harmonized gene expression dataframe

    - covar_df: a harmonized covariate matrix dataframe


    Outputs
    -------
    - three files consist of the input files which have mutual samples

    """

    start_time = datetime.now()

    # load in matrices into pandas dataframe for each input file
    vaf_df = pd.read_table(var_dir, sep="\t")

    gene_express_df = pd.read_table(gex_dir, sep="\t")

    covar_df = pd.read_table(cov_dir, sep="\t")

    # extract samples columns for each dataframe
    vaf_cols = vaf_df.iloc[:, 1:].columns

    gene_express_cols = gene_express_df.iloc[:, 1:].columns

    covar_cols = covar_df.iloc[:, 1:].columns

    # get the samples in all 3 matrices with pandas intersection
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
                                             extract_filename_tag(var_dir),
                                             '_harmonized.txt')

    vaf_df.to_csv(filename_var, sep="\t", index=False, na_rep='NA')

    filename_genexp = output_filename_generator(output,
                                                extract_filename_tag(gex_dir),
                                                '_harmonized.txt')

    gene_express_df.to_csv(filename_genexp, sep="\t", index=False, na_rep='NA')

    filename_cov = output_filename_generator(output,
                                             extract_filename_tag(cov_dir),
                                             '_harmonized.txt')

    covar_df.to_csv(filename_cov, sep="\t", index=False, na_rep='NA')

    print(f"\nfile for variant matrix is saved to {filename_var}")

    print(f"\nfile for gene expression matrix is saved to {filename_genexp}")

    print(f"\nfile is for gene covariate matrix is saved to {filename_cov}\n")

    if cli:
        print(f"Analysis took "
              f"{(datetime.now() - start_time).total_seconds()} sec")

    else:
        return vaf_df, gene_express_df, covar_df


def main() -> None:
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
                        the module 'vaf_matrix.py'""")
    parser.add_argument('-g',
                        dest='gex_dir',
                        required=True,
                        help="""the path to the gene expression matrix 
                        generated from the module 'gene_matrix.py'""")
    parser.add_argument('-cov',
                        dest='cov_dir',
                        help="""the path to the covariate. In case this is not 
                        supplemented you have to remove the references to the 
                        covariate matrix in this very module""")
    parser.add_argument("-c",
                        dest="cli",
                        default=False,
                        type=bool_conv_args,
                        help="""Whether the function is been executed with the 
                        command line. Default is False!""")

    args = parser.parse_args()
    var_dir = args.var_dir
    gex_dir = args.gex_dir
    cov_dir = args.cov_dir
    cli = args.cli

    try:
        harmonize_matrices(var_dir, gex_dir, cov_dir, cli)
    except KeyboardInterrupt:
        sys.exit('\nthe user ends the program')


if __name__ == '__main__':
    main()
