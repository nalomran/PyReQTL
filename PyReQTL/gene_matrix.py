#!/usr/bin/env python3

"""Transforming the raw expression files into matrix

Created on Aug, 29 2020

@author: Nawaf Alomran

This module transforms the raw gene expression files into matrix.

Inputs + Options
-----------------

    + -i: a directory contains all raw gene expression files from Stringtie

    + -o: the prefix for the gene expression and gene location files

    + -c: Whether the function is been executed with the command line. Default
    is False.


Return
------

In case cli argument is left unchanged "False" or set explicitly to False, the
function build_gene_exp_mat will return 2 dataframes for gene expression matrix and gene locus
information.


Output
-------
    + a file with the gene expression values to be fed into MatrixEQTL

    + a file with the gene location for MatrixEQTL



How to Run
-----------
    python -m PyReQTL.build_gen_exp_matrix -i data -o ReQTL_test -c True

* Python runtime with time command 2.31s user 0.39s system 293% cpu 0.921 total
* R time command line 5.06s user 0.49s system 94% cpu 5.858 total
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from scipy.stats import norm  # type: ignore

try:
    from common import create_output_dir, output_filename_generator
except ModuleNotFoundError:
    from PyReQTL.common import create_output_dir, output_filename_generator


def build_gene_matrix(gene_dir: str,
                      prefx_out: str = "ReQTL_test",
                      cli: bool = False):
    """Transforms the raw expression files into matrix

    Parameters
    ----------
    gene_dir: a directory consists of all the raw gene expression files from
            Stringtie

    prefx_out: the prefix for the gene expression and gene location files

    cli: if the function is been executed with the command line or not.
    Default is False.

    Return
    ------
    None

    Outputs
    -------
    - gene expression values matrix file to be fed in MatrixEQTL

    - gene location file for MatrixEQTL

    """
    # start time of analysis run
    start_time = datetime.now()

    # a suffix for the gene expression files name
    gene_exp_suffix = "_gene_abund.tab"

    # list of dataframes from the gene expression files
    dataframes = []

    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#
    # -----parsing and loading gene expression files into single dataframe----#
    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#

    # the code block below performs the following:
    #     1. reading & parsing all gene expression files
    #     2. selecting the necessary columns from each dataframe
    #     3. creating new column 'sample' denoting the sample name
    #     4. excluding those genes that lack annotations
    #     5. appending every gene expression dataframe into a list

    for filename in Path(gene_dir).glob("*" + gene_exp_suffix):
        gene_exps = pd.read_table(filename, sep="\t")

        select_cols = gene_exps[['Gene ID', 'Gene Name',
                                 'Reference', 'Start',
                                 'End', 'TPM']]

        clean_sample = select_cols.assign(sample=re.sub(gene_exp_suffix,
                                                        "",
                                                        filename.name))

        gene_exp_df = clean_sample[clean_sample['Gene Name'] != "-"]

        dataframes.append(gene_exp_df)

    # concatenate list of dataframes into one large dataframe
    gene_express_df = pd.concat(dataframes)

    # building gene coordinate table from the combined dataframe
    # define a new dataframe from certain columns from full_gene_exp_df df
    gene_loc = gene_express_df[['Gene ID', 'Reference', 'Start', 'End']]

    gene_loc = gene_loc.drop_duplicates()  # removing duplicates data

    # update the Reference column to include chr label right before its number
    gene_loc['Reference'] = gene_loc['Reference'].apply(lambda x: f"chr{x}")
    gene_loc.rename(columns={"Gene ID": "GeneID"}, inplace=True)

    gene_express_df = gene_express_df[['Gene ID', 'Gene Name',
                                       'TPM', 'sample']]

    # obtain the number of samples by measuring the length
    # of unique sample labels
    num_samples = len(gene_express_df['sample'].unique())

    # keep only genes that only expressed (TPM >= 1)
    # in at least 20% of the samples
    gene_express_df['count_zero'] = (gene_express_df['TPM'] < 1).groupby(
        gene_express_df['Gene ID']).transform('sum')

    gene_express_df['perc_zero'] = gene_express_df['count_zero'] / num_samples

    # remove duplicate entries for each gene if necessary
    gene_express_df = gene_express_df[gene_express_df['perc_zero'] < 0.8]
    gene_express_df = gene_express_df[['sample', 'TPM', 'Gene ID']]
    gene_express_df.groupby(['sample', 'Gene ID'])

    gen_expr_pivot = pd.pivot_table(gene_express_df,
                                    index='Gene ID',
                                    columns='sample',
                                    values="TPM")

    # reset the index to reverting to column-based
    gene_exp_revert_df = gen_expr_pivot.reset_index()

    # quantile normalize gene expression values code adopted from MatrixEQTL
    # sample code
    # http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/faq.html
    # making a replicated code from MatrixEQTL sample code but in python
    # with pandas which is much faster and efficient

    genes_exp_samples_df = gene_exp_revert_df.iloc[:, 1:]

    # rank the dataframe by column in the order they appear
    # for any values except na which will be assigned the highest rank
    order_genes_exp_df = genes_exp_samples_df.rank(axis=1,
                                                   method="first",
                                                   na_option="bottom")

    # handling the 0 values
    mask = (genes_exp_samples_df == 0)

    update_orders_df = np.broadcast_to(
        np.expand_dims(order_genes_exp_df[mask].mean(axis=1), 1),
        genes_exp_samples_df.shape)

    order_genes_exp_df.mask(mask, update_orders_df, inplace=True)

    matrix_computed = norm.ppf(
        order_genes_exp_df / (len(genes_exp_samples_df.columns) + 1))

    genes_exp_samples_df[:] = matrix_computed

    # adding a gene_id column to the genes_exp_samples_df data frame
    genes_exp_samples_df['gene_id'] = gene_exp_revert_df['Gene ID']

    # arrange the columns so that "gene_id" col will be placed to the
    # beginning of the dataframe
    cols = list(genes_exp_samples_df.columns)
    genes_exp_samples_df = genes_exp_samples_df[[cols[-1]] + cols[:-1]]

    # create a directory if it doesn't exists
    output_dir = create_output_dir("output")

    # write result to files
    filename_gen_exp = output_filename_generator(output_dir,
                                                 prefx_out,
                                                 "_gene-exp_matrix.txt")

    genes_exp_samples_df.to_csv(filename_gen_exp, sep="\t", index=False)

    file_gen_loc = output_filename_generator(output_dir,
                                             prefx_out,
                                             "_gene-exp-loc_matrix.txt")

    gene_loc.to_csv(file_gen_loc, sep="\t", index=False)

    print(f"\nfile is saved for gene expression in "
          f"{output_dir}/{prefx_out}_gene-exp_matrix.txt")

    print(f"file is saved for gene location in "
          f"{output_dir}/{prefx_out}_gene-exp-loc_matrix.txt\n")

    if cli:
        # calculate the time it took for the analysis from start to end
        print(f"Analysis took "
              f"{(datetime.now() - start_time).total_seconds()} sec")
    else:
        return genes_exp_samples_df, gene_loc


def main() -> None:
    """Parses the command line arguments entered by the user

    Parameters
    ----------
    None

    Return
    -------
    None

    """

    # parsing arguments from the command line with argparse.ArgumentParser()
    parser = argparse.ArgumentParser(
        description="build gene expression matrix")
    parser.add_argument("-i",
                        dest='gene_dir',
                        required=True,
                        help="""a directory with all raw gene expression files 
                        from Stringtie.[REQUIRED]!""")
    parser.add_argument("-o",
                        dest="prefx_out",
                        default="ReQTL_test",
                        help="""the prefix of the gene expression and gene 
                        location files. Default is 'ReQTL_test'!""")
    parser.add_argument("-c",
                        dest="cli",
                        default=False,
                        type=bool_conv_args,
                        help="""Whether the function is been executed with the 
                        command line. Default is False!""")

    args = parser.parse_args()
    gene_dir = args.gene_dir
    prefx_out = args.prefx_out
    cli = args.cli

    try:
        build_gene_matrix(gene_dir, prefx_out, cli)
    except KeyboardInterrupt:
        sys.exit('\nthe user ends the program')


if __name__ == '__main__':
    main()
