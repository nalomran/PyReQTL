#!/usr/bin/env python3

"""Transforms the read counts into a variant fraction matrix

Created on Aug, 29 2020

@author: Nawaf Alomran

This module transforms the read counts into a variant fraction matrix.


Input + Options
----------------

    + -r: a directory containing the ".csv" files from the output of readCounts
    python module

    + -o: the prefix for the SNV matrix and SNV location files


Output
------
    + a file containing the SNV locations for MatrixEQTL

    + a file with the SNV variant allele fraction matrix for
    MatrixEQTL


How to Run
-----------
    python -m PyReQTL.vaf_matrix -r data -o ReQTL_test -c True

* Python runtime with time command 1.32s user 0.23s system 254% cpu 0.610 total
* R time command line 1.65s user 0.18s system 91% cpu 2.011 total

"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
from numpy import savetxt  # type: ignore

try:
    from common import (create_output_dir, output_filename_generator,
                        bool_conv_args, remove_vars)
except ModuleNotFoundError:
    from PyReQTL.common import (create_output_dir, output_filename_generator,
                                bool_conv_args, remove_vars)


def build_vaf_matrix(read_dir: str,
                     prefx_out: str = "ReQTL_test",
                     cli: bool = False):

    """Transforms the read counts into a variant fraction matrix

    Parameters
    ----------
    read_dir: a directory with the .csv files from the output of readCounts

    prefx_out: the prefix for the SNV matrix and SNV location files

    cli: Whether the function is been executed with the command line.
    Default is False!

    Returns
    -------
    In case cli argument is left unchanged "False" or set explicitly to False,
    the following will be retuned:

    - read_counts_df: numpy array with the SNV variant allele fraction matrix
    for MatrixEQTL

    - snv_loc: numpy array with the SNV locations for MatrixEQTL

    Output
    ------
    - SNV locations matrix file for MatrixEQTL

    - SNV variant allele fraction matrix file for MatrixEQTL

    """

    start_time = datetime.now()

    # list of dataframes from the readCounts module files
    dataframes = []

    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#
    # -------------------parsing and loading readCounts files-----------------#
    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#

    # the code block below performs the following:
    #     1. reading & parsing all readCounts files
    #     2. selecting the necessary columns from each dataframe
    #     3. appending every gene expression dataframe into a list
    # ------------------------------------------------------------------------#
    #        parsing and loading readCount files into single dataframe        #
    # ------------------------------------------------------------------------#

    for filename in Path(read_dir).glob("*csv"):
        read_counts = pd.read_table(filename, sep=",")
        select_cols = read_counts[['CHROM', 'POS',
                                   'REF', 'ALT',
                                   'AlignedReads', 'R']]
        dataframes.append(select_cols)

    read_counts_df = pd.concat(dataframes)

    # number of unique sample labels
    num_samples = len(read_counts_df['AlignedReads'].unique())

    # assign a new SNV coordinate column to serve as SNV identifier
    read_counts_df['SNV'] = read_counts_df['CHROM'].astype(str) + ":" + \
                            read_counts_df['POS'].astype(str) + "_" + \
                            read_counts_df['REF'].astype(str) + ">" + \
                            read_counts_df['ALT'].astype(str)

    # remove variants that are homozygous var in more that 80% of the samples
    cols_var = ["count_homo_var", 'perc_homo_var']
    read_counts_df = remove_vars(read_counts_df, 1, num_samples, cols_var)

    # remove variants that are homozygous ref in more that 80% of the samples
    cols_ref = ["count_homo_ref", 'perc_homo_ref']
    read_counts_df = remove_vars(read_counts_df, 0, num_samples, cols_ref)

    # remove variants that are NA in more than 80% of the samples
    read_counts_df['count'] = read_counts_df.SNV. \
        groupby(read_counts_df['SNV']).transform('count')
    read_counts_df['perc_non_na'] = read_counts_df['count'] / num_samples

    read_counts_df = read_counts_df[read_counts_df['perc_non_na'] > 0.2]

    # select relevant columns and create SNV location matrix
    snv_loc = read_counts_df[['SNV', 'CHROM', 'POS']].copy()
    snv_loc['CHROM'] = snv_loc['CHROM'].apply(lambda x: f"chr{x}")
    snv_loc = np.unique(snv_loc.to_records(index=False))

    # select only the following columns for the dataframe "read_counts_df"
    read_counts_df = read_counts_df[["SNV", "AlignedReads", 'R']]

    # building the SNV matrix with pivot
    read_counts_df = read_counts_df.pivot(index="SNV",
                                          columns="AlignedReads",
                                          values="R")

    # SNV column first then the rest of the columns
    header = ["SNV"] + list(read_counts_df.columns)
    header_loc = ['SNV', 'CHROM', 'POS']

    # converting from dataframe to matrix (numpy array)
    read_counts_df = read_counts_df.to_records()

    # create a directory if it doesn't exists
    output = create_output_dir("output")

    # format the output files based on their columns types including nan types
    dig_fmt = '%f,' * num_samples
    fmt_str = '%s,' + dig_fmt[:-1]  # string + the rest of the columns to float
    fmt = fmt_str.split(",")  # generating a list

    # format column for SNV locations file
    fmt_loc = "%s,%s,%.0f".split(",")

    snv_mat_file = output_filename_generator(output,
                                             prefx_out,
                                             '_VAF_matrix.txt')

    # added comments to remove "#" at the beginning of the header
    savetxt(snv_mat_file,
            read_counts_df,
            delimiter="\t",
            comments='',
            fmt="\t".join(fmt),
            header=str("\t".join([str(hdr) for hdr in header])))

    snv_loc_file = output_filename_generator(output,
                                             prefx_out,
                                             '_VAF_loc_matrix.txt')

    savetxt(snv_loc_file,
            snv_loc,
            delimiter="\t",
            comments='',
            fmt="\t".join(fmt_loc),
            header=str("\t".join([str(hdr) for hdr in header_loc])))

    print(f"\nfile is saved for variant matrix to {output}/{prefx_out}"
          f"_VAF_matrix.txt")

    print(f"\nfile is saved for gene location to {output}/{prefx_out}"
          f"_VAF_loc_matrix.txt\n")

    if cli:
        print(f"Analysis took "
              f"{(datetime.now() - start_time).total_seconds()} sec")
    else:
        return read_counts_df, snv_loc


def run_command_lines() -> None:
    """Parses the command line arguments entered by the user

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    parser = argparse.ArgumentParser(description="build variant allele "
                                                 "fraction expression matrix")
    parser.add_argument("-r",
                        dest="read_dir",
                        required=True,
                        help="""a directory with the .csv files from the 
                        output of readCounts module. [REQUIRED]!""")
    parser.add_argument("-o",
                        dest="prefx_out",
                        default="ReQTL_test",
                        help="""the prefix for the SNV matrix and SNV location 
                        files. [OPTIONAL]!""")
    parser.add_argument("-c",
                        dest="cli",
                        default=False,
                        type=bool_conv_args,
                        help="""Whether the function is been executed with the 
                        command line. Default is False!""")

    args = parser.parse_args()

    read_dir = args.read_dir
    prefx_out = args.prefx_out
    cli = args.cli

    try:
        build_vaf_matrix(read_dir, prefx_out, cli)
    except KeyboardInterrupt:
        sys.exit('\nthe user ends the program')


if __name__ == '__main__':
    run_command_lines()
