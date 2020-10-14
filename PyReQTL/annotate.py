#!/usr/bin/env python3

"""Annotate the output of ReQTL as cis or trans

Created on Aug, 29 2020

@author: Nawaf Alomran

This module annotates the output of ReQTL as cis or trans based on whether the
SNVs resides within its paired gene.

Input + Options
----------------

    + -r: the path to the ReQTL analysis result file

    + -ga: the path to the file gene location annotations

    + -o: the prefix for the output annotated result


Output
------
    + a file with the ReQTLs annotated as cis or trans


How to Run
----------
    python -m PyReQTL.annotate \
        -r output/ReQTL_test_all_ReQTLs.txt \
        -ga ../data/gene_locations_hg38.txt \
        -o ReQTL_test \
        -c True


* Python runtime via time command  8.19s user 0.61s system 112% cpu 7.838 total
* R time command line 3.15s user 0.22s system 99% cpu 3.383 total
* Note that the speed after the importr statements Python is faster than than R

"""

import argparse
import sys
from datetime import datetime

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import rpy2.robjects.packages as rpackages  # type: ignore
from rpy2.robjects import pandas2ri  # type: ignore
from rpy2.robjects.packages import importr  # type: ignore

try:
    from common import (create_output_dir, output_filename_generator,
                        bool_conv_args)
except ModuleNotFoundError:
    from PyReQTL.common import (create_output_dir, output_filename_generator,
                                bool_conv_args)

# install the R package GenomicFeatures from within Python
if not rpackages.isinstalled('GenomicFeatures'):
    print("installing GenomicFeatures package ...")
    bioc_manager = rpackages.importr('BiocManager')
    bioc_manager.install('GenomicFeatures')
    print("Done installing the package.")

# importing the following required R packages to be used within Python
print("Kindly wait for the required R packages to be imported into Python...")

g_ranges = importr('GenomicRanges')
print("GenomicRanges package is imported.")

g_alignments = importr('GenomicAlignments')
print("GenomicAlignments package is imported.")

iranges = importr('IRanges')
print("IRanges package is imported.")

print("Done importing.")

# This needs to be activated in order to perform pandas conversion
pandas2ri.activate()


def cis_trans_annotator(rqt_rst: str,
                        gene_ann: str,
                        out_prefx: str,
                        cli: bool = False) -> None:
    """Annotate the output of ReQTL as cis or trans based on whether the
    SNVs resides within its paired gene

    Parameter
    ---------
    rqt_rst: the path to the ReQTL analysis result file

    gene_ann: the path to the file gene location annotation

    out_prefx: the prefix for the output annotated result

    cli: Whether the function is been executed with the command line.
    Default is False.

    Return
    ------
    reqtl_reslt_arranged: dataframe ReQTLs annotated as cis or trans

    Output
    ------
    - file with the ReQTLs annotated as cis or trans

    """

    start_time = datetime.now()

    # reading the ReQTL result file from run_matrix_ReQTL
    reqtl_result = pd.read_table(rqt_rst, sep="\t")

    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#
    # -----------------annotate which gene harbors the snp--------------------#
    # classify ReQTLs in which the two members of the pair are in the same----#
    # gene as cis and classify all others as trans----------------------------#
    # ------------------------------------------------------------------------#
    # ------------------------------------------------------------------------#

    reqtl_reslt_arranged = reqtl_result.assign(new_SNP=reqtl_result.SNP)

    # split them into four columns based on the pattern "[:_>]"
    reqtl_reslt_arranged = reqtl_reslt_arranged.new_SNP.str.split('[:_>]',
                                                                  expand=True)
    reqtl_reslt_arranged.columns = ['chrom', 'start', 'ref', 'alt']

    # concatenating the re-arranged dataframe with the original dataframe
    reqtl_reslt_arranged = pd.concat([reqtl_result, reqtl_reslt_arranged],
                                     axis=1)
    # making the new end column the same as the start column
    reqtl_reslt_arranged = reqtl_reslt_arranged.assign(
        end=reqtl_reslt_arranged.start)

    # convert Python Pandas DataFrame to R-dataframe
    reqtl_result_df_r = pandas2ri.py2rpy(reqtl_reslt_arranged)

    # read gene location file and then convert to R dataframe
    gene_locs_py_df = pd.read_table(gene_ann, sep="\t")
    gene_locs_df_r = pandas2ri.py2rpy(gene_locs_py_df)

    # storing the location of genomic features for both R dataframes
    reqtl_reslt_granges_r = g_ranges.GRanges(reqtl_result_df_r)
    gene_loc_granges_r = g_ranges.GRanges(gene_locs_df_r)

    # finding the overlap between the ranges
    overlaps = iranges.findOverlaps(reqtl_reslt_granges_r,
                                    gene_loc_granges_r,
                                    select="last",
                                    type="within")

    # ignore the Pycharm warning later
    overlaps = np.where(overlaps == -2147483648, None, overlaps)
    overlaps = overlaps.tolist()

    # reindex the gene_locs dataframe by the overlaps
    genes_snp = gene_locs_py_df.ensembl_gene.reindex(overlaps)

    reqtl_reslt_arranged['genes_snp'] = pd.Series(genes_snp.values.tolist())

    # if genes_snp == gene in reqtl_reslt_arranged dataframe then it cis
    # otherwise it will be trans
    reqtl_reslt_arranged['class'] = np.where(
        reqtl_reslt_arranged.genes_snp == reqtl_reslt_arranged.gene,
        'cis',
        'trans')

    reqtl_reslt_arranged.loc[reqtl_reslt_arranged['genes_snp'].isna(),
                             'class'] = reqtl_reslt_arranged['genes_snp']

    # drop the unneeded columns
    reqtl_reslt_arranged.drop(
        ['chrom',
         'end',
         'ref',
         'alt',
         'start'], axis=1, inplace=True)

    out_dir = create_output_dir("output")

    annotated_file = output_filename_generator(out_dir,
                                               out_prefx,
                                               "_ReQTLs_cistrans_ann.txt")

    reqtl_reslt_arranged.to_csv(annotated_file, sep="\t", index=False,
                                na_rep='NULL')

    print(f"\nCis/trans annotated ReQTLs saved in {annotated_file}\n")

    if cli:
        print(f"Analysis took after importing the required packages "
              f"{(datetime.now() - start_time).total_seconds()} sec")
    else:
        return reqtl_reslt_arranged


def main() -> None:
    """Parses the command line arguments entered by the user

    Parameters
    ---------
    None

    Return
    -------
    None

    """

    USAGE = """Annotate the output of ReQTL as cis or trans based on whether 
    the SNV resides within its paired gene"""

    parser = argparse.ArgumentParser(description=USAGE)
    parser.add_argument('-r',
                        dest="rqt_rst",
                        required=True,
                        help="the path to the ReQTL analysis result file")
    parser.add_argument('-ga',
                        dest='gene_ann',
                        required=True,
                        help="the path to the file gene location annotations")
    parser.add_argument('-o',
                        dest="out_prefx",
                        required=True,
                        help="the prefix for the output annotated result")
    parser.add_argument("-c",
                        dest="cli",
                        default=False,
                        type=bool_conv_args,
                        help="""Whether the function is been executed with the 
                        command line. Default is False!""")

    args = parser.parse_args()

    rqt_rst = args.rqt_rst
    gene_ann = args.gene_ann
    out_prefx = args.out_prefx
    cli = args.cli

    try:
        cis_trans_annotator(rqt_rst, gene_ann, out_prefx, cli)
    except KeyboardInterrupt:
        sys.exit('\nthe user ends the program')


if __name__ == '__main__':
    main()
