#!/usr/bin/env python3

"""Runs the ReQTL analysis using MatrixEQTL

Created on Aug, 29 2020

@author: Nawaf Alomran

This module is based off the sample code from Shabalin, et al (2012) which is
an R package "designed for fast eQTL analysis on large datasets that test for
association between genotype and gene expression using linear regression
including ANOVA genotype effects". For more information about the package,
please consider visiting the package's page at:
http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

|--------------|
|Important Note|
|--------------|

Due to the lack of an equivalent python library to R package "MatrixEQTL" and
to my knowledge I believe that one good alternative to overcome this is by
using rpy2 library to interface with R codes, objects or even packages within
Python. It is noteworthy to cite the documentation of the rpy2 library that
rpy2 is "more efficient, better integrated with Python" than using subprocess.
More information is found in:
https://rpy2.github.io/doc/latest/html/introduction.html

Inputs + Options
-----------------

    + -s: the SNV or variant matrix file created from harmonize_matrices

    + -sl: the SNV location matrix file created from build_VAF_matrix

    + -ge: the gene expression matrix file created from build_gene-exp_matrix

    + -gl: the gene locations file created from build_gene-exp_matrix

    + -c: the covariates matrix file created from "harmonize_matrices".
    [OPTIONAL]! you can also get the file under data if you wish

    + -o: the prefix for the path to the output files

    + -ct: logical (T or F) specifying whether to split the output into cis
    or trans

    + -pcis: p-value thresholds for the cis output files

    + -ptran: p-value thresholds for the trans output files

    + -p: p-value thresholds for the unified output file


Output
-------
    + one output could be cis ReQTLs and trans ReQTLs or one with all of the
    unified ReQTLs. This depends on what you choose for the value of parameter
    "-ct"

     + one QQ plot of p-values


How to Run
----------
     # execute it by splitting cis and trans
     python run_matrix_ReQTL.py \
      -s output/ReQTL_test_VAF_matrix_harmonized.txt \
      -sl output/ReQTL_test_VAF_loc_matrix.txt \
      -ge output/ReQTL_test_gene-exp_matrix_harmonized.txt \
      -gl output/ReQTL_test_gene-exp-loc_matrix.txt \
      -c output/covariates_matrix_harmonized.txt \
      -ct T \
      -o "ReQTL_test" \
      -pcis 0.001 \
      -ptra 0.00001

     # execute by unified cis and trans
     python run_matrix_ReQTL.py \
       -s output/ReQTL_test_VAF_matrix_harmonized.txt \
       -sl output/ReQTL_test_VAF_loc_matrix.txt \
       -ge output/ReQTL_test_gene-exp_matrix_harmonized.txt \
       -gl output/ReQTL_test_gene-exp-loc_matrix.txt \
       -c output/covariates_matrix_harmonized.txt \
       -ct F \
       -o "ReQTL_test" \
       -p 0.001

* Python runtime (T) with time 2.41s user 0.36s system 146% cpu 1.890 total
* R time command line 2.09s user 0.22s system 85% cpu 2.695 total

* Python runtime of (F) via time 2.24s user 0.37s system 151% cpu 1.730 total
* R time command line 1.70s user 0.18s system 88% cpu 2.131 total

"""

import argparse
import sys
from datetime import datetime

import rpy2.robjects as robjects  # type: ignore
import rpy2.robjects as ro
from rpy2.robjects.packages import importr  # type: ignore

from common import create_output_dir, output_filename_generator

# use the following R operators to get and set R attributes
get_r_attribute = ro.baseenv['$']
set_r_attribute = ro.baseenv['$<-']


class MapTOS4(ro.methods.RS4):
    """Mapping SR4 class to Python class which will extend rpy2’s RS4.

    This class will allow to access attributes or fields and method of
    SlicedData class.

    """

    def __init__(self, r_obj, file_slice_size=2000):
        super().__init__(r_obj)
        self.file_slice_size = file_slice_size

    def load_file(self, filename):
        """Access the LoadFile method of SlicedData class

        Parameters
        ----------
        filename: the name of file to be loaded into SlicedData class

        Return
        -------
        get the R method which load the filename

        """

        return get_r_attribute(self, 'LoadFile')(filename)

    @property
    def file_slice_size(self):
        """Access the fileSliceSize field or attribute of SlicedData class

        Parameters
        ----------
        None

        Return
        -------
        get the R attribute fileSliceSize

        """
        return get_r_attribute(self, 'fileSliceSize')

    @file_slice_size.setter
    def file_slice_size(self, value):
        """Access the fileSliceSize field or attribute of SlicedData class

        Parameters
        ----------
        value: the value to be set for fileSliceSize field

        Return
        -------
        None

        """

        set_r_attribute(self, 'fileSliceSize', value)


def main(args) -> None:
    """This function will be based off the sample code from Shabalin,
    et al (2012) of the R package MatrixEQTL.

    Parameters
    ----------
    args: please read the above docstring (comments at the beginning of the
    module) for more information about the arguments used.

    Returns
    -------
    None


    Output
    ------
    - file either cis ReQTLs and trans ReQTLs or with all of the unified
    ReQTLs.

    - file QQ plot of p-values

    """

    start_time = datetime.now()

    # check for installed package or install it, installing MatrixEQTL
    r_str_download = """
          testPkg <- function(x){
              if (!require(x,character.only = TRUE))
              {
                install.packages("MatrixEQTL",dep=TRUE)
                 if(!require(x,character.only = TRUE)) stop("missing package!")
              }
            }
          testPkg('MatrixEQTL')
      """

    robjects.r(r_str_download)

    # import MatrixEQTL package
    mql = importr("MatrixEQTL")
    # import utils package
    utils = importr("utils")
    # import grDevices package
    gr_devices = importr('grDevices')
    # import base package
    base = importr('base')

    snv_filename = args.snv
    snvs_data = MapTOS4(mql.SlicedData())
    # load snv/genotype data into SlicedData class
    snvs_data.load_file(snv_filename)

    # load gene expression data file into SlicedData class
    gene_express_filename = args.gen_exp
    gene_exp_data = MapTOS4(mql.SlicedData(), file_slice_size=2000)
    gene_exp_data.load_file(gene_express_filename)

    # load Covariates data
    # covar_filename = args.cov_mt
    covar_data = MapTOS4(mql.SlicedData(), file_slice_size=1000)

    snv_loc_filename = args.snv_loc

    # need the utils package to read table for the downstream analysis
    snv_pos = utils.read_table(snv_loc_filename, header=True,
                               stringsAsFactors=False)

    gene_loc_filename = args.gen_loc
    gene_pos = utils.read_table(gene_loc_filename, header=True,
                                stringsAsFactors=False)

    # value of either C for cis and T for trans
    cis_or_trans = args.ct

    output = create_output_dir("output")

    output_trans_file = output_filename_generator(output,
                                                  args.out_dir,
                                                  "_trans_ReQTLs.txt")

    output_cis_file = output_filename_generator(output,
                                                args.out_dir,
                                                "_cis_ReQTLs.txt")

    # call the matrix_eQTL_main of MatrixEQTL package in case of trans case
    if cis_or_trans == "T":

        mat_eqtl = mql.Matrix_eQTL_main(
            snps=snvs_data,
            gene=gene_exp_data,
            cvrt=covar_data,
            output_file_name=output_trans_file,
            pvOutputThreshold=float(args.ptra),
            useModel=117348,
            verbose=False,
            output_file_name_cis=output_cis_file,
            pvOutputThreshold_cis=float(args.pcis),
            snpspos=snv_pos,
            genepos=gene_pos,
            cisDist=1e6,
            pvalue_hist="qqplot",
            min_pv_by_genesnp=False,
            noFDRsaveMemory=False
        )

    else:

        mat_eqtl = mql.Matrix_eQTL_main(
            snps=snvs_data,
            gene=gene_exp_data,
            output_file_name=output_trans_file,
            useModel=117348,
            verbose=False,
            pvOutputThreshold=float(args.ptra),
            snpspos=snv_pos,
            genepos=gene_pos,
            pvalue_hist="qqplot",
            min_pv_by_genesnp=False,
            noFDRsaveMemory=False
        )

    ggplot_file = output_filename_generator(output, args.out_dir,
                                            "_qqplot.tiff")

    gr_devices.tiff(filename=ggplot_file)

    base.plot(mat_eqtl)

    gr_devices.dev_off()

    print(f"Analysis took {(datetime.now() - start_time).total_seconds()} sec")


def run_command_lines() -> None:
    """Parses the command line arguments entered by the user

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    USAGE = """Runs the ReQTL analysis using MatrixEQTL package"""

    parser = argparse.ArgumentParser(description=USAGE)
    parser.add_argument('-s',
                        required=True,
                        dest='snv',
                        help="the SNV or variant matrix file from "
                             "harmonize_matrices")

    parser.add_argument('-sl',
                        required=True,
                        dest='snv_loc',
                        help="the SNV location matrix file from build_"
                             "VAF_matrix")

    parser.add_argument('-ge',
                        required=True,
                        dest="gen_exp",
                        help="gene expression file matrix from "
                             "build_gene-exp_matrix")

    parser.add_argument('-gl',
                        required=True,
                        dest="gen_loc",
                        help="gene locations file from build_gene-exp_matrix")

    parser.add_argument('-c',
                        dest="cov_mt",
                        help="""the covariates matrix file from "
                             "harmonize_matrices. [OPTIONAL]!""")

    parser.add_argument('-o',
                        dest="out_dir",
                        required=True,
                        help="the prefix for the path to the output files")

    parser.add_argument('-ct',
                        required=True,
                        help="logical (T or F) specifying whether to split "
                             "the output into cis or trans")

    parser.add_argument('-pcis',
                        help="p-value thresholds for the cis output files")

    parser.add_argument('-ptra',
                        help="p-value thresholds for the cis")

    parser.add_argument('-p',
                        help="p-value thresholds for the unified output file")

    args = parser.parse_args()

    try:
        main(args)
    except KeyboardInterrupt:
        sys.exit('\nthe user ends the program')


if __name__ == '__main__':
    run_command_lines()
