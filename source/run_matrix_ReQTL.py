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
from rpy2.robjects import globalenv  # type: ignore
from rpy2.robjects.packages import importr  # type: ignore

from common import create_output_dir, output_filename_generator


def load_data_from_MatrixEQTL(input_file: str):
    """Load data including snvs and gene expression
    using MatrixEQTL R package within Python using rpy2.robjects.r.

    Parameters
    -----------
    input_file: snv and gene expression


    Returns
    -------
    load_data: a RS4 object which is a dataframe that is stored in slice object

    """

    # create a variable name in R's global environment for each input file
    globalenv['input_file'] = input_file

    r_string_load_data = """
        load_data <- SlicedData$new()
        load_data$fileDelimiter = "\t"
        load_data$fileOmitCharacters = "NA"
        load_data$fileSkipRows = 1
        load_data$fileSkipColumns = 1
        load_data$fileSliceSize = 2000
        load_data$LoadFile(input_file)
    """

    load_data = robjects.r(r_string_load_data)

    return load_data


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

    # load the library MatrixEQTL in R from within Python and load it
    # into a namespace
    robjects.r("library('MatrixEQTL')")

    # import utils package
    utils = importr("utils")

    # load snv/genotype data from MatrixEQTL using rpy2.robjects.r
    # and then store it to a variable in R's global environment
    snv_filename = args.snv

    # load_data_from_MatrixEQTL does the magic show R-to-Python interfacing
    snvs = load_data_from_MatrixEQTL(snv_filename)
    globalenv['snvs'] = snvs

    # load gene expression data file same the was as the snv data
    gene_express_filename = args.gen_exp
    gene = load_data_from_MatrixEQTL(gene_express_filename)
    globalenv['gene'] = gene

    # load Covariates data
    covar_filename = args.cov_mt
    globalenv['covar_filename'] = covar_filename

    r_string_load_covar = '''
            cvrt <- SlicedData$new()
            cvrt$fileDelimiter = "\t"      
            cvrt$fileOmitCharacters = "NA" 
            cvrt$fileSkipRows = 1          
            cvrt$fileSkipColumns = 1       
        
            if(length(covar_filename) > 1) {
                cvrt$LoadFile(covar_filename)
            } else{cvrt}
    '''

    cov_r_obj = robjects.r(r_string_load_covar)
    globalenv['cvrt'] = cov_r_obj

    snv_loc_filename = args.snv_loc

    # need the utils package to read table for the downstream analysis

    snv_pos = utils.read_table(snv_loc_filename, header=True,
                               stringsAsFactors=False)

    gene_loc_filename = args.gen_loc

    gene_pos = utils.read_table(gene_loc_filename, header=True,
                                stringsAsFactors=False)

    # value of either C for cis and T for trans
    cis_or_trans = args.ct

    # globalenv['cis_or_trans'] = cis_or_trans

    output = create_output_dir("output")

    globalenv['output_trans_file'] = output_filename_generator(output,
                                                               args.out_dir,
                                                               "_trans_ReQTLs.txt")

    globalenv['output_cis_file'] = output_filename_generator(output,
                                                             args.out_dir,
                                                             "_cis_ReQTLs.txt")

    globalenv['snv_pos'] = snv_pos
    globalenv['gene_pos'] = gene_pos

    globalenv['output_file_name'] = output_filename_generator(output,
                                                              args.out_dir,
                                                              "_all_ReQTLs.txt")

    # call the matrix_eQTL_main of MatrixEQTL package in case of trans case

    # code inherited from R ReQTL toolkit
    r_str_T_option = '''
        useModel <<- modelLINEAR
        errorCovariance <<- numeric()
        cisDist <<- 1e6
        
        Mat_eQTL_m <- Matrix_eQTL_main(
            snps = snvs,
            gene = gene,
            cvrt = cvrt,
            output_file_name = output_trans_file,
            pvOutputThreshold  = as.numeric(pval_tra),
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            output_file_name.cis = output_cis_file,
            pvOutputThreshold.cis = as.numeric(pval_cis),
            snpspos = snv_pos,
            genepos = gene_pos,
            cisDist = cisDist,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE
        )
    '''

    r_str_F_option = '''
        useModel <<- modelLINEAR
        errorCovariance <<- numeric()
    
        Mat_eQTL_m <- Matrix_eQTL_main(
        snps = snvs,
        gene = gene,
        output_file_name = output_file_name,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = FALSE,
        pvOutputThreshold= as.numeric(pval),
        snpspos = snv_pos,
        genepos = gene_pos,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    );
    
    '''

    if cis_or_trans == "T":
        globalenv['pval_cis'] = args.pcis
        globalenv['pval_tra'] = args.ptra
        meql_main = robjects.r(r_str_T_option)
        globalenv['Mat_eQTL_m'] = meql_main

    else:
        globalenv['pval'] = args.p
        meql_main = robjects.r(r_str_F_option)
        globalenv['Mat_eQTL_m'] = meql_main

    # prepare the ggplot file
    globalenv['ggplot_file'] = output_filename_generator(output,
                                                         args.out_dir,
                                                         "_qqplot.tiff")

    robjects.r('''
        tiff(filename = ggplot_file)
    ''')

    robjects.r('''
        plot(Mat_eQTL_m)
    ''')

    robjects.r('''dev.off''')

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
