# PyReQTL: Collections of Python Modules Equivalent to R ReQTL Toolkit

## Introduction

I present to you a collection of python modules which were developed for individuals or research 
groups who prefer to perform ReQTL analysis using python instead of R. Hence, this Python modules are replica of 
the R toolkit **"ReQTL"** but implemented in Python (please see the section References for direct link of the
ReQTL R toolkit GitHub page). The original R toolkit requires a number of R packages such as 
MatrixEQTL and among others Python lacks equivalent libraries. Furthermore, despite that many of the code-base in this 
particular project were smoothly translated into python, there were consistence needs to interface with R and its 
environment including importing R packages to utilizing their functions for instance from within Python and thus 
the Python library "rpy2" and R are requirements to these modules.


Briefly, ReQTL does identify the correlation between expressed SNVs and their gene expression using RNA-sequencing data 
and each module transforms the sequencing files into ReQTL input files which will eventually be used 
as inputs for the R package MatrixEQTL to find the significant variation-expression relationships. 


## Asking for Help

If you have an issue (error or bug) when executing any of the individual module or via the prepared shell script, 
or have any other enquires please feel free to send me an [email](mailto:nawafalomran@hotmail.com)


## Release
This is the first release version of PyReQTL. 


## Getting Started

Please follow the instructions in the sections that follow. You need to have a terminal to be able to execute the 
modules via the command line. In case you are getting an error or bug please refer to the section **Asking for Help**.


## Operating System Requirement

PyReQTL modules were developed on *macOS Catalina v.10.15.2* and should work also on *Linux OS*.


## Python Version 
Both Python version *3.7* and *3.8* were successfully tested on PyReQTL. You need Python version >=3.5.

**Note:**

- I would recommend checking your default python version which it may be python *version 2.7*  so running ```python```
may start version *2.6* or *2.7* and therefore you need to run ```python3``` to start python 3 version 
or you may set your python3 as your default python version.


- Note that PyReQTL were not tested on Windows and it may or may not work.

## Prerequisites

* The following are the core Python libraries (dependencies) required for PyReQTL modules:

  ```
  pandas
  numpy
  scipy
  rpy2
  ```
  
* Optional Python Library:
    ```
  mypy
  ```

* You will also need to install the latest version of *R (4.0.2)* from the the official language site 
at: https://cloud.r-project.org
  
* #### **Important notes:** 
  1. Python will install the required R packages for you but you may be asked to choose a secure mirror for installing 
  those R packages.
  
  2. You may also be asked when executing the modules that require R packages to upgrade some R dependencies, you may 
  ignore it by typing```no```.
  
  3. Some bioconductor packages will also be required for the bioinformatics related analysis.


## Install 
1.  First, you need to download the compressed PyReQLT modules or clone it from GitHub

2.  If you choose to download the modules from the zip file, then upzip it and change to the directory *(cd)* 
where the modules live (at the *source* directory):
      ```bash
    unzip PyReQTL-master.zip
    cd PyReQTL-master/source
      ```    
3. Its highly recommended to create conda or python virtual environment to install the dependencies within such 
environment (OPTIONAL):

     - Creating conda virtual environment:
     
        ```bash
        conda create -n your_env_name
        conda activate your_env_name
       ```
     - Creating python virtual environment:
     
       ```bash
        pip install virtualenv
        virtualenv your_env_name
        source your_env_name/bin/activate
       ```
       
     - **Note:** please don't forget to deactivate your environment of choice once you're done running the modules:
       ```bash
            # for conda env
            conda deactivate
       
            # for python env
            deactivate
       ```

4. You may now try to install the dependencies by running the command:
    ```bash
    sudo pip install -r requirements.txt
   ```   
   
5. Then, issue the  ```ls ``` command to list the files that end with .py (i.e. the PyReQTL modules) (OPTIONAL)
   ```bash
    ls -1 *.py    
   ```  
   You should see the following list of python scripts:
   ```                                                                                                                                                  ✔  PyReQTL 🐍  09:08:09 
    annotate_cis_trans.py
    build_VAF_matrix.py
    build_gen_exp_matrix.py
    commons.py
    harmonize_matrices.py
    run_matrix_ReQTL.py
    ```
   
## Running the Modules

You are given two options to run the modules:
 1. Running the prepared shell script *"run_all.sh"* found at the parent's directory of *source*. This script is 
 prepared to automatically run all modules sequentially with just a single command.
 
 2. Running each individual python module separately.


### - Option One (run all the modules sequentially)
   
Run the whole analysis at once:

1. type ```cd ..``` to go back to the main directory in which the shell script lives

    - You can modify the script based on your needs
 
2. run the following command:

```bash
sh run_all.sh
```

### - Option Two (run each module separately)


#### 1. build_gen_exp_matrix module

##### **Note**: 

1. please refer to the link in the *References* section for the description of each module
2. make sure you are in the ```source``` directory in order to run each individual python module

- Run the following command to go the source directory:
 ```bash
cd source/
```
    
##### Inputs + Options

  + -i: a directory contain all the raw gene expression files from Stringtie

  + -o: the prefix for the gene expression and gene location files

##### Output
  * the gene expression values to be fed into MatrixEQTL

  * the gene location for MatrixEQTL

#### Command Example 
```bash
python build_gen_exp_matrix.py -i ../data -o ReQTL_test
```

***

#### 2. build_VAF_matrix module

##### Inputs + Options

  + -r: a directory containing the .csv files from the output of readCounts
    module

  + -o: the prefix for the SNV matrix and SNV location files

##### Output

  + the SNV locations for MatrixEQTL

  + the SNV variant allele fraction matrix for MatrixEQTL

##### Command Example
```bash
python build_VAF_matrix.py -r ../data -o ReQTL_test
```

***


#### 3. harmonize_matrices module

##### Input

 + -v: the path to the variant matrix generated from the build_VAF_matrix.py

 + -g: the path to the gene expression matrix generated from the build_gen_exp_matrix.py

 + -c: the path to the covariate matrix located in the data folder

##### Output
 *  three output files consists of the three input files that have mutual samples
 

##### Command Example
```bash
python harmonize_matrices.py \
    -v output/ReQTL_test_VAF_matrix.txt \
    -g output/ReQTL_test_gene-exp_matrix.txt \
    -c ../data/covariates_matrix.txt
```
##### Note: 
1. the covariates matrix file is supplied in the "data" folder
2. this is an optional module but it is recommended to avoid possible error or bug

***

#### 4. run_matrix_ReQTL module

This module requires the R package MatrixEQTL. Please refer to the section References for more 
information about this package 

##### Input

  + -s: the SNV or variant matrix file created from harmonize_matrices

  + -sl: the SNV location matrix file created from build_VAF_matrix

  + -ge: the gene expression file matrix created from build_gene-exp_matrix

  + -gl: the gene locations file created from build_gene-exp_matrix

  + -c: the covariates matrix file created from "harmonize_matrices.
    [OPTIONAL]! you can access the file under data directory

  + -o: the prefix for the path to the output files

  + -ct: logical (T or F) specifying whether to split the output into cis
    or trans

  + -pcis: p-value thresholds for the cis output files

  + -ptran: p-value thresholds for the trans output files

  + -p: p-value thresholds for the unified output file


##### Output

  + one output could be cis ReQTLs and trans ReQTLs or one with all of the
    unified ReQTLs. This depends on what you choose for the value of parameter
    "-ct"

  + QQ plot of p-values

##### Commands Example

* splitting *cis* and *trans*
```bash
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
```

* unified *cis* and *trans*
```bash
python run_matrix_ReQTL.py \
       -s output/ReQTL_test_VAF_matrix_harmonized.txt \
       -sl output/ReQTL_test_VAF_loc_matrix.txt \
       -ge output/ReQTL_test_gene-exp_matrix_harmonized.txt \
       -gl output/ReQTL_test_gene-exp-loc_matrix.txt \
       -c output/covariates_matrix_harmonized.txt \
       -ct F \
       -o "ReQTL_test" \
       -p 0.001
```
***

#### 5. annotate_cis_trans module

#### Input

  + -r: the path to the ReQTL analysis result file

  + -ga: the path to the file gene location annotations

  + -o: the prefix for the output annotated result



#### Output

  * the ReQTLs annotated as cis or trans

##### Commands Example
```bash
python annotate_cis_trans.py \
    -r output/ReQTL_test_all_ReQTLs.txt \
    -ga ../data/gene_locations_hg38.txt \
    -o ReQTL_test
```
***


## Features
   - Faster than the R ReQTL version 
   - For efficiency and performance, it uses both Pandas DataFrames and numpy arrays
   - rpy2 is more efficient and better integrated than using subprocess to interface with R 
   

## Authors & Acknowledgment for the ReQTL R Toolkit


- Liam Spurr, **Nawaf Alomran**, Pavlos Bousounis, Dacian Reece-Stremtan, Prashant N M, Hongyu Liu, Piotr Słowiński, Muzi Li, Qianqian Zhang, Justin Sein, Gabriel Asher, Keith A. Crandall, Krasimira Tsaneva-Atanasova, and Anelia Horvath 

- We would like to thank the Matrix EQTL team (Shabalin, et al. 2012) for their sample code and R package 
upon which *run\_matrix_ReQTL.R* is based.


## Paper
- https://academic.oup.com/bioinformatics/article/36/5/1351/5582649

## References

1. The ReQTL R toolkit github repository: https://github.com/HorvathLab/ReQTL

2. MatrixEQTL R package: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

2. readCounts module: https://github.com/HorvathLab/NGS/tree/master/readCounts
