# GSEM Double Subtraction
A package to perform double genomic subtraction through Genomic SEM.

## Setup
This README file provides instructions on how to set up the project environment by installing the required dependencies for both R and Python. Having at least 8 GB of RAM memory is recommended, although this depends on the size of the summary statistics input files.

First of all, download the Github package and navigate in the package directory:

```console
git clone https://github.com/gpmerola/GSEM-DoubleSub.git
cd GSEM-DoubleSub
```

Then, install the R and Python dependecies.

### R Dependencies

To install the required R packages, run the following command:

```console
Rscript requirements.R
```

### Python Dependencies
To install the required Python packages, run the following command:

```console
pip install -r requirements.txt
```

### Usage

  1) Prepare Input Files: Ensure your input files are in the correct format and located in the appropriate directory.

  2) Choose Input Settings: The first argument contains the name of the directory you want your output files to be stored in. Add numbers from "1" to "6", or "r", to the input as shown below to run specific parts of the code:

     ```
      1: Preprocessing and Preparation - Performs munging, LD score regression, and prepares SNP files.
      
      2: Model Fitting - Fits the specified structural equation model using the input data.
      
      3: GSEM - Performs the most computationally intensive step, running a synthethic GWAS on the latent variables specified in the model. Including "r" runs the GWAS in a test mode, using only chromosome 2.
      
      4: Plots and Post-Munging - Generates Manhattan and QQ plots, performs post-munging, and computes LD score regression for the new phenotype.
      
      5: Genetic Correlation - Computes genetic correlation between the new phenotype and the input traits.
      
      6: Matrix Generation - Generates a matrix of genetic correlations and performs significance testing between the new phenotype and input traits.```

The parts must be run in order, with "3" being the most computationally intensive step.


  3) Run the Script: Execute the main script to perform the subtraction by navigating to the scripts directory and running the following command ("r" is optional):

```console
cd scripts
Rscript main.R  "working_directory" 1 2 3 4 5 6 "r"
```

  4) Plotting: If needed, the plotting script for genetic correlation can be run:

```console
python plot.py
```

### Settings

These variables are located at the top of the main.R file and can be edited to modify them:

```r
      1: files_input: List of file paths to cleaned GWAS summary statistics files. The third one represents the phenotype from which the subtraction is conducted.

      2: ref_file: File path to the reference panel.

      3: hm3: File path to the HapMap 3 SNP list.

      4: paths_corr: Directory path for the files for the correlation (see the "Correlation_input.csv" file, and the section below).

      5: ld: Directory path to the linkage disequilibrium (LD) reference data.

      6: wld: Directory path the weighted linkage disequilibrium (LD) reference data, if relevant. Otherwise set equal to "ld".

      7: traitnames: Names of traits for analysis. Has to follow the same order as "files_input".

      8: latentnames: Names of latent variables corresponding to the traits.  Has to follow the same order as "files_input".

      9: output_name: Name of the output for the synthetic phenotype file.

      10: infofilter: Information score filter threshold.

      11: maffilter: Minor allele frequency filter threshold.

      12: sample.prev: Vector of sample prevalence for each trait.  Has to follow the same order as "files_input".

      13: population.prev: Vector of population prevalence for each trait.  Has to follow the same order as "files_input".

      14: se.logit_vector: Logical vector indicating if standard error of logit transformation should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference).  Has to follow the same order as "files_input".

      15: OLS_vector: Logical vector indicating if Ordinary Least Squares (OLS) regression should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference).  Has to follow the same order as "files_input".

      16: linprob_vector: Logical vector indicating if linear probability model should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference).  Has to follow the same order as "files_input".

      17: ncores: Number of CPU cores to use for computation.
```

