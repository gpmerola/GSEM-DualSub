# GSEM Dual Subtraction
A package to perform dual genomic subtraction on sumstats through [Genomic SEM](https://github.com/GenomicSEM/GenomicSEM).

## Overview of scripts
  - main.R: main script, with a customizable part at the beginning to set up variables and inputs.

  - plot.py: optional part of the code, to be executed after main.R, to plot the correlation graph.

  - Correlation_input.csv: input file (with available template) for the genetic correlation matrix.

## Prerequisites
  - R 3.6 or higher
  
  - Python 3.6 or higher

  - Git

## Setup
This README file provides instructions on how to set up the project environment by installing the required dependencies for both R and Python. Having at least 8 GB of RAM memory is recommended, although this depends on the size of the summary statistics input files.

First of all, download the Github package and navigate in the package directory:

```console
git clone https://github.com/gpmerola/GSEM-DualSub.git
cd GSEM-DualSub
```

Then, install the R and Python dependecies. To install the required R packages, run the following command:

```console
Rscript requirements.R
```

To install the required Python packages, run the following command:

```console
pip install -r requirements.txt
```

## Usage
  1) Prepare Input Files: Ensure your input files are in the correct format and located in the appropriate directory.

  2) Choose Input Settings: The first argument contains the name of the directory you want your output files to be stored in. Add numbers from "1" to "6", or "r", to the input as shown below to run specific parts of the code:

  - Preprocessing and Preparation - Performs munging, LD score regression, and prepares SNP files.
      
  - Model Fitting - Fits the specified structural equation model using the input data.
      
  - GSEM - Performs the most computationally intensive step, running a synthethic GWAS on the latent variables specified in the model. Including "r" runs the GWAS in a test mode, using only chromosome 2.
      
  - Plots and Post-Munging - Generates Manhattan and QQ plots, performs post-munging, and computes LD score regression for the new phenotype.
      
  - Genetic Correlation - Computes genetic correlation between the new phenotype and the input traits.
      
  - Matrix Generation - Generates a matrix of genetic correlations and performs significance testing between the new phenotype and input traits.

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

  - files_input: List of file paths to cleaned GWAS summary statistics files, 3 elements in the vector. The third one represents the phenotype from which the subtraction is conducted.

  - ref_file: File path to the reference panel.

  - hm3: File path to the HapMap 3 SNP list.

  - paths_corr: Directory path for the files for the correlation (see the "Correlation_input.csv" file, and the section below).

  - ld: Directory path to the linkage disequilibrium (LD) reference data.

  - wld: Directory path the weighted linkage disequilibrium (LD) reference data, if relevant. Otherwise set equal to "ld".

  - traitnames: Names of traits for analysis, 3 elements in the vector. Has to follow the same order as "files_input". Make sure that the third file in "traitnames" correspond to the first element in the "trait" column in "Correlation_input.csv"). 

  - latentnames: Names of latent variables corresponding to the traits, 3 elements in the vector. Has to follow the same order as "files_input".

  - output_name: Name of the output for the synthetic phenotype file.

  - infofilter: Information score filter threshold.

  - maffilter: Minor allele frequency filter threshold.

  - sample.prev: Vector of sample prevalence for each trait, 3 elements in the vector. Has to follow the same order as "files_input".

  - population.prev: Vector of population prevalence for each trait, 3 elements in the vector. Has to follow the same order as "files_input".

  - se.logit_vector: Logical vector indicating if standard error of logit transformation should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference), 3 elements in the vector. Has to follow the same order as "files_input".

  - OLS_vector: Logical vector indicating if Ordinary Least Squares (OLS) regression should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference), 3 elements in the vector. Has to follow the same order as "files_input".

  - linprob_vector: Logical vector indicating if linear probability model should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference), 3 elements in the vector. Has to follow the same order as "files_input".

  - ncores: Number of CPU cores to use for computation.

Additionally, the Correlation_input.csv file contains important information necessary for the genetic correlation analysis. Each column in the file is used to define specific parameters for the traits being analyzed. An example of the Correlation_input.csv file might look like this:

```csv
trait,code,sampleprev,popprev,Cluster
SCZ,SCHI06,0.425,0.01,A
MDD,DEPR14,0.346,0.10,B
BD,BIPO03,0.1013,0.02,C
```

  1) trait: The name of the trait.

  2) code: A unique identifier or code for each trait used, for file path construction (see "paths_corr" in the section above).

  3) sampleprev: The sample prevalence.

  4) popprev: The population prevalence.

  5) Cluster: This column lists the cluster information for each trait, used to group traits into clusters based on their genetic correlation patterns. The cluster information only affects the graphical representation in the analysis output.

Ensure that the data in Correlation_input.csv matches the order wanted in the final plot (with the exception of the first trait, which must be the one from which the subtraction is conducted). The file should be placed in the same directory as the .R file and the .py file.

## Output
The following output files are generated by the script during the analysis workflow. Each file serves a specific purpose in the overall analysis process, from preprocessing and modeling to visualization and correlation analysis. These files are saved in the new working directory created and set at the beginning of the script, provided through the console command.

  - dir.txt: Contains the working directory name. It helps share the information between the .R and the .py script.

  - LDSC_main.rds: Contains the LDSC output data.

  - SNP_files.rds: Contains the sumstats for the SNPs, harmonized.

  - model_LD.csv: Contains the factor loadings, SEs and p-values from the GSEM model.

  - SEMplot.png: The GSEM plot saved as a PNG image.

  - outputMANN.png: The Manhattan plot saved as a PNG image.

  - outputQQ.png: The QQ plot saved as a PNG image.

  - [latentname3]_gwas.gz: The GSEM results file for the synthetic phenotype.

  - trait_correlation_data.csv: Contains the trait correlation data. It helps share the information between the .R and the .py script.

  - Correlation_output.rds: Contains the correlation output data. It helps share the information between the .R and the .py script.

  - trait_correlation_plot_with_clusters_background.png: PNG image of the correlation plot with cluster backgrounds.

## Troubleshooting
  - Memory Issues: Ensure that you have sufficient RAM available. Try running the scripts on a machine with more memory if you encounter memory errors.
    
  - Dependency Issues: Verify that all dependencies are correctly installed by checking their versions and reinstalling if necessary.
    
  - File Paths: Double-check all file paths provided in the settings to ensure they are correct.
    
  - Permissions: Ensure you have read and write permissions for all directories and files used by the scripts.

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Citations
  - Genomic SEM:
    - Grotzinger, A.D., Mallard, T.T., Aivadyan, C., et al. (2019). Genomic SEM provides insights into the multivariate genetic architecture of complex traits. Nature Human Behaviour, 3(5), 513-525. DOI: 10.1038/s41562-019-0566-x.

  - LDSC (LD Score Regression):
    - Bulik-Sullivan, B.K., Loh, P.R., Finucane, H.K., et al. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nature Genetics, 47, 291-295. DOI: 10.1038/ng.3211.

  - HapMap 3:
    - International HapMap Consortium. (2010). Integrating common and rare genetic variation in diverse human populations. Nature, 467, 52-58. DOI: 10.1038/nature09298.

  - Reference Panel:
    - The 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. Nature, 526, 68-74. DOI: 10.1038/nature15393.
