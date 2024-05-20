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

```
      1: files_input: List of file paths to cleaned GWAS summary statistics files, 3 elements in the vector. The third one represents the phenotype from which the subtraction is conducted.

      2: ref_file: File path to the reference panel.

      3: hm3: File path to the HapMap 3 SNP list.

      4: paths_corr: Directory path for the files for the correlation (see the "Correlation_input.csv" file, and the section below).

      5: ld: Directory path to the linkage disequilibrium (LD) reference data.

      6: wld: Directory path the weighted linkage disequilibrium (LD) reference data, if relevant. Otherwise set equal to "ld".

      7: traitnames: Names of traits for analysis, 3 elements in the vector. Has to follow the same order as "files_input". Make sure that the third file in "traitnames" correspond to the first element in the "trait" column in "Correlation_input.csv"). 

      8: latentnames: Names of latent variables corresponding to the traits, 3 elements in the vector. Has to follow the same order as "files_input".

      9: output_name: Name of the output for the synthetic phenotype file.

      10: infofilter: Information score filter threshold.

      11: maffilter: Minor allele frequency filter threshold.

      12: sample.prev: Vector of sample prevalence for each trait, 3 elements in the vector. Has to follow the same order as "files_input".

      13: population.prev: Vector of population prevalence for each trait, 3 elements in the vector. Has to follow the same order as "files_input".

      14: se.logit_vector: Logical vector indicating if standard error of logit transformation should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference), 3 elements in the vector. Has to follow the same order as "files_input".

      15: OLS_vector: Logical vector indicating if Ordinary Least Squares (OLS) regression should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference), 3 elements in the vector. Has to follow the same order as "files_input".

      16: linprob_vector: Logical vector indicating if linear probability model should be used (https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information for reference), 3 elements in the vector. Has to follow the same order as "files_input".

      17: ncores: Number of CPU cores to use for computation.
```

Additionally, the Correlation_input.csv file contains important information necessary for the genetic correlation analysis. Each column in the file is used to define specific parameters for the traits being analyzed.

```csv
trait: This column lists the names of the traits being analyzed. Each row corresponds to a different trait.

code: This column provides unique codes associated with each trait, used for file path construction.

sampleprev: This column specifies the sample prevalence for each trait. The values should be provided in numeric format and represent the proportion of samples exhibiting the trait within the study population.

popprev: This column specifies the population prevalence for each trait. The values should be provided in numeric format and represent the proportion of the general population exhibiting the trait.

Cluster: This column lists the cluster information for each trait, used to group traits into clusters based on their genetic correlation patterns.
```

An example of the Correlation_input.csv file might look like this:

```csv
Copia codice
trait,code,sampleprev,popprev,Cluster
SCZ,SCHI06,0.425,0.01,A
MDD,DEPR14,0.346,0.10,B
BD,BIPO03,0.1013,0.02,C
trait: The name of the trait (e.g., SCZ, MDD, BD).
code: A unique identifier or code for each trait (e.g., SCHI06, DEPR14, BIPO03).
sampleprev: The sample prevalence (e.g., 0.425, 0.346, 0.1013).
popprev: The population prevalence (e.g., 0.01, 0.10, 0.02).
Cluster: Cluster grouping for each trait (e.g., A, B, C).
```

Ensure that the data in Correlation_input.csv matches the order and values required for the analysis in the script. The file should be correctly formatted and saved in the appropriate directory for the script to access and utilize the data effectively.








