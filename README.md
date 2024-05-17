# GSEM Double Subtraction
A package to perform double genomic subtraction through Genomic SEM.

## Setup
This README file provides instructions on how to set up the project environment by installing the required dependencies for both R and Python. Having at least 10 GB of RAM memory, is recommended, although this depends on the size of the summary statistics input files.

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


