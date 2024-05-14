# GSEM Double Subtraction
A package to perform double genomic subtraction through Genomic SEM.

module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6

## Setup
This README file provides instructions on how to set up the project environment by installing the required dependencies for both R and Python.

First of all, download the Github package and navigate there:

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

Follow these steps to perform double genomic subtraction:

  1) Prepare Input Files: Ensure your input files are in the correct format and located in the appropriate directory.

  2) Run the Script: Execute the main script to perform the subtraction:


```console
cd scripts
Rscript main.R "1" "2" "3" "4" "5" "6" "r" 
```


