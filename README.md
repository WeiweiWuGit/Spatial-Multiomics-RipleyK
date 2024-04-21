
# Ripley's K Based Statistic for RNA-Protein Colocalization Analysis

## Overview

This repository is dedicated to calculating both mark-weighted and unweighted Ripley's K functions and Ripley's K based statistic to quantify spatial patterns of RNA gene and protein expressions for colocalization analysis, using public [10x Visium melanoma data](https://www.10xgenomics.com/datasets/human-melanoma-if-stained-ffpe-2-standard) as demo.

## Repository Structure

```
|- code/
|   |-- data_preparation.py           # Prepares and preprocesses the input data
|   |-- Ripleys_K_function.py         # Core function calculate_RipleyK for weighted/unweighted Ripley's K value
|   |-- Ripleys_K_based_statistics.py # Function RipleyK_statistic for Ripley's K based statistic
|   |-- use_functions.py              # Demonstrates how to use the functions
|
|- data/
|   |-- raw/                          # Raw data from 10x Genomics
|   |-- processed/                    # Processed data files
|       |-- IF_matrix.csv             # Immunofluorescence data matrix
|       |-- protein_gene_matrix.csv   # Matrix combining protein and gene data
|
|- output/                            # Output directory for analysis results
|
README.md                             # Project documentation
```

## Setup

### Cloning the Repository

Clone this repository to your local machine using:
```bash
git clone https://github.com/WeiweiWuGit/Spatial-Multiomics-RipleyK.git
```

### Dependencies

Ensure Python 3.6 or newer is installed and install the required packages:
```bash
pip install numpy pandas scipy pyprojroot
```

### Prepare the Environment

Ensure that the script and data paths are set correctly, according to the project structure outlined above.

## Demo Data

This project uses publicly available 10x Visium data for melanoma [Human Melanoma, IF Stained (FFPE)](https://www.10xgenomics.com/datasets/human-melanoma-if-stained-ffpe-2-standard), which includes spatial gene expression profiles. The raw demo data are downloaded from the 10x Genomics website and provided in the `data/raw/` directory includes processed files used directly for analysis.

## Usage

### Data Preparation

`data_preparation.py` formats and preprocesses raw demo data from 10x Visium suitable for Ripley's K calculations. This script processes raw data from `data/raw/` and outputs to `data/processed/`, making it suitable for Ripley's K calculations.

### Functions: calculate_RipleyK( ) and RipleyK_statistic( )

- **calculate_RipleyK( )**: Computes both weighted and unweighted Ripley's K values for multiple RNAs and one protein.  

  **Parameters**:  
  - `x`, `y` (array): Arrays of x and y coordinates that represent the spatial locations of each observation point.
  - `protein_weights` (array): An array where each element corresponds to the weight (importance) of the protein measurement at the respective spatial location.
  - `gene_matrix` (array): A matrix containing gene expression data, where each row corresponds to a spatial location and each column to a different gene.
  - `r_min`, `r_max`, `r_step` (float, optional): Specify the minimum, maximum, and step size of radii (in units corresponding to those of `x` and `y`) to be used in the calculation of Ripley's K. Defaults to the entire range if not specified.
  - `area` (float, optional): The total area over which the points are distributed. This parameter influences the normalization of the Ripley's K statistic.
  - `weighted` (bool): If `True`, calculations consider the weights assigned to proteins; if `False`, the analysis is performed unweighted.
  - `n_genes` (int, optional): The number of genes to include in the analysis. Defaults to all columns in `gene_matrix` if not specified.
  - `return_df` (bool): Determines whether the function returns a pandas DataFrame (if `True`) or a numpy array (if `False`).
  - `gene_names` (list of str, optional): Custom labels for genes in the output; if not provided, genes are named based on their column indices in `gene_matrix`.
  
   
   

- **RipleyK_statistic( )**: Calculates Ripley's K-based statistic to quantify colocalization of multiple RNAs and one protein using data output from `calculate_RipleyK( )`.  

  **Parameters**:  
  - `weighted_RipleyK` (DataFrame): A pandas DataFrame containing the Ripley's K values computed with weighting, indexed by radius.
  - `unweighted_RipleyK` (DataFrame): A pandas DataFrame containing the Ripley's K values computed without weighting, also indexed by radius.
  - `radius_low` (int): The minimum radius value to include in the analysis, filtering out any results below this threshold.
  - `radius_up` (int, optional): The maximum radius value to include in the analysis, filtering out any results above this threshold. If not specified, all radii greater than `radius_low` are included.

### Demonstration

Refer to `use_functions.py` for examples on how to integrate and use these functions with the example data provided. This script will guide you through executing the functions and viewing the results. The preprocessed data used here is located at `data/processed/`.
