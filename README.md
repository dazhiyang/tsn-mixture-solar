# Solar Radiation Distribution and Classification

This repository contains the code and supporting files for a study on solar radiation distribution modeling using mixture models (specifically Truncated Skew-Normal distributions) and radiation climate classification.

## Project Structure

### `Code/`
This directory contains the R scripts used for analysis and visualization:
- **`0.Functions.R`**: Core utility functions, including probability density functions for Skew-Normal (SN) and Truncated Skew-Normal (TSN) distributions, and Wasserstein distance calculations.
- **`1.MainLoop.R`**: The primary analysis script that iterates through station data to perform model fitting.
- **`3.PerfWorld.R`**: Scripts to evaluate and compare model performance across stations.
- **`4.Classification.R`**: Implements the Radiation Climate Classification (RCC) using clustering techniques based on station distribution similarity.
- **`5.Sensitivity.R`**: Analysis of model sensitivity to sample size and parameter stability.
- **`Fig.126.R`**: Generates summary plots for the 126 stations analyzed.
- **`Fig.PAL.R`**: Creates palette visualizations showing the individual sub-components of the mixture models.
- **`Fig.SN.R`**: Generates figures specifically for the Skew-Normal distribution model.
- **`Fig.TSN.R`**: Generates figures specifically for the Truncated Skew-Normal distribution model.

### `Data/`
Input data and reference files:
- **`location.csv`**: Metadata (station IDs, names, coordinates) for the 126 measurement sites.
- **`Beck_KG_V1_present_0p083.tif`**: Raster dataset for the Koppen-Geiger climate classification system.
- **`climatology variables of 126 sites.csv`**: Summary of climatological variables for the studied locations.
- **`Sample/`**: Contains original R object data for one sample station (ADE). Due to the large total size of the dataset, only this representative sample is provided to demonstrate the code's functionality.
- **`FigDist/`**: Generated distribution plots for individual stations.
- **`FigPIT/`**: Generated Probability Integral Transform (PIT) plots for model validation.
- *Note: `Data/Results` and `Data/Sensitivity` are excluded from version control.*

### `tex/`
Manuscript and output files:
- *Note: The `tex/` folder, containing the research paper LaTeX source and generated publication figures, is excluded from version control.*

## Getting Started

1. Ensure the `dir0` path in the R scripts is set to your local project root.
2. The analysis relies on data stored in `Data/Sample`, but intermediate results should be stored locally in `Data/Results` (which is not tracked by Git).
3. Source `0.Functions.R` before running other analysis scripts.
