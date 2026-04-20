# MetaboPeak <img src="man/figures/logo.png" align="right" height="134" alt="" />

MetaboPeak is an R package for preprocessing and exploring LC-MS metabolomics data. It focuses on filtering, quality control, and visualization of peak data to support reproducible workflows.

## You can install the development version from GitHub:
install.packages("devtools")

devtools::install_github("rennadiv/MetaboPeak")

## What problem does it solve?
LC-MS peak tables often require multiple preprocessing steps:
- filtering peaks with high missing values
- removing noisy features (high CV)
- selecting relevant peaks
- visual inspection of peak intensities

These steps are often done manually or inconsistently.
MetaboPeak provides a structured and reproducible workflow for these tasks.

## Example workflow and description of functions
<img width="521" height="753" alt="workflow_function" src="https://github.com/user-attachments/assets/1846086d-587c-4918-906a-beb3e148ef15" />

## Input data format

The input data frame should contain:

- identificator of peaks
- intensity values for each sample (numeric columns)
- a column with m/z values (e.g. m.z)
- retention time (RT)

Example structure:

<img width="311" height="227" alt="dataframe" src="https://github.com/user-attachments/assets/9b86e926-3a19-4642-a9cc-88a127c42dc4" />

## Testing

The package includes unit tests using testthat to ensure correct functionality and robustness.

## Citation

If you use MetaboPeak in your work, please cite the associated publication (in preparation).
