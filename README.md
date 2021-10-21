# MSeQUiP
Identify peptide-spectrum matches (PSMs) from LC-MS/MS data using the method described in https://arxiv.org/abs/2110.07818. The main function is , which uses  described in  to infer PSMs.

# Functions
The main functions are:

## PerformSearch
Use spectral libraries and our Bayes factor scoring function to infer PSMs.

## Deisotope
Use our novel MS/MS deisotoping and charge assignment algorithm to clean MS/MS spectra.

## SimulateDataPaper
Simulate realistic LC-MS/MS data using the method described in https://arxiv.org/abs/2110.07818.

## AnalyzeSimulatedData
Analyze the above simulated.

# Reproduce simulation results
To reproduce simulation results from https://arxiv.org/abs/2110.07818:
1) Download the data from https://drive.google.com/file/d/1Iz6WSTIQV-8nj5qPIKyrGv7NkVTW5YI1/view?usp=sharing
2) Unzip the folder and follow the instructions in Simulations.Rmd


