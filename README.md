# Analysis of LC-MS data from microfluidic effluent
LC-MS data and analysis workflow for idenitification of vitamins secreted by *Sphingomonas* Leaf257 when cultivated in a microfluidic device.

## General
The main script `workflow.py` determines absolute concentrations of the vitamins nicotinamide, dethiobiotin, pantothenate, biotin, and thiamine based on calibration curves for each compound. All required LC-MS/MS data sets are provided within this repository as mzML files (./Data/mzML_files_1to10, and ./Data/mzML_files_1to1). 

Samples were analyzed using a Thermo Ultimate 3000 UPLC instrument coupled to a QExactive plus mass spectrometer in parallel reaction monitoring (PRM) mode. Thermo raw files were
converted to mZML files using [MSConvert](https://proteowizard.sourceforge.io/download.html).

The workflow performs the following analysis steps: 
- Extraction of peaks based on provided retention time m/z windows (defined in ./Data/peaks_table.xlsx).
- Baseline subtraction of peaks based on extracted ion chromatogram using [pybaselines](https://pypi.org/project/pybaselines/).
- Fitting of s^2 weighted calibration curves for each vitamin (found in ./Data/Calibration_results).
- Calculation of absolute metabolite concentrations from calibration curve fitting parameters. Compounds whose estimated concentrations fell below the lowest measurable concentration of the corresponding standard were considered to be undetected.
- Generation of bar plots for vitamin concentrations in each condition. 

Outputs:
- A table with measured concentrations of compounds is generated and saved in './Data/results/peaks_table.xlsx'.
- A detailed results table containing the chromatograms for each measurement is saved in './Data/results/results_samples.table'. This table can be interactively explored using [emzed-gui](https://emzed.ethz.ch/get_started.html).

## Getting started
1. Install [emzed3](https://emzed.ethz.ch/get_started.html)
2. Open emzed-spyder app and install the following packages:
```python
pip install seaborn
pip install pybaselines
pip install uncertainties
pip instal PyYAML
```
3. Close the current IPython console and open a new IPython console
4. Clone this project and run the command 
```python
runfile('path/to/microfluidic-effluent-analysis/src/workfllow.py', wdir='path/to/microfluidic-effluent-analysis/src')
```
