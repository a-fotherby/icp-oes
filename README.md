# Concentration Analyser

A Python-based tool for processing spectrometer data, computing calibration statistics, estimating measurement errors, and aggregating results by chemical element.

---

## Overview

The Concentration Analyser is designed to ingest spectrometer data from CSV files, clean and format it into an xarray dataset, and execute an error analysis pipeline. The primary goal is to:
- Calculate calibration statistics (mean and standard deviation in parts-per-million, PPM) from a subset of calibration samples.
- Estimate measurement errors for each wavelength (or species) by comparing the calibration mean to certified values.
- Filter and aggregate the data by chemical element, retaining only those wavelengths that pass strict quality control.

---

## Features

- **Data Import and Cleaning:**  
  - Reads CSV files containing spectrometer measurements and metadata.
  - Cleans data by stripping extra whitespace, converting numeric columns, handling missing values, and merging element labels.
  - Converts the cleaned DataFrame into an xarray dataset with meaningful coordinates (`sample_name` and `species`).

- **Calibration Statistics:**  
  - Selects calibration samples based on a naming convention (e.g., samples starting with `check_std_`).
  - Computes the mean and standard deviation (SD) of concentration measurements for each wavelength.

- **Error Estimation:**  
  - Calculates the percentage difference between calibration means and provided certified values.
  - Determines two types of error estimates:
    - **Percentage-based error:** Scaled from the calibration mean.
    - **Variation-based error:** Two times the calibration SD.
  - Chooses the maximum of these two error estimates for each calibration standard and species.
  - Applies a quality control check by only considering check standards with an absolute percentage difference of less than 5%.
  - Broadcasts the maximum valid error for each species across all samples.

- **Aggregation by Element:**  
  - Filters wavelengths to include only those with valid error estimates across all samples.
  - Groups wavelengths by chemical element (extracted from species names) and averages their concentrations and errors.
  - Records the list of contributing wavelengths for each element as metadata.

- **Output:**  
  - Combines raw concentration measurements and associated errors into a final results dataset.
  - Exports the final dataset to a NetCDF file.

---

## Requirements

- Python 3.x
- [xarray](https://xarray.pydata.org/)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [matplotlib](https://matplotlib.org/)

Install the dependencies via pip:

```bash
pip install xarray numpy pandas matplotlib
```

---

## Installation

Clone this repository or download the script file to your local machine.

```bash
git clone https://github.com/yourusername/concentration-analyser.git
cd concentration-analyser
```

---

## Usage

The script is executed from the command line and expects two arguments:
1. **Input CSV file path** – The CSV file containing spectrometer data.
2. **Output NetCDF file name** – The name for the output file that will store the final results.

Run the script as follows:

```bash
./concentration_analyser.py input_file.csv output_file.nc
```

### Example

```bash
./concentration_analyser.py data/spectrometer_data.csv results/output.nc
```

---

## Code Structure

- **Data Importing Functions:**  
  - `_load_esws_csv()`, `_clean_data()`, `_combine_element_labels()`:  
    Functions to load and clean the CSV data, ensuring it is formatted correctly.

- **Dataset Construction:**  
  - `_to_xarray_dataset()`, `_update_label_names()`, `_remove_extraneous_variables()`, `_move_unit_to_concentration()`, `_sanitize_netcdf_names()`:  
    Functions to convert the cleaned DataFrame into an xarray dataset with proper coordinates and variable names.

- **Calibration and Error Analysis Pipeline:**  
  - `compute_calibration_stats()`:  
    Identifies calibration samples and computes the mean and SD of concentrations.
  
  - `compute_error_estimates()`:  
    Computes the percentage difference between calibration means and certified values, calculates two error estimates (percentage-based and variation-based), and chooses the maximum of these. It then filters out unreliable calibration data based on a 5% threshold and broadcasts the final error for each species.
  
  - `get_results()`:  
    Combines the concentration measurements with the computed error estimates to generate the final dataset.
  
  - `aggregate_by_element()`:  
    Aggregates the data by averaging wavelengths (species) corresponding to the same chemical element and records the contributing wavelengths.

- **Main Execution Block:**  
  - Parses command-line arguments, creates an analyser instance, runs the full pipeline (data import, calibration, error computation, aggregation), and writes the final results to a NetCDF file.

---

## Contributing

Contributions, bug reports, and feature requests are welcome! Please create an issue or submit a pull request on the GitHub repository.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---
