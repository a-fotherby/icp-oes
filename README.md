# ICP-OES data processing

A Python-based data processing and analysis toolkit for handling and visualizing scientific data, specifically focused on ICPOES analysis.

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Modules](#modules)
- [Dependencies](#dependencies)
- [License](#license)

## Introduction
This project provides a set of tools for processing, analyzing, and visualizing data related to ICPOES (Inductively Coupled Plasma Optical Emission Spectroscopy). It includes functions for data import, unit conversion, standard checking, error handling, and graphical representation of results.

## Features
- **Data Import:** Import ICP-OES data for processing.
- **Standard Checking:** Validate data against predefined standards.
- **Error Handling:** Detect and manage errors in data processing.
- **Plotting:** Generate visual representations of data.
- **Unit Conversion:** Convert between different measurement units.

## Installation
To install the required dependencies, run:

```bash
pip install -r requirements.txt
```

Ensure you have Python 3.x installed.

## Usage
1. Run `processing.ipynb` to import data from raw and save as `netcdf`.
2. Run `errors.ipynb` do to error analysis based on check standards. Overwrites previous `.nc`. 

Example:

```python
from import_icpoes import load_data
from process import process_data
from plotting import plot_results

data = load_data("datafile.csv")
processed_data = process_data(data)
plot_results(processed_data)
```

## Modules
### `check_stds.py`
Handles validation of data against predefined standards.

### `errors.py`
Contains error handling mechanisms.

### `import_icpoes.py`
Responsible for importing ICPOES data from various file formats.

### `plotting.py`
Provides functions for visualizing processed data.

### `process.py`
Processes raw data, applies necessary transformations, and prepares it for analysis.

### `unit_conv.py`
Performs unit conversions necessary for standardizing data.

### `utils.py`
Contains helper functions used across different modules.

## Dependencies
- Python 3.x
- pandas
- matplotlib
- numpy

## License
This project is licensed under the MIT License.


