#!/usr/bin/env python3
import argparse
import numpy as np
import xarray as xr

# Atomic weights in g/mol for common elements.
ATOMIC_WEIGHTS = {
    "Al": 26.9815,
    "Ba": 137.327,
    "Ca": 40.078,
    "Fe": 55.845,
    "K": 39.0983,
    "Mg": 24.305,
    "Na": 22.9898,
    "Si": 28.0855,
    "Li": 6.9400,
    "Mn": 54.9380,
    "S": 32.0600,
    "Sr": 87.6200
}

def convert_species_ppm_to_mM(species_name, ppm_values):
    """
    Convert concentration from PPM to milli molar (mM) for a given species.
    
    Parameters:
      species_name (str): The species name (e.g., "Ba 455.403 Ax"). The chemical
                          element is assumed to be the first token (here, "Ba").
      ppm_values (array-like): Concentration values in PPM.
      
    Returns:
      numpy array: The converted concentration values in mM.
    
    The conversion is performed using:
         mM = PPM / (atomic weight in g/mol)
    """
    element = species_name.split()[0]
    if element in ATOMIC_WEIGHTS:
        return ppm_values / ATOMIC_WEIGHTS[element]
    else:
        # If the element is not in our dictionary, return an array of NaNs.
        return np.full_like(ppm_values, np.nan)

def add_concentration_mM(ds):
    """
    For each species (from the 'species' coordinate) in the dataset,
    convert the corresponding column in the "Concentration" variable (in PPM)
    to milli molar (mM) and add the result as a new variable "Concentration_mM".
    
    Assumes the dataset has:
      - A variable "Concentration" with dimensions ("sample_name", "species")
      - A coordinate "species" containing species names.
      
    Returns:
      The dataset with the new variable "Concentration_mM".
    """
    # Extract the original concentration array (in PPM).
    conc_ppm = ds["Concentration"].values  # shape: (n_samples, n_species)
    n_samples, n_species = conc_ppm.shape
    conc_mM = np.empty_like(conc_ppm, dtype=float)
    
    # Loop over each species (i.e. each column) and perform conversion.
    for i, species in enumerate(ds["species"].data):
        conc_mM[:, i] = convert_species_ppm_to_mM(species, conc_ppm[:, i])
    
    # Add the new variable to the dataset.
    ds["Concentration_mM"] = (("sample_name", "species"), conc_mM)
    return ds

def convert_units(ds):
    """
    Main workflow:
      - Iterates over each species in the dataset and converts the "Concentration"
        variable (in PPM) to milli molar (mM) for that species.
      - Adds the converted data as a new variable "Concentration_mM" in the dataset.
    
    Returns:
      The updated xarray Dataset.
    """
    # Convert the concentration from PPM to mM and add as a new variable.
    ds = add_concentration_mM(ds)
    
    # Return the updated dataset.
    return ds


def compare_to_standards(ds):
    """
    Analyzes an xarray dataset to determine good and possible wavelengths based on certified standards.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to be analyzed.

    Returns
    -------
    tuple of dict
        A tuple containing two dictionaries:
            - good_wavelengths: mapping each standard to a list of species with percentage differences < 5%.
            - possible_wavelengths: mapping each standard to a list of species with percentage differences between 5% and 10%.
    """
    from check_stds import check_standards
    from utils import select_vars_by_element as sel_element
    from utils import select_by_fuzzy_sample as search_samples

    chk_stds = check_standards()
    good_wavelengths = {}
    possible_wavelengths = {}

    for standard in chk_stds:
        print(f'*** {standard} ***')
        good_wavelengths[standard] = []
        possible_wavelengths[standard] = []

        for element in ATOMIC_WEIGHTS.keys():
            print(f'*** {element} ***')
            try:
                certified = chk_stds[standard][0][element]
                uncertainty = chk_stds[standard][1][element]
                print(f'{standard} certified value: {certified:.3g} ± {uncertainty:.3g} mM')
            except KeyError:
                print(f'{standard} does not have a certified value for {element}')
                continue

            # Select the dataset for the current element and standard.
            na_ds = sel_element(ds, element)
            na_ds = search_samples(na_ds, standard)

            # Iterate over each species in the dataset.
            for specie in na_ds.coords['species'].values:
                concentration = na_ds.sel(species=specie)['Concentration_mM']
                mean = concentration.mean()
                two_sd = concentration.std() * 2
                percentage_diff = (mean - certified) / certified * 100
                print(f'{specie}: {mean:.3g} ± {two_sd:.3g} mM, {percentage_diff:.2f}%')

                if abs(percentage_diff) < 5:
                    good_wavelengths[standard].append(specie)
                elif 5 < abs(percentage_diff) < 10:
                    possible_wavelengths[standard].append(specie)

    return good_wavelengths, possible_wavelengths
