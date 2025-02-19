import xarray as xr

def append_stats_to_dataarray(dataarray: xr.DataArray, sample_subset: list, label_prefix: str) -> xr.DataArray:
    """
    Calculate the mean and 2* standard deviation for a subset of the 'sample_name'
    coordinate (given by sample_subset) along the 'species' dimension, and append these
    as new entries on the 'sample_name' coordinate. The new sample labels are formed by
    appending '_mean' and '_2sd' to the provided label_prefix.

    Parameters
    ----------
    dataarray : xr.DataArray
        The input DataArray. It should have at least the coordinates 'sample_name' and 'species'.
    sample_subset : list
        A list of labels from the 'sample_name' coordinate that will be used to compute the stats.
    label_prefix : str
        A string used to form the names of the new sample_name entries (by appending '_mean'
        and '_2sd').

    Returns
    -------
    xr.DataArray
        The updated DataArray with two additional entries along the 'sample_name' coordinate.
    """
    # Select the subset of the data array for the given sample names.
    subset = dataarray.sel(sample_name=sample_subset)
    
    # Compute the mean and standard deviation along the 'sample_name' dimension, per species.
    mean_vals = subset.mean(dim='sample_name')
    two_sd_vals = 2 * subset.std(dim='sample_name')
    
    # Create new DataArrays for the computed statistics.
    # We expand the dimensions so that each statistic has its own 'sample_name' coordinate.
    mean_da = mean_vals.expand_dims(dim={'sample_name': [f"{label_prefix}_mean"]})
    two_sd_da = two_sd_vals.expand_dims(dim={'sample_name': [f"{label_prefix}_2sd"]})
    
    # Concatenate the new DataArrays with the original one along the 'sample_name' dimension.
    updated_dataarray = xr.concat([dataarray, mean_da, two_sd_da], dim='sample_name')
    
    return updated_dataarray

import xarray as xr
import numpy as np

def append_percentage_difference(dataarray: xr.DataArray,
                                 certified_values: dict,
                                 base_sample: str) -> xr.DataArray:
    """
    Calculate the percentage difference between a computed mean (located at the sample
    label "{base_sample}_mean" in the 'sample_name' coordinate) and a dictionary of certified
    values for each species. The species coordinate labels are assumed to be formatted as
    "{ELEMENT_SYMBOL}_average", so that the element symbol can be extracted (by removing
    the '_average' suffix) to match the keys in the certified_values dictionary.

    The percentage difference is computed as:
        (computed_mean - certified_value) / certified_value * 100

    This function then appends the calculated percentage differences as a new sample (with
    label "{base_sample}_pct_diff") along the 'sample_name' coordinate.

    Parameters
    ----------
    dataarray : xr.DataArray
        The input DataArray containing a computed mean sample.
    certified_values : dict
        Dictionary of certified values with element symbols as keys.
    base_sample : str
        The base name used to construct the computed mean label ("{base_sample}_mean") and
        the percentage difference label ("{base_sample}_pct_diff").

    Returns
    -------
    xr.DataArray
        The updated DataArray with the new percentage difference sample appended.
    """
    # Generate the labels automatically.
    mean_sample_label = f"{base_sample}_mean"
    pct_diff_label = f"{base_sample}_pct_diff"
    
    # Extract the computed mean sample using the generated mean_sample_label.
    mean_sample = dataarray.sel(sample_name=mean_sample_label)
    
    # Get the species labels.
    species_labels = mean_sample.coords['species'].values

    # Build a list of certified values corresponding to the species.
    cert_values_list = []
    for species_label in species_labels:
        # Extract the element symbol by removing the '_average' suffix.
        element_symbol = species_label.replace('_average', '')
        if element_symbol not in certified_values:
            print(f"Certified value for element '{element_symbol}' not found in the provided dictionary.")
            cert_values_list.append(np.nan)
        else:
            cert_values_list.append(certified_values[element_symbol])
    
    # Create a DataArray of certified values aligned with the species coordinate.
    cert_da = xr.DataArray(cert_values_list, coords={'species': species_labels}, dims=['species'])
    
    # Compute the percentage difference:
    #   (computed mean - certified value) / certified value * 100
    pct_diff = (mean_sample - cert_da) / cert_da * 100

    # Expand dimensions so that the result has a new 'sample_name' coordinate.
    pct_diff = pct_diff.expand_dims(dim={'sample_name': [pct_diff_label]})
    
    # Append the new percentage difference sample along the 'sample_name' coordinate.
    updated_dataarray = xr.concat([dataarray, pct_diff], dim='sample_name')
    
    return updated_dataarray


def append_error_estimate(dataarray: xr.DataArray, base_sample: str) -> xr.DataArray:
    """
    Calculate the error estimate for a given base sample by comparing the absolute 
    concentration error derived from the percentage difference with the 2* standard 
    deviation value. In this updated version, the absolute error from the percentage 
    difference is calculated using the BASE_SAMPLE_mean value instead of a certified value.
    
    The error from the percentage difference is computed as:
        error_from_pct = (abs(pct_diff) / 100) * BASE_SAMPLE_mean

    For each species, the final error is taken as the maximum of error_from_pct and the 
    2* standard deviation value (BASE_SAMPLE_2sd). The new error sample is then appended 
    with the label "{BASE_SAMPLE}_error" along the 'sample_name' coordinate.

    Expected input samples in the dataarray:
      - A computed mean sample labeled "{BASE_SAMPLE}_mean"
      - A percentage difference sample labeled "{BASE_SAMPLE}_pct_diff"
      - A 2* standard deviation sample labeled "{BASE_SAMPLE}_2sd"

    Parameters
    ----------
    dataarray : xr.DataArray
        Input DataArray containing the required samples.
    base_sample : str
        The base prefix for the sample labels. For example, if base_sample is "sample1",
        the function will look for "sample1_mean", "sample1_pct_diff", and "sample1_2sd", 
        and will create "sample1_error".

    Returns
    -------
    xr.DataArray
        Updated DataArray with the new error sample appended along the 'sample_name' coordinate.
    """
    # Automatically generate the labels
    pct_diff_label = f"{base_sample}_pct_diff"
    twosd_label    = f"{base_sample}_2sd"
    mean_label     = f"{base_sample}_mean"
    error_label    = f"{base_sample}_error"

    # Extract the relevant samples from the dataarray
    pct_diff_sample  = dataarray.sel(sample_name=pct_diff_label)
    twosd_sample     = dataarray.sel(sample_name=twosd_label)
    base_mean_sample = dataarray.sel(sample_name=mean_label)

    # Convert the percentage difference to an absolute concentration error
    # using the base sample mean
    error_from_pct = (np.abs(pct_diff_sample) / 100) * base_mean_sample

    # For each species, select the larger error between error_from_pct and twosd_sample.
    error_final = xr.apply_ufunc(np.maximum, error_from_pct, twosd_sample)

    # Expand dimensions to add a new 'sample_name' coordinate for the error sample.
    error_sample = error_final.expand_dims(dim={'sample_name': [error_label]})

    # Append the new error sample along the 'sample_name' coordinate.
    updated_dataarray = xr.concat([dataarray, error_sample], dim='sample_name')
    
    return updated_dataarray

import re

def add_error_variable(sample_dict, data):
    """
    Creates a new dataset that includes a new variable 'error' computed from the provided dictionary.
    
    Parameters
    ----------
    sample_dict : dict
        A dictionary where keys are sample names (e.g., "sample1") and values are lists of strings.
        Each string is formatted like "Ba 455.403 Ax" (an element symbol, a number, and then either "R" or "Ax").
        If an element appears more than once for a given key, it is counted only once.
    
    data : xr.DataArray
        An xarray DataArray with two coordinates:
            - 'species': labels formatted as "Element_average" (e.g., "Ba_average")
            - 'sample_name': labels formatted as "{sample_key}_error" (e.g., "sample1_error")
        The data array represents a variable named "Concentration_mM".
    
    Returns
    -------
    xr.Dataset
        A dataset containing two variables:
            - "Concentration_mM": the original data array.
            - "error": a new variable that, for each species, contains:
                * the corresponding error value (from the sample error entry) if the element is mentioned in one dictionary key,
                * the maximum of the corresponding error values if the element is mentioned in more than one key,
                * NaN if the element is not mentioned in any dictionary entry.
        The "error" variable is indexed only by the 'species' coordinate.
    """
    
    # Regular expression to parse strings like "Ba 455.403 Ax" (allows for extra spaces)
    pattern = re.compile(r'^\s*([A-Z][a-z]?)\s*\d+(?:\.\d+)?\s*(R|Ax)\s*$')
    
    # Build a mapping: each sample key -> set of unique element symbols found in its list of strings.
    sample_to_elements = {}
    for sample_key, lst in sample_dict.items():
        elems = set()
        for s in lst:
            match = pattern.match(s)
            if match:
                elem = match.group(1)
                elems.add(elem)
        sample_to_elements[sample_key] = elems

    # Invert the mapping: each element -> list of sample keys in which it appears.
    element_to_samples = {}
    for sample_key, elems in sample_to_elements.items():
        for elem in elems:
            element_to_samples.setdefault(elem, []).append(sample_key)
    
    # Compute an error value for each species.
    # For each species label (formatted like "Ba_average"), extract the element symbol ("Ba").
    error_values = []
    for species_label in data.coords['species'].values:
        element = species_label.split('_')[0]
        if element in element_to_samples:
            sample_keys = element_to_samples[element]
            error_candidates = []
            for sample_key in sample_keys:
                sample_err_label = f"{sample_key}_error"
                # Only consider error values from existing sample_name coordinates.
                if sample_err_label in data.coords['sample_name'].values:
                    # Retrieve the error value from data at the intersection of species and the sample error label.
                    err_val = data.sel(species=species_label, sample_name=sample_err_label).values.item()
                    error_candidates.append(err_val)
            if error_candidates:
                # Use the maximum error if there is more than one candidate.
                computed_error = max(error_candidates) if len(error_candidates) > 1 else error_candidates[0]
            else:
                computed_error = np.nan
        else:
            computed_error = np.nan
        
        error_values.append(computed_error)
    
    # Create a DataArray for the error variable with dimension 'species'
    error_da = xr.DataArray(
        np.array(error_values),
        coords={'species': data.coords['species'].values},
        dims=['species'],
        name='error'
    )
    
    # Ensure the input DataArray is named 'Concentration_mM'
    if data.name is None:
        data = data.rename('Concentration_mM')
    
    # Convert the original data array to a dataset.
    ds = data.to_dataset()
    # Add the error variable.
    ds['error'] = error_da
    
    return ds