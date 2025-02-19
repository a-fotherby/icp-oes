import xarray as xr
import numpy as np

def find_dilutions(ds: xr.Dataset, variable_name: str, reference_sample: str = "H") -> dict:
    """
    For each sample (other than the reference sample 'H') and for each species,
    this function checks if the value for the species in the reference sample is smaller
    than that in the sample being compared. If so, it computes the dilution ratio:
    
        ratio = (sample value) / (reference sample value)
    
    For each sample, the function then identifies:
      - The species with the greatest ratio (the calibration species)
      - All other species that also have values greater than in the reference sample
    
    The function prints a formatted output for each sample with:
      - The sample name, calibration species, and its dilution ratio
      - A list of "Other species outside of calibration"
    
    It also returns a dictionary mapping each sample name (other than the reference) to a tuple:
              (calibration_species, max_ratio)
    
    Parameters:
        ds (xr.Dataset): An xarray.Dataset containing a data variable with dimensions
                         'species' and 'sample_name'.
        variable_name (str): The name of the data variable in the dataset to operate on.
        reference_sample (str): The label in the 'sample_name' coordinate to compare against.
                                Defaults to "H".
    
    Returns:
        dict: A dictionary mapping each sample name (other than the reference) to a tuple:
              (calibration_species, max_ratio).
    """
    results = {}
    # This dictionary will store, for each sample, a list of species that qualify (i.e. species with values greater than H)
    other_species_dict = {}
    
    # Extract the data variable from the dataset.
    data = ds[variable_name]
    
    # Get the list of species and sample names from the data variable's coordinates.
    species_list = data.coords['species'].values
    sample_names = data.coords['sample_name'].values

    # Loop over every sample except the reference sample.
    for sample in sample_names:
        if sample == reference_sample:
            continue  # Skip the reference sample
        
        max_ratio = None
        calibration_species = None
        qualifying_species = []  # All species that have sample value > H
        
        # Compare every species for this sample against the reference sample.
        for sp in species_list:
            sample_val = data.loc[dict(species=sp, sample_name=sample)].item()
            ref_val = data.loc[dict(species=sp, sample_name=reference_sample)].item()
            
            # Check the condition (avoid division by zero)
            if ref_val != 0 and sample_val > ref_val:
                ratio = sample_val / ref_val
                qualifying_species.append(sp)
                if max_ratio is None or ratio > max_ratio:
                    max_ratio = ratio
                    calibration_species = sp
        
        # Record the results if at least one species qualifies.
        if calibration_species is not None:
            results[sample] = (calibration_species, max_ratio)
            # Exclude the calibration species from the list of qualifying species for the "other species" list.
            others = [sp for sp in qualifying_species if sp != calibration_species]
            other_species_dict[sample] = others

    # Print the formatted output for each sample.
    for sample, (cal_species, ratio) in results.items():
        print('Sample:', sample)
        print(f'Largest species outside calibration: {cal_species}, Required dilution ratio: {ratio:.2f}')
        others = other_species_dict.get(sample, [])
        print(f"Other species outside of calibration: {others}")
        print('-----------------------------------')
    
    return results