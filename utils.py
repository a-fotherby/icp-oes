import numpy as np
import xarray as xr

def select_vars_by_element(ds, element, dim='species'):
    """
    Subsets an xarray.Dataset by selecting only those entries along a given 
    coordinate dimension that match a specified chemical element.

    This function assumes that each coordinate value is a string with the 
    chemical element as its first token (e.g., "Al 396.152 Ax").

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset that contains the coordinate dimension.
    element : str
        The chemical element to filter by (e.g., 'Al', 'Mg', etc.).
    dim : str, optional
        The name of the coordinate dimension to filter on. Default is 'variable'.

    Returns
    -------
    xarray.Dataset
        A subset of the original dataset containing only the entries for which 
        the coordinate (in dimension `dim`) starts with the given chemical element.
    """
    # Extract the coordinate values from the dataset.
    coord_values = ds.coords[dim].values
    
    # Build a boolean mask: True if the first token of the coordinate matches the element.
    mask = np.array([val.split()[0] == element for val in coord_values])
    
    # If no coordinates match, you might want to raise an error or return an empty dataset.
    if not np.any(mask):
        raise ValueError(f"No coordinates found for chemical element '{element}'")
    
    # Use the boolean mask to filter the coordinate values.
    filtered_values = coord_values[mask]
    
    # Select and return the subset of the dataset using .sel()
    return ds.sel({dim: filtered_values})


    import xarray as xr


import difflib
import numpy as np

def select_by_fuzzy_sample(ds, sample_query, dim='sample_name', cutoff=0.6):
    """
    Fuzzy selects a subset of an xarray.Dataset based on a fuzzy match against 
    the values in the specified coordinate dimension (default 'sample_name').

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset containing the coordinate to be fuzzy matched.
    sample_query : str
        The query string to fuzzy match against the coordinate values.
    dim : str, optional
        The coordinate dimension to match on. Default is 'sample_name'.
    cutoff : float, optional
        A value between 0 and 1 that determines the similarity threshold for 
        matching. Higher values require a closer match. Default is 0.6.

    Returns
    -------
    xarray.Dataset
        A subset of the original dataset containing only those entries whose
        coordinate value in `dim` fuzzy matches the query.

    Raises
    ------
    ValueError
        If no coordinate values match the fuzzy query with the given cutoff.
    """
    # Ensure that coordinate values are strings.
    coord_values = ds.coords[dim].values.astype(str)
    
    # Use difflib.get_close_matches to find values that match the query.
    matches = difflib.get_close_matches(sample_query, coord_values, n=len(coord_values), cutoff=cutoff)
    
    if not matches:
        raise ValueError(f"No sample names fuzzy matched to query '{sample_query}' with cutoff {cutoff}")
    
    # Return the subset of the dataset with matching sample names.
    return ds.sel({dim: matches})


def find_common_elements(d: dict) -> set:
    """
    Takes a dictionary whose values are all lists,
    and returns the set of elements common to all lists.
    
    Example:
    d = {
        'first': [1, 2, 3, 4],
        'second': [2, 3, 5],
        'third': [0, 2, 3, 6]
    }
    The function will return {2, 3}.
    """
    # Start by converting the first list in the dictionary to a set
    # and set that as our "running intersection".
    all_keys = list(d.keys())
    if not all_keys:
        return set()  # Edge case: if the dictionary is empty, return empty set
    
    # Convert the first list's values into a set to start the intersection
    common_set = set(d[all_keys[0]])
    
    # Intersect with the sets of the remaining lists
    for key in all_keys[1:]:
        common_set = common_set.intersection(d[key])
        # If at any point the common_set is empty, we can stop early
        if not common_set:
            return set()
    
    return common_set

def find_missing_elements(required_list, target_list):
    """
    Checks which elements in 'required_list' do not appear in 'target_list'.
    Returns a list of those missing elements (preserving duplicates if present).

    :param required_list: List of elements you require to be in 'target_list'.
    :param target_list:   List in which you check for the elements.
    :return: A list of elements from 'required_list' that are not in 'target_list'.
    """
    missing = []
    # Convert target_list to a set for O(1) lookups
    target_set = set(target_list)

    # Check each element of the required_list
    for elem in required_list:
        if elem not in target_set:
            missing.append(elem)

    return missing


import xarray as xr

def average_wavelengths(ds: xr.Dataset, species_set: set) -> xr.Dataset:
    """
    For each chemical element (determined by the first token of each species label)
    average the values of the "Concentration_mM" variable over all species that contain that
    element, and then add the result as a new species coordinate with a label of the form
    "{element}_average".
    
    Parameters:
      ds (xr.Dataset): Input dataset that contains:
                         - a coordinate "species" whose values are strings like "Al 396.152 R"
                         - a data variable "Concentration_mM" defined along the "species" dimension.
      species_set (set): A set of species names (strings) that are contained in ds.coords["species"].
    
    Returns:
      xr.Dataset: A new dataset with an updated "Concentration_mM" variable that includes additional
                  entries along the "species" dimension for each element average.
    """
    
    # Group species names by chemical element (assumed to be the first token of the string)
    element_to_species = {}
    for sp in species_set:
        element = sp.split()[0]  # e.g. "Al" from "Al 396.152 R"
        element_to_species.setdefault(element, []).append(sp)
    
    # List to store the DataArrays for each element's average, and the new species labels.
    avg_dataarrays = []
    new_species_labels = []
    
    for element, sp_list in element_to_species.items():
        # Select the data corresponding to the species in this group and average over that dimension.
        # Note: if Concentration_mM has additional dimensions (e.g. time), they are preserved.
        avg = ds["Concentration_mM"].sel(species=sp_list).mean(dim="species")
        # Because the "species" dimension was reduced by the mean(), we add it back so that
        # we can concatenate later. (This creates a new axis of length 1.)
        avg = avg.expand_dims("species")
        
        avg_dataarrays.append(avg)
        new_species_labels.append(f"{element}_average")
    
    # Concatenate all the new average DataArrays along the species dimension.
    avg_da = xr.concat(avg_dataarrays, dim="species")
    # Set the species coordinate for these new entries.
    avg_da = avg_da.assign_coords(species=new_species_labels)
    
    return avg_da


def match_sample_name(ds: xr.Dataset, labels, threshold=0.6) -> xr.Dataset:
    """
    Subset an xarray.Dataset by fuzzy matching a list of labels against the "sample_name" coordinate.

    Parameters:
    -----------
    ds : xr.Dataset
        The input dataset, which must contain a coordinate named "sample_name".
    labels : list of str
        A list of strings to fuzzy match against the "sample_name" coordinate.
    threshold : float, optional
        The minimum similarity ratio (between 0 and 1) required for a match.
        Defaults to 0.6.

    Returns:
    --------
    xr.Dataset
        A subset of the original dataset containing only the samples whose "sample_name"
        matches at least one of the provided labels based on the specified threshold.
    """
    if "sample_name" not in ds.coords:
        raise KeyError('The dataset does not have a coordinate named "sample_name".')
    
    # Ensure all sample names are treated as strings.
    sample_names = [str(s) for s in ds.coords["sample_name"].values]
    matching_samples = []

    # Loop over each sample name and check if it fuzzy-matches any label.
    for sample in sample_names:
        for label in labels:
            similarity = difflib.SequenceMatcher(None, sample, label).ratio()
            if similarity >= threshold:
                matching_samples.append(sample)
                break  # If one label matches, move to the next sample.

    if not matching_samples:
        print("Warning: No samples matched the provided labels with the given threshold.")
        # Return an empty subset along the 'sample_name' dimension.
        return ds.isel(sample_name=[])

    # Return the subset of the dataset based on the matching sample names.
    return ds.sel(sample_name=matching_samples)
