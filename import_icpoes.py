#!/usr/bin/env python3
import argparse
import pandas as pd
import xarray as xr
from io import StringIO
import re  # Needed for sanitization

def load_esws_csv(file_path):
    """
    Load an ESWS CSV file into a DataFrame.
    Skips initial metadata by finding the first line that starts with '"Label",'.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    header_index = next((i for i, line in enumerate(lines) if line.startswith('"Label",')), None)
    if header_index is None:
        raise ValueError("Header row not found in the file.")
    csv_data = ''.join(lines[header_index:])
    return pd.read_csv(StringIO(csv_data))

def clean_data(df):
    """
    Clean the DataFrame:
      - Trim column names.
      - Replace '-' and 'N/A' with pd.NA.
      - Convert "Date Time" (and "Sample Date Time" if present) to datetime.
      - Coerce specified columns to numeric.
      - Drop columns that are entirely NA.
    """
    df.columns = df.columns.str.strip()
    df = df.replace({'-': pd.NA, 'N/A': pd.NA})
    if "Date Time" in df.columns:
        df["Date Time"] = pd.to_datetime(df["Date Time"], errors='coerce')
    if "Sample Date Time" in df.columns:
        df["Sample Date Time"] = pd.to_datetime(df["Sample Date Time"], errors='coerce')
    
    numeric_columns = [
        "Unadjusted Data", "Concentration", "Intensity", "Concentration SD", 
        "Concentration % RSD", "Intensity SD", "Intensity % RSD", "%RSE", "Weight", 
        "Volume", "Dilution", "Correlation coefficient", "%RSE limit", 
        "Calibration Coefficient [1]", "Calibration Coefficient [2]", "Calibration Coefficient [3]", 
        "Replicates", "Concentration Replicate 1", "Concentration Replicate 2", 
        "Concentration Replicate 3", "Concentration Replicate 4", "Concentration Replicate 5", 
        "Intensity Replicate 1", "Intensity Replicate 2", "Intensity Replicate 3", 
        "Intensity Replicate 4", "Intensity Replicate 5"
    ]
    for col in numeric_columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df.dropna(axis=1, how='all')

def combine_element_labels(df):
    """
    Combine the 'Element' and 'Element Label' columns into a new 'Combined Element' column.
    Example: "Ba 455.403" and "Ba Ax" -> "Ba 455.403 Ax".
    Then drop the original 'Element' and 'Element Label' columns.
    """
    for col in ["Element", "Element Label"]:
        if col not in df.columns:
            raise ValueError(f"DataFrame must contain '{col}' column.")
    
    def combine_row(row):
        element = row["Element"]
        element_label = row["Element Label"]
        # Use the element symbol (first token) to remove duplicate information.
        symbol = element.split()[0] if isinstance(element, str) else ""
        if isinstance(element_label, str) and element_label.startswith(symbol):
            remainder = element_label[len(symbol):].strip()
        else:
            remainder = element_label
        return f"{element} {remainder}".strip() if remainder else element

    df["Combined Element"] = df.apply(combine_row, axis=1)
    return df.drop(columns=["Element", "Element Label"])

def to_xarray_dataset(df):
    """
    Convert the DataFrame into an xarray.Dataset.
      - Create a unique sample identifier ("UniqueLabel") by appending the timestamp
        (formatted as YYYYMMDDHHMMSS) to the original "Label".
      - The 'label' dimension will be the sorted list of UniqueLabel.
      - The 'combined_element' dimension is the sorted list of unique "Combined Element" names.
      - Measurement columns (i.e. non-identifier columns) are pivoted into data variables.
    """
    if "Label" not in df.columns or "Date Time" not in df.columns:
        raise ValueError("DataFrame must contain 'Label' and 'Date Time' columns.")
    df['UniqueLabel'] = df['Label'] + "_" + df['Date Time'].dt.strftime('%Y%m%d%H%M%S')
    
    if "Combined Element" not in df.columns:
        raise ValueError("DataFrame must contain 'Combined Element' column.")
    combined_elements = sorted(df["Combined Element"].dropna().unique())
    unique_labels = sorted(df["UniqueLabel"].unique())
    
    # Identify measurement columns (exclude identifier columns)
    id_vars = {"UniqueLabel", "Combined Element", "Label", "Rack:Tube", "Type", "Date Time"}
    measurement_columns = [col for col in df.columns if col not in id_vars]
    
    ds = xr.Dataset()
    for var in measurement_columns:
        pivot = df.pivot_table(index="UniqueLabel", columns="Combined Element", values=var, aggfunc='first')
        pivot = pivot.reindex(index=unique_labels, columns=combined_elements)
        ds[var] = (("label", "combined_element"), pivot.values)
    ds = ds.assign_coords(label=unique_labels, combined_element=combined_elements)
    return ds

def update_label_names(ds):
    """
    Update the 'label' coordinate:
      - Remove the appended timestamp (split on the last underscore).
      - For duplicate base names, append a sequential suffix (_1, _2, etc.).
    """
    original_labels = list(ds.coords['label'].data)
    # Remove timestamp by splitting on the last underscore.
    bases = [lbl.rsplit('_', 1)[0] if '_' in lbl else lbl for lbl in original_labels]
    
    # Count occurrences of each base.
    counts = {}
    for base in bases:
        counts[base] = counts.get(base, 0) + 1
    
    # Disambiguate duplicate base names by appending a number.
    occurrence = {}
    new_labels = []
    for base in bases:
        occurrence[base] = occurrence.get(base, 0) + 1
        if counts[base] > 1:
            new_labels.append(f"{base}_{occurrence[base]}")
        else:
            new_labels.append(base)
    return ds.assign_coords(label=new_labels)

def remove_extraneous_variables(ds):
    """
    Remove extraneous variables from the dataset.
    Variables removed (if present):
      Concentration SD, Concentration % RSD, Weight, Volume, Dilution,
      Internal Standard, Background correction, %RSE limit, LIMS ID,
      Customer, Customer Ref., Sample Date Time, Sample Site, Description.
    """
    vars_to_remove = [
        "Concentration SD", "Concentration % RSD", "Weight", "Volume", "Dilution",
        "Internal Standard", "Background correction", "%RSE limit", "LIMS ID",
        "Customer", "Customer Ref.", "Sample Date Time", "Sample Site", "Description"
    ]
    return ds.drop_vars(vars_to_remove, errors='ignore')

def move_unit_to_concentration(ds):
    """
    If a variable called 'Unit' exists in the dataset, extract its unique value(s) and add it as a 
    'units' attribute to the 'Concentration' variable. Then remove the 'Unit' variable from the dataset.
    """
    if "Unit" in ds.data_vars and "Concentration" in ds.data_vars:
        # Extract unit values and flatten the array.
        unit_vals = ds["Unit"].values
        unit_vals_flat = unit_vals.flatten()
        # Get unique non-null unit values.
        unique_units = pd.unique(unit_vals_flat)
        unique_units = [u for u in unique_units if pd.notna(u)]
        if len(unique_units) == 1:
            ds["Concentration"].attrs["units"] = unique_units[0]
        elif len(unique_units) > 1:
            # In case of multiple unique units, join them into a comma-separated string.
            ds["Concentration"].attrs["units"] = ','.join(map(str, unique_units))
        # Drop the 'Unit' variable.
        ds = ds.drop_vars("Unit")
    return ds

def sanitize_netcdf_names(ds):
    """
    Sanitize variable and coordinate names in the xarray.Dataset to ensure compatibility with the netCDF format.
    The netCDF naming conventions typically require names to start with a letter or underscore and contain only
    letters, digits, or underscores.
    """
    def sanitize(name):
        # Replace any character that is not a letter, digit, or underscore with an underscore.
        new_name = re.sub(r'[^A-Za-z0-9_]', '_', name)
        # If the first character is not a letter or underscore, prefix with an underscore.
        if not re.match(r'^[A-Za-z_]', new_name):
            new_name = '_' + new_name
        return new_name

    # Sanitize data variable names.
    data_var_mapping = {var: sanitize(var) for var in ds.data_vars}
    ds = ds.rename(data_var_mapping)
    
    # Sanitize coordinate names.
    coord_mapping = {coord: sanitize(coord) for coord in ds.coords}
    ds = ds.rename(coord_mapping)
    
    return ds

def import_raw(file_path):
    """
    Main workflow:
      1. Load the CSV file.
      2. Clean the DataFrame.
      3. Combine element labels.
      4. Convert the DataFrame to an xarray.Dataset.
      5. Update the label names (remove timestamps and disambiguate duplicates).
      6. Rename coordinates: 'label' -> 'sample_name' and 'combined_element' -> 'species'.
      7. Remove extraneous variables.
      8. Move the 'Unit' variable's information into the metadata for 'Concentration'.
      9. Sanitize all variable and coordinate names to be netCDF-compatible.
    Returns the final xarray.Dataset.
    """
    df = load_esws_csv(file_path)
    df = clean_data(df)
    df = combine_element_labels(df)
    ds = to_xarray_dataset(df)
    ds = update_label_names(ds)
    ds = ds.rename({'label': 'sample_name', 'combined_element': 'species'})
    ds = remove_extraneous_variables(ds)
    ds = move_unit_to_concentration(ds)  # Move 'Unit' into 'Concentration' metadata and delete it.
    ds = sanitize_netcdf_names(ds)  # Sanitize variable names for netCDF compatibility.
    return ds