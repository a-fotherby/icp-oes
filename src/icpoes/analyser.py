#!/usr/bin/env python3
import argparse
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import re
import difflib
from io import StringIO

class ConcentrationAnalyser:
    def __init__(self, dataset: xr.Dataset = None):
        """
        Initialize with an xarray.Dataset.
        The dataset is expected to have:
          - A variable "Concentration" (in PPM)
          - Coordinates "sample_name" and "species"
        """
        self.raw = dataset.copy() if dataset is not None else None
        # For analysis, we use a lean dataset containing only "Concentration".
        self.ds = dataset[['Concentration']] if dataset is not None else None
        # Dictionary to store downstream results.
        self.results = {}

    @classmethod
    def from_csv(cls, file_path: str):
        """
        Import an ESWS CSV file, perform cleaning and conversion,
        and return a ConcentrationAnalyser instance.
        """
        ds = cls._import_raw(file_path)
        # For analysis we keep only the Concentration variable.
        ds = ds[['Concentration']]
        return cls(dataset=ds)

    # ------------------------
    # Data Importing Functions (unchanged)
    # ------------------------
    @staticmethod
    def _load_esws_csv(file_path: str) -> pd.DataFrame:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        header_index = next((i for i, line in enumerate(lines) if line.startswith('"Label",')), None)
        if header_index is None:
            raise ValueError("Header row not found in the file.")
        csv_data = ''.join(lines[header_index:])
        return pd.read_csv(StringIO(csv_data))

    @staticmethod
    def _clean_data(df: pd.DataFrame) -> pd.DataFrame:
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

    @staticmethod
    def _combine_element_labels(df: pd.DataFrame) -> pd.DataFrame:
        for col in ["Element", "Element Label"]:
            if col not in df.columns:
                raise ValueError(f"DataFrame must contain '{col}' column.")
        def combine_row(row):
            element = row["Element"]
            element_label = row["Element Label"]
            symbol = element.split()[0] if isinstance(element, str) else ""
            if isinstance(element_label, str) and element_label.startswith(symbol):
                remainder = element_label[len(symbol):].strip()
            else:
                remainder = element_label
            return f"{element} {remainder}".strip() if remainder else element
        df["Combined Element"] = df.apply(combine_row, axis=1)
        return df.drop(columns=["Element", "Element Label"])

    @staticmethod
    def _to_xarray_dataset(df: pd.DataFrame) -> xr.Dataset:
        if "Label" not in df.columns or "Date Time" not in df.columns:
            raise ValueError("DataFrame must contain 'Label' and 'Date Time' columns.")
        df['UniqueLabel'] = df['Label'] + "_" + df['Date Time'].dt.strftime('%Y%m%d%H%M%S')
        if "Combined Element" not in df.columns:
            raise ValueError("DataFrame must contain 'Combined Element' column.")
        combined_elements = sorted(df["Combined Element"].dropna().unique())
        unique_labels = sorted(df["UniqueLabel"].unique())
        id_vars = {"UniqueLabel", "Combined Element", "Label", "Rack:Tube", "Type", "Date Time"}
        measurement_columns = [col for col in df.columns if col not in id_vars]
        ds = xr.Dataset()
        for var in measurement_columns:
            pivot = df.pivot_table(index="UniqueLabel", columns="Combined Element", values=var, aggfunc='first')
            pivot = pivot.reindex(index=unique_labels, columns=combined_elements)
            ds[var] = (("label", "combined_element"), pivot.values)
        ds = ds.assign_coords(label=unique_labels, combined_element=combined_elements)
        return ds

    @staticmethod
    def _update_label_names(ds: xr.Dataset) -> xr.Dataset:
        original_labels = list(ds.coords['label'].data)
        bases = [lbl.rsplit('_', 1)[0] if '_' in lbl else lbl for lbl in original_labels]
        counts = {}
        for base in bases:
            counts[base] = counts.get(base, 0) + 1
        occurrence = {}
        new_labels = []
        for base in bases:
            occurrence[base] = occurrence.get(base, 0) + 1
            if counts[base] > 1:
                new_labels.append(f"{base}_{occurrence[base]}")
            else:
                new_labels.append(base)
        return ds.assign_coords(label=new_labels)

    @staticmethod
    def _remove_extraneous_variables(ds: xr.Dataset) -> xr.Dataset:
        vars_to_remove = [
            "Concentration SD", "Concentration % RSD", "Weight", "Volume", "Dilution",
            "Internal Standard", "Background correction", "%RSE limit", "LIMS ID",
            "Customer", "Customer Ref.", "Sample Date Time", "Sample Site", "Description"
        ]
        return ds.drop_vars(vars_to_remove, errors='ignore')

    @staticmethod
    def _move_unit_to_concentration(ds: xr.Dataset) -> xr.Dataset:
        if "Unit" in ds.data_vars and "Concentration" in ds.data_vars:
            unit_vals = ds["Unit"].values
            unit_vals_flat = unit_vals.flatten()
            unique_units = pd.unique(unit_vals_flat)
            unique_units = [u for u in unique_units if pd.notna(u)]
            if len(unique_units) == 1:
                ds["Concentration"].attrs["units"] = unique_units[0]
            elif len(unique_units) > 1:
                ds["Concentration"].attrs["units"] = ','.join(map(str, unique_units))
            ds = ds.drop_vars("Unit")
        return ds

    @staticmethod
    def _sanitize_netcdf_names(ds: xr.Dataset) -> xr.Dataset:
        def sanitize(name):
            new_name = re.sub(r'[^A-Za-z0-9_]', '_', name)
            if not re.match(r'^[A-Za-z_]', new_name):
                new_name = '_' + new_name
            return new_name
        data_var_mapping = {var: sanitize(var) for var in ds.data_vars}
        ds = ds.rename(data_var_mapping)
        coord_mapping = {coord: sanitize(coord) for coord in ds.coords}
        ds = ds.rename(coord_mapping)
        return ds

    @classmethod
    def _import_raw(cls, file_path: str) -> xr.Dataset:
        df = cls._load_esws_csv(file_path)
        df = cls._clean_data(df)
        df = cls._combine_element_labels(df)
        ds = cls._to_xarray_dataset(df)
        ds = cls._update_label_names(ds)
        ds = ds.rename({'label': 'sample_name', 'combined_element': 'species'})
        ds = cls._remove_extraneous_variables(ds)
        ds = cls._move_unit_to_concentration(ds)
        ds = cls._sanitize_netcdf_names(ds)
        return ds

    # ------------------------
    # Error Analysis Pipeline (working in PPM until final conversion)
    # ------------------------
    def drop_extraneous_vars(self):
        """
        Ensure the dataset contains only the Concentration variable.
        (The from_csv method already limits the dataset, but this is provided for pipeline clarity.)
        """
        self.ds = self.ds[['Concentration']]
        return self

    
    def compute_calibration_stats(self, check_std_names):
        """
        Compute calibration statistics (mean and SD in PPM) for the calibration samples identified 
        by the provided check standard stem names. Sample names are assumed to be of the form 
        "{check_std}_{number}". For each check_std, all samples with sample_name starting with 
        "{check_std}_" are selected, and the mean and standard deviation are computed (collapsing 
        the sample_name dimension). The results are stored in a calibration statistics dataset that 
        has a new dimension "check_std" with coordinates corresponding to each check standard.
        
        The resulting dataset contains:
        - calibration_mean_ppm (dimensions: check_std, species)
        - calibration_sd_ppm   (dimensions: check_std, species)
        
        Parameters
        ----------
        check_std_names : str or list of str
            One or more check standard names (e.g. "SLRS-6") used to identify calibration samples.
        
        Returns
        -------
        self
            The instance with self.results['calibration_stats'] updated.
        """
        # Ensure we have a list of check standard names.
        if isinstance(check_std_names, str):
            check_std_names = [check_std_names]
        
        # Get the full sample_name coordinate from the Concentration DataArray.
        all_samples = self.ds["Concentration"].coords["sample_name"].values.astype(str)
        
        calib_stats_list = []
        
        # Loop over each check standard name.
        for cs in check_std_names:
            # Identify samples whose name starts with "{cs}_"
            mask = np.char.startswith(all_samples, cs + "_")
            selected_samples = all_samples[mask]
            if len(selected_samples) == 0:
                print(f"Warning: No samples found for check standard '{cs}'.")
                continue
            
            # Select the calibration subset.
            subset = self.ds["Concentration"].sel(sample_name=selected_samples)
            
            # Optionally, warn if some species have missing values.
            total_count = len(selected_samples)
            valid_counts = subset.count(dim='sample_name')
            missing_mask = valid_counts < total_count
            if missing_mask.any():
                missing_species = valid_counts.where(missing_mask, drop=True).coords['species'].values
                print(f"Warning for {cs}: The following species have missing values in the calibration subset: {missing_species}")
            
            # Compute mean and standard deviation (in PPM) over the calibration subset.
            calibration_mean = subset.mean(dim='sample_name', skipna=True)
            calibration_sd = subset.std(dim='sample_name', skipna=True)
            
            # Create a temporary dataset for this check standard.
            ds_cs = xr.Dataset({
                'calibration_mean_ppm': calibration_mean,
                'calibration_sd_ppm': calibration_sd
            }, coords={'species': calibration_mean.coords['species']})
            
            # Expand dims to add a new dimension "check_std" with coordinate value cs.
            ds_cs = ds_cs.expand_dims(check_std=[cs])
            calib_stats_list.append(ds_cs)
        
        if len(calib_stats_list) == 0:
            raise ValueError("No calibration statistics computed: none of the provided check standards matched any samples.")
        
        # Concatenate all calibration statistics along the new "check_std" dimension.
        calib_ds = xr.concat(calib_stats_list, dim='check_std')
        self.results['calibration_stats'] = calib_ds
        return self
    

    def compute_error_estimates(self, certified_values: dict):
        """
        Using the calibration_stats dataset (dimensions: check_std, species) computed from the calibration subset,
        for each check_std and species, compute the percentage difference between the calibration mean (in PPM) and the
        certified value. The certified value is obtained from a dictionary of dictionaries, where the top-level keys are
        the check_std values (matching those in calibration_stats) and the inner dictionary maps element symbols to
        certified values (in PPM). Then compute the error for each (check_std, species) as the maximum of:
            - (error from percentage difference) and
            - (2 * calibration SD).
        However, only check_std values where the absolute percentage difference is less than 5% are considered.
        For each species, take the maximum error (across the check_std dimension) among those that meet the criterion,
        and broadcast that error to all samples. The final error remains in PPM.
        
        Parameters
        ----------
        certified_values : dict
            A dictionary of dictionaries. For each check_std key (e.g. "SLRS-6"), the value is another dictionary
            mapping element symbols (e.g. "Al", "Ba", "Fe") to the certified value (in PPM).
        
        Returns
        -------
        self
        """
        calib_ds = self.results.get('calibration_stats')
        if calib_ds is None:
            raise ValueError("Calibration statistics not computed. Call compute_calibration_stats() first.")
        
        # Get the list of calibration standards (check_std) and species.
        cs_list = calib_ds.coords['check_std'].values
        species = calib_ds.coords['species'].values
        
        # Build a DataArray of certified values, shape (len(cs_list), len(species)).
        # For each check_std (cs) and species, extract the element and lookup the certified value.
        cert_array = np.empty((len(cs_list), len(species)))
        for i, cs in enumerate(cs_list):
            for j, sp in enumerate(species):
                element = sp.split()[0]
                cert_array[i, j] = certified_values.get(cs, {})[0].get(element, np.nan)
        certified_ppm_da = xr.DataArray(cert_array, coords={'check_std': cs_list, 'species': species},
                                        dims=['check_std', 'species'])
        
        # Compute percentage difference for each (check_std, species).
        pct_diff = (calib_ds['calibration_mean_ppm'] - certified_ppm_da) / certified_ppm_da * 100

        # Compute the error from the percentage difference.
        error_from_pct = (abs(pct_diff) / 100) * calib_ds['calibration_mean_ppm']
        
        # For each (check_std, species), compute the error as the maximum of error_from_pct and 2 * calibration_sd_ppm.
        error_per_cs = xr.apply_ufunc(np.fmax, error_from_pct, 2 * calib_ds['calibration_sd_ppm'])
        
        # Only consider check_std values where the absolute percentage difference is less than 5%.
        valid_mask = abs(pct_diff) < 10
        error_per_cs = error_per_cs.where(valid_mask, np.nan)
        
        # For each species, take the maximum error across check_std (ignoring NaNs).
        final_error_ppm = error_per_cs.max(dim="check_std", skipna=True)
        
        # Broadcast the final error (per species) to all samples.
        all_samples = self.ds["Concentration"].coords["sample_name"].values
        final_error_broadcast = final_error_ppm.expand_dims({"sample_name": all_samples})
        final_error_broadcast = final_error_broadcast.broadcast_like(self.ds["Concentration"])
        
        # Save the final calibration error (in PPM) in the calibration stats dataset.
        self.results['calibration_stats'] = calib_ds.assign(calibration_error_ppm_max=final_error_broadcast)
        
        return self

    def aggregate_by_element(self):
        """
        Aggregate the unaggregated results by averaging over wavelengths (species) that share the same element.
        Only wavelengths with non-NaN error values in all samples (i.e. passing quality control) are included.
        The element is determined by taking the first token from each species string.
        
        This method assumes that self.results["unaggregated"] is a dataset with:
        - Variables "Concentration_ppm" and "error"
        - Coordinates "sample_name" and "species"
        
        It returns a new dataset with dimensions ("sample_name", "element") where the variables are averaged
        over species belonging to each element. In addition, a dictionary mapping each element to the list
        of species that went into the average is attached as an attribute "contributing_wavelengths" on the
        "Concentration_ppm" variable.
        """
        ds = self.results.get("unaggregated")
        if ds is None:
            raise ValueError("Unaggregated results not available. Ensure self.results['unaggregated'] is set.")
        
        # Determine valid species: those with non-NaN error across all samples.
        valid_mask = ds["error"].notnull().all(dim="sample_name")
        valid_species = ds.coords["species"].values[valid_mask.values]
        
        # Determine the full set of elements (from all species) and valid elements.
        all_species = ds.coords["species"].values
        all_elements = {sp.split()[0] for sp in all_species}
        valid_elements = {sp.split()[0] for sp in valid_species}
        
        # Warn for any element that has no wavelengths passing quality control.
        for elem in all_elements:
            if elem not in valid_elements:
                print(f"Warning: For element '{elem}', no wavelengths pass quality control.")
        
        # Select only the valid species.
        ds_valid = ds.sel(species=valid_species)
        
        # Create a new coordinate "element" by extracting the first token from each valid species name.
        element_coords = [sp.split()[0] for sp in ds_valid.coords["species"].values]
        ds_valid = ds_valid.assign_coords(element=("species", element_coords))
        
        # Build a dictionary mapping each element to the list of species that are valid.
        contributing = {}
        for sp, elem in zip(ds_valid.coords["species"].values, ds_valid.coords["element"].values):
            contributing.setdefault(elem, []).append(sp)
        # Convert the lists to strings for storage as an attribute. 
        contributing = str(contributing)

        # Group by the "element" coordinate and average over the species dimension.
        aggregated_ds = ds_valid.groupby("element").mean(dim="species")
        
        # Attach the contributing wavelengths as an attribute on the "Concentration_ppm" variable.
        if "Concentration_ppm" in aggregated_ds:
            aggregated_ds["Concentration_ppm"].attrs["contributing_wavelengths"] = contributing
        else:
            aggregated_ds.attrs["contributing_wavelengths"] = contributing
        
        self.results["aggregated"] = aggregated_ds
        return self



    def get_results(self):
        """
        Combine the raw concentration measurements (in ppm) with the calculated maximum calibration error 
        (in ppm) into a final results dataset. This method assumes that the calibration error has already been 
        computed and broadcast in self.results['calibration_stats'] under the variable 'calibration_error_ppm_max'.
        
        Returns
        -------
        xr.Dataset
            A dataset with dimensions ("sample_name", "species") containing:
            - "Concentration": raw concentration measurements in ppm.
            - "error": the calibration error (in ppm) for each species applied to all samples.
        """
        calib_stats = self.results.get('calibration_stats')
        if calib_stats is None or 'calibration_error_ppm_max' not in calib_stats:
            raise ValueError("Calibration error not computed. Run compute_error_estimates() first.")
        
        conc_ppm = self.ds["Concentration"]
        error_ppm = calib_stats["calibration_error_ppm_max"]
        
        final_ds = xr.Dataset(
            {
                "Concentration_ppm": (("sample_name", "species"), conc_ppm.values),
                "error": (("sample_name", "species"), error_ppm.values)
            },
            coords={
                "sample_name": conc_ppm.coords["sample_name"].values,
                "species": conc_ppm.coords["species"].values
            }
        )
        
        self.results['unaggregated'] = final_ds
        return self


    def save_data(self, file_path: str):
        """
        Save the final results dataset to a NetCDF file.
        """
        final_ds = self.results.get('aggregated')
        if final_ds is None:
            raise ValueError("No results to save. Run the error analysis pipeline first.")
        final_ds.to_netcdf(file_path)
        print(f"Results saved to {file_path}")
        return self


# ------------------------
# Main Entry Point
# ------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyser for spectrometer data with error analysis pipeline")
    parser.add_argument("file_path", help="Path to the input CSV file")
    parser.add_argument("output_name", help="Output NetCDF file name")
    args = parser.parse_args()
    
    # Create an analyser instance from the CSV file (data in PPM).
    analyser = ConcentrationAnalyser.from_csv(args.file_path)
    
    # Run the error analysis pipeline:
    # 1. Compute calibration statistics using only the provided subset of samples.
    # 2. Compute error estimates based on the calibration statistics relative to certified values.
    # 3. Convert the calibration stats (mean, SD, error) from PPM to mM.
    # Here we assume the calibration subset is ["sample1", "sample2", "sample3"].
    analyser.drop_extraneous_vars()\
            .compute_calibration_stats(sample_subset=["sample1", "sample2", "sample3"])\
            .compute_error_estimates(certified_values_mM={"Al": 1.25e-3, "Ba": 1.04e-4, "Fe": 1.51e-3},
                                     conversion_factor={"Al": 26.98, "Ba": 137.33, "Fe": 55.85})\
            .final_conversion(conversion_factor={"Al": 26.98, "Ba": 137.33, "Fe": 55.85})
    
    # Retrieve final results and save them.
    results_ds = analyser.get_results()
    results_ds.to_netcdf(args.output_name)
    print(f"Error analysis complete. Results saved to {args.output_name}")
