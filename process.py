#!/usr/bin/env python3
import argparse
from import_icpoes import import_raw
from unit_conv import convert_units

def main():
    parser = argparse.ArgumentParser(
        description="Import raw data, convert units, and save as a NetCDF file."
    )
    parser.add_argument(
        "file_path", 
        help="Path to the input CSV file (e.g., /path/to/data.csv)"
    )
    parser.add_argument(
        "output_name", 
        help="Name of the output NetCDF file (e.g., output.nc)"
    )
    args = parser.parse_args()

    # Import the raw dataset from the file
    ds = import_raw(args.file_path)
    # Convert the dataset units
    ds = convert_units(ds)
    # Save the dataset to a NetCDF file
    ds.to_netcdf(args.output_name)
    print(f"Data processed and saved to {args.output_name}")

if __name__ == "__main__":
    main()
