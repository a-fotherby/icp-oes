import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import math

def bar_by_sample(avg_ds, comp_sample1="SLRS-6_1", comp_sample2="SPS-SW2 10%_1"):
    """
    Create a figure with subplots comparing species abundances for each sample in avg_ds.
    Each subplot is a grouped bar chart showing:
      - The current sample (left-shifted)
      - Comparison sample 1 (centered)
      - Comparison sample 2 (right-shifted)
    
    Parameters
    ----------
    avg_ds : xarray.Dataset
        The dataset containing your data. It must have:
          - A coordinate 'sample_name' (used for looping over samples)
          - A coordinate 'species' (used for x-axis labels)
        The data variable(s) are assumed to be directly accessible via .values.
    comp_sample1 : str, optional
        The name of the first comparison sample (default is "SLRS-6_1").
    comp_sample2 : str, optional
        The name of the second comparison sample (default is "SPS-SW2 10%_1").
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the subplots.
    """
    # Get all sample names from the dataset
    sample_names = avg_ds.sample_name.values
    n_samples = len(sample_names)
    
    # Create a square-ish grid for subplots
    ncols = int(np.ceil(np.sqrt(n_samples)))
    nrows = int(np.ceil(n_samples / ncols))
    
    # Create the figure and an array of axes
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
    
    # Flatten axes array for easy iteration (if only one subplot, make it a list)
    if n_samples == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Extract species names and create x positions for the bars
    species = avg_ds.species.values
    x = np.arange(len(species))
    bar_width = 0.25  # width for each bar
    
    # Pre-extract the comparison samples' data
    data_comp1 = avg_ds.sel(sample_name=comp_sample1).values
    data_comp2 = avg_ds.sel(sample_name=comp_sample2).values
    
    # Loop over each sample and plot its grouped bar chart in a subplot
    for idx, sample in enumerate(sample_names):
        ax = axes[idx]
        
        # Extract the data for the current sample
        curr_data = avg_ds.sel(sample_name=sample).values
        
        # Plot three sets of bars:
        # - current sample: left-shifted
        # - comparison sample 1: centered
        # - comparison sample 2: right-shifted
        ax.bar(x - bar_width, curr_data, width=bar_width, label=str(sample), color='skyblue')
        ax.bar(x, data_comp1, width=bar_width, label=comp_sample1, color='salmon')
        ax.bar(x + bar_width, data_comp2, width=bar_width, label=comp_sample2, color='limegreen')
        
        # Customize each subplot
        ax.set_xlabel("Species")
        ax.set_ylabel("Concentration / mM")
        ax.set_title(f"Comparison for Sample: {sample}")
        ax.set_xticks(x)
        ax.set_xticklabels(species, rotation=45)
        ax.legend()
    
    # Remove any unused subplots if the grid is larger than the number of samples
    for j in range(idx + 1, nrows * ncols):
        fig.delaxes(axes[j])
    
    fig.tight_layout()
    return fig


def bar_by_species(ds, var_name=None, highlight_samples=None):
    """
    Plot bar charts for each species from an xarray dataset and return the figure.
    
    Parameters
    ----------
    ds : xarray.Dataset
        A dataset with two coordinates: 'sample_name' and 'species'. It should contain
        a data variable with dimensions ('species', 'sample_name').
    var_name : str, optional
        The name of the variable in ds to plot. If not provided, the first data variable is used.
    highlight_samples : list of str, optional
        List of sample_name values to highlight. Bars corresponding to these sample names
        will be drawn in a different color.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the subplots.
    """
    # If no variable name is provided, pick the first data variable.
    if var_name is None:
        var_name = list(ds.data_vars.keys())[0]
    
    # If no samples to highlight are provided, default to an empty list.
    if highlight_samples is None:
        highlight_samples = []
    
    # Get the list of species.
    species_list = ds.coords['species'].values
    n_species = len(species_list)
    
    # Create one subplot per species.
    fig, axes = plt.subplots(n_species, 1, figsize=(10, 4 * n_species), squeeze=False)
    
    for i, species in enumerate(species_list):
        ax = axes[i, 0]
        
        # Select data for the current species.
        data = ds.sel(species=species)
        
        # Extract sample names and their corresponding values.
        sample_names = data.coords['sample_name'].values
        values = data.values
        
        # Set bar colors: highlight the bar if its sample_name is in the list.
        colors = ['red' if sample in highlight_samples else 'blue' for sample in sample_names]
        
        # Create the bar chart.
        ax.bar(sample_names, values, color=colors)
        ax.set_title(f"Species: {species}")
        ax.set_xlabel("Sample Name")
        ax.set_ylabel(var_name)
        ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    return fig