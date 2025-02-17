import numpy as np
import matplotlib.pyplot as plt
import itertools

def bar_by_sample(avg_ds, comp_samples=["SLRS-6_1", "SPS-SW2 10%_1"]):
    """
    Create a figure with subplots comparing species abundances for each sample in avg_ds.
    Each subplot is a grouped bar chart showing:
      - The current sample (plotted as the left-most bar)
      - Each comparison sample (in the order provided in comp_samples)
    Error bars for each species are applied using the values contained in the 'error'
    variable of the dataset.
    
    Parameters
    ----------
    avg_ds : xarray.Dataset
        The dataset containing your data. It must have:
          - A coordinate 'sample_name' (used for looping over samples) in the 'Concentration_mM' variable.
          - A coordinate 'species' (used for x-axis labels) in the 'Concentration_mM' variable.
          - Two variables:
              'Concentration_mM': containing concentration data.
              'error': containing error values for each species.
    comp_samples : list of str, optional
        A list of sample names to use for comparison (default is ["SLRS-6_1", "SPS-SW2 10%_1"]).
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the subplots.
    """
    # Access the concentration and error variables
    conc = avg_ds['Concentration_mM']
    errors = avg_ds['error']
    
    # Remove the error sample (if it exists) from the list of samples to plot
    sample_names = [s for s in conc.sample_name.values if s != "error"]
    n_samples = len(sample_names)
    
    # Extract error values (assumed to correspond to each species)
    error = errors.values
    
    # Create a square-ish grid for subplots
    ncols = int(np.ceil(np.sqrt(n_samples)))
    nrows = int(np.ceil(n_samples / ncols))
    
    # Create the figure and an array of axes
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
    
    # Flatten the axes array for easy iteration (if only one subplot, wrap it in a list)
    if n_samples == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Extract species names and create x positions for the bars
    species = conc.species.values
    x = np.arange(len(species))
    bar_width = 0.25  # width for each bar
    
    # Determine the total number of bars in each group:
    # current sample + len(comp_samples)
    total_bars = 1 + len(comp_samples)
    # Calculate offsets so that the left-most bar (current sample) is at:
    #   x + ( - (total_bars - 1)/2 * bar_width )
    # and subsequent bars are spaced by bar_width.
    offsets = [ -((total_bars - 1) / 2) * bar_width + i * bar_width for i in range(total_bars) ]
    
    # Define colors for the bars.
    # We'll use a fixed color for the current sample and cycle through colors for comparisons.
    current_color = 'skyblue'
    comp_color_cycle = itertools.cycle(['salmon', 'limegreen', 'orange', 'violet'])
    
    # Loop over each sample (excluding the error sample) and plot its grouped bar chart
    for idx, sample in enumerate(sample_names):
        ax = axes[idx]
        
        # Extract the concentration data for the current sample
        curr_data = conc.sel(sample_name=sample).values
        
        # Plot the current sample bar at its offset
        ax.bar(x + offsets[0], curr_data, width=bar_width, label=str(sample),
               color=current_color, yerr=error, capsize=3)
        
        # Plot each comparison sample at its respective offset
        for i, comp_sample in enumerate(comp_samples):
            comp_data = conc.sel(sample_name=comp_sample).values
            color = next(comp_color_cycle)
            ax.bar(x + offsets[i+1], comp_data, width=bar_width, label=comp_sample,
                   color=color, yerr=error, capsize=3)
        # Reset the color cycle for the next subplot
        comp_color_cycle = itertools.cycle(['salmon', 'limegreen', 'orange', 'violet'])
        
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


def bar_by_species(ds, var_name="Concentration_mM", comp_samples=None):
    """
    Plot bar charts for each species from an xarray dataset and return the figure.
    
    This function expects the dataset to have two coordinates: 'sample_name' and 'species', and
    two variables:
      - 'Concentration_mM': the concentration data with dimensions ('species', 'sample_name')
      - 'error': the error values (with dimensions matching 'Concentration_mM')
    
    Parameters
    ----------
    ds : xarray.Dataset
        A dataset with coordinates 'sample_name' and 'species' and the variables 'Concentration_mM'
        and 'error'.
    var_name : str, optional
        The name of the variable in ds to plot for concentrations. Default is 'Concentration_mM'.
    comp_samples : list of str, optional
        List of sample_name values to highlight as comparison samples. Bars corresponding to these
        sample names will be drawn in a different color. If not provided, no samples are highlighted.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the subplots.
    """
    # If no comparison samples are provided, default to an empty list.
    if comp_samples is None:
        comp_samples = []
    
    # Get the list of species from the dataset coordinates.
    species_list = ds.coords['species'].values
    n_species = len(species_list)
    
    # Create one subplot per species.
    fig, axes = plt.subplots(n_species, 1, figsize=(10, 4 * n_species), squeeze=False)
    
    for i, species in enumerate(species_list):
        ax = axes[i, 0]
        
        # Select the concentration data for the current species.
        # Data is assumed to have dimensions ('sample_name',) after selecting a specific species.
        data = ds[var_name].sel(species=species)
        # Also select the error data for the current species.
        error_data = ds['error'].sel(species=species).values
        
        # Extract the sample names and their corresponding values.
        sample_names = data.coords['sample_name'].values
        values = data.values
        
        # Determine bar colors: if a sample is in comp_samples, highlight it in red; otherwise, use blue.
        colors = ['red' if sample in comp_samples else 'blue' for sample in sample_names]
        
        # Create the bar chart for the current species with error bars.
        ax.bar(sample_names, values, color=colors, yerr=error_data, capsize=3)
        ax.set_title(f"Species: {species}")
        ax.set_xlabel("Sample Name")
        ax.set_ylabel(var_name)
        ax.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    return fig