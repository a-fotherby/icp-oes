import numpy as np
import matplotlib.pyplot as plt
import itertools

def bar_by_sample(results_ds, comp_samples=["SLRS-6_1", "SPS-SW2 10%_1"]):
    """
    Create a grouped bar chart figure comparing species abundances for each sample in results_ds.
    For each sample (from the full dataset) the plot shows:
      - The concentration for the current sample (as the left-most bar)
      - Comparison bars for the specified comp_samples (if present in the dataset)
    Error bars for each species are taken from the 'error' variable.
    
    Parameters
    ----------
    results_ds : xarray.Dataset
        The dataset containing your analysis results. It must have:
          - A coordinate 'sample_name' in the 'Concentration_ppm' variable.
          - A coordinate 'species' in the 'Concentration_ppm' variable.
          - Variables 'Concentration_ppm' (with dimensions (sample_name, species))
            and 'error' (with the same dimensions).
    comp_samples : list of str, optional
        A list of sample names to use for comparison (default is ["SLRS-6_1", "SPS-SW2 10%_1"]).
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the subplots.
    """
    conc = results_ds['Concentration_ppm']
    err = results_ds['error']
    
    # Get all sample names
    sample_names = conc.coords["sample_name"].values
    n_samples = len(sample_names)
    
    # Extract error values (assumed to be a 2D array with shape (n_samples, n_species))
    error_vals = err.values
    
    # Create a grid for subplots (one per sample)
    ncols = int(np.ceil(np.sqrt(n_samples)))
    nrows = int(np.ceil(n_samples / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
    axes = np.atleast_1d(axes).flatten()
    
    # x positions for species bars
    species = conc.coords["element"].values
    x = np.arange(len(species))
    bar_width = 0.25
    
    total_bars = 1 + len(comp_samples)
    offsets = [ -((total_bars - 1)/2)*bar_width + i*bar_width for i in range(total_bars) ]
    
    current_color = 'skyblue'
    comp_color_cycle = itertools.cycle(['salmon', 'limegreen', 'orange', 'violet'])
    
    for idx, sample in enumerate(sample_names):
        ax = axes[idx]
        # Current sample data
        curr_data = conc.sel(sample_name=sample).values
        ax.bar(x + offsets[0], curr_data, width=bar_width, label=sample,
               color=current_color, yerr=error_vals[idx], capsize=3)
        
        # Plot comparison samples if present
        for i, comp_sample in enumerate(comp_samples):
            if comp_sample in sample_names:
                comp_idx = np.where(sample_names == comp_sample)[0][0]
                comp_data = conc.sel(sample_name=comp_sample).values
                color = next(comp_color_cycle)
                ax.bar(x + offsets[i+1], comp_data, width=bar_width, label=comp_sample,
                       color=color, yerr=error_vals[comp_idx], capsize=3)
        comp_color_cycle = itertools.cycle(['salmon', 'limegreen', 'orange', 'violet'])
        ax.set_xlabel("Species")
        ax.set_ylabel("Concentration (mM)")
        ax.set_title(f"Sample: {sample}")
        ax.set_xticks(x)
        ax.set_xticklabels(species, rotation=45)
        ax.legend(fontsize='small')
    
    # Remove extra axes if any
    for j in range(idx+1, len(axes)):
        fig.delaxes(axes[j])
    
    fig.tight_layout()
    return fig

def bar_by_species(results_ds, var_name="Concentration_ppm", comp_samples=None):
    """
    Plot bar charts for each species from a results dataset and return the figure.
    For each species (one subplot per species), the function plots a bar for each sample,
    using error bars from the 'error' variable. Optionally, samples in comp_samples are highlighted.
    
    Parameters
    ----------
    results_ds : xarray.Dataset
        The dataset with coordinates 'sample_name' and 'species', and with variables
        'Concentration_ppm' and 'error'.
    var_name : str, optional
        The variable to plot (default is 'Concentration_ppm').
    comp_samples : list of str, optional
        A list of sample names to highlight (bars drawn in a different color).
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the subplots (one per species).
    """
    if comp_samples is None:
        comp_samples = []
    
    species_list = results_ds.coords["element"].values
    n_species = len(species_list)
    
    fig, axes = plt.subplots(n_species, 1, figsize=(10, 4 * n_species), squeeze=False)
    
    for i, sp in enumerate(species_list):
        ax = axes[i, 0]
        data = results_ds[var_name].sel(element=sp)
        error_data = results_ds["error"].sel(element=sp).values
        sample_names = data.coords["sample_name"].values
        values = data.values
        
        # Set bar colors: highlight if in comp_samples
        colors = ['red' if sample in comp_samples else 'blue' for sample in sample_names]
        ax.bar(sample_names, values, color=colors, yerr=error_data, capsize=3)
        ax.set_title(f"Species: {sp}")
        ax.set_xlabel("Sample Name")
        ax.set_ylabel(var_name)
        ax.tick_params(axis="x", rotation=45)
    
    plt.tight_layout()
    return fig
