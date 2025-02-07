"""All check standard values given in mM"""

def slrs6():
    # Concentration values in mM
    slrs6_concentrations = {
        "Al": 1.25e-3,
        "Sb": 2.77e-6,
        "As": 7.61e-6,
        "Ba": 1.04e-4,
        "Be": 7.32e-7,
        "Cd": 5.61e-8,
        "Ca": 0.218,
        "Cr": 4.85e-6,
        "Co": 9.00e-7,
        "Cu": 3.76e-4,
        "Fe": 1.51e-3,
        "Pb": 8.21e-7,
        "Mg": 0.0877,
        "Mn": 3.86e-5,
        "Mo": 2.24e-6,
        "Ni": 1.05e-5,
        "K": 1.67e-2,
        "Na": 0.120,
        "Sr": 4.64e-4,
        "U": 2.93e-7,
        "V": 6.89e-6,
        "Zn": 2.69e-5
    }

    # Uncertainty values in mM
    slrs6_uncertainties = {
        "Al": 8.2e-5,
        "Sb": 4.8e-8,
        "As": 1.07e-6,
        "Ba": 3.5e-6,
        "Be": 2.44e-7,
        "Cd": 1.25e-8,
        "Ca": 0.0050,
        "Cr": 2.3e-7,
        "Co": 2.04e-7,
        "Cu": 2.8e-6,
        "Fe": 6.4e-5,
        "Pb": 1.26e-7,
        "Mg": 0.0024,
        "Mn": 1.8e-6,
        "Mo": 1.9e-7,
        "Ni": 3.8e-7,
        "K": 1.4e-3,
        "Na": 0.0096,
        "Sr": 3.7e-6,
        "U": 1.4e-9,
        "V": 1.2e-7,
        "Zn": 1.8e-6
    }

    return slrs6_concentrations, slrs6_uncertainties

def sls_sw2_10per():
    # Original dictionaries
    mm_values = {
        "Al": 0.00926,
        "As": 0.000667,
        "B": 0.0231,
        "Ba": 0.00182,
        "Ca": 0.250,
        "Cd": 2.23e-05,
        "Ce": 1.78e-05,
        "Co": 0.000170,
        "Cr": 0.000192,
        "Cs": 7.52e-05,
        "Cu": 0.00157,
        "Dy": 1.54e-05,
        "Er": 1.49e-05,
        "Eu": 1.65e-05,
        "Fe": 0.00179,
        "Gd": 1.59e-05,
        "Ho": 1.52e-05,
        "K": 0.0256,
        "La": 1.80e-05,
        "Lu": 1.43e-05,
        "Mg": 0.0823,
        "Mn": 0.000910,
        "Mo": 0.000521,
        "Na": 0.435,
        "Nd": 1.73e-05,
        "Ni": 0.000853,
        "P": 0.0161,
        "Pb": 0.000121,
        "Pr": 1.78e-05,
        "Rb": 0.000586,
        "S": 0.312,
        "Sc": 5.56e-05,
        "Se": 0.000127,
        "Si": 0.1781,
        "Sm": 1.66e-05,
        "Sr": 0.00285,
        "Tb": 1.57e-05,
        "Th": 1.08e-05,
        "Tl": 1.22e-05,
        "Tm": 1.48e-05,
        "U": 1.05e-05,
        "V": 0.000982,
        "Y": 2.81e-05,
        "Yb": 1.44e-05,
        "Zn": 0.00153
    }

    mm_uncertainties = {
        "Al": 3.71e-05,
        "As": 4.00e-06,
        "B": 0.0,
        "Ba": 7.28e-06,
        "Ca": 0.00125,
        "Cd": 1.78e-07,
        "Ce": 1.43e-07,
        "Co": 8.49e-07,
        "Cr": 9.62e-07,
        "Cs": 3.76e-07,
        "Cu": 1.57e-05,
        "Dy": 1.23e-07,
        "Er": 1.20e-07,
        "Eu": 1.32e-07,
        "Fe": 1.79e-05,
        "Gd": 1.27e-07,
        "Ho": 1.21e-07,
        "K": 1.28e-04,
        "La": 1.44e-07,
        "Lu": 1.14e-07,
        "Mg": 4.11e-04,
        "Mn": 5.46e-06,
        "Mo": 3.13e-06,
        "Na": 0.00217,
        "Nd": 1.39e-07,
        "Ni": 5.11e-06,
        "P": 9.68e-05,
        "Pb": 4.83e-07,
        "Pr": 1.42e-07,
        "Rb": 3.51e-06,
        "S": 0.00156,
        "Sc": 4.45e-07,
        "Se": 6.33e-07,
        "Si": 0.00107,
        "Sm": 1.33e-07,
        "Sr": 1.14e-05,
        "Tb": 1.26e-07,
        "Th": 8.62e-08,
        "Tl": 9.79e-08,
        "Tm": 1.18e-07,
        "U": 8.40e-08,
        "V": 5.89e-06,
        "Y": 2.25e-07,
        "Yb": 1.16e-07,
        "Zn": 3.06e-05
    }

    # Create new dictionaries with values scaled by 1/10
    mm_values_div10 = {element: val / 10 for element, val in mm_values.items()}
    mm_uncertainties_div10 = {element: val / 10 for element, val in mm_uncertainties.items()}

    return mm_values_div10, mm_uncertainties_div10

def check_standards():
    # Concentration values in mM
    # Keys must match name in dataset
    standards = {'SLRS-6': slrs6(),
                 'SPS-SW2 10%': sls_sw2_10per()}
    return standards
