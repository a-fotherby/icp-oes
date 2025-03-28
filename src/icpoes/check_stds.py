"""All check standard values given in mM"""

"""
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
"""

# Define atomic weights (g/mol) for the elements involved
atomic_weights = {
    # SLRS-6 and common elements:
    "Al": 26.98,
    "Sb": 121.76,
    "As": 74.92,
    "Ba": 137.33,
    "Be": 9.0122,
    "Cd": 112.41,
    "Ca": 40.08,
    "Cr": 52.00,
    "Co": 58.93,
    "Cu": 63.55,
    "Fe": 55.85,
    "Pb": 207.2,
    "Mg": 24.31,
    "Mn": 54.94,
    "Mo": 95.95,
    "Ni": 58.69,
    "K": 39.10,
    "Na": 22.99,
    "Sr": 87.62,
    "U": 238.03,
    "V": 50.94,
    "Zn": 65.38,
    # Additional elements for SPS-SW2 10%
    "B": 10.81,
    "Ce": 140.12,
    "Cs": 132.91,
    "Dy": 162.50,
    "Er": 167.26,
    "Eu": 151.96,
    "Gd": 157.25,
    "Ho": 164.93,
    "La": 138.91,
    "Lu": 174.97,
    "Nd": 144.24,
    "P": 30.97,
    "Pr": 140.91,
    "Rb": 85.47,
    "S": 32.06,
    "Sc": 44.96,
    "Se": 78.96,
    "Si": 28.09,
    "Sm": 150.36,
    "Tb": 158.93,
    "Th": 232.04,
    "Tl": 204.38,
    "Tm": 168.93,
    "Y": 88.91,
    "Yb": 173.04
}


def slrs6():
# SLRS-6 (converted from mM to ppm)
    slrs6 = {
        "Al": 0.03373,
        "Sb": 0.0003374,
        "As": 0.0005701,
        "Ba": 0.01428,
        "Be": 6.60e-06,
        "Cd": 6.30e-06,
        "Ca": 8.74,
        "Cr": 0.0002522,
        "Co": 5.30e-06,
        "Cu": 0.02390,
        "Fe": 0.08432,
        "Pb": 0.0001701,
        "Mg": 2.133,
        "Mn": 0.00212,
        "Mo": 0.000215,
        "Ni": 0.0006157,
        "K": 0.653,
        "Na": 2.7588,
        "Sr": 0.04062,
        "U": 6.97e-05,
        "V": 0.000350,
        "Zn": 0.001759
    }

    slrs6_uncertainties_ppm = {
        "Al": 0.0022136,
        "Sb": 5.84e-06,
        "As": 8.02e-05,
        "Ba": 0.00048066,
        "Be": 2.20e-06,
        "Cd": 1.41e-06,
        "Ca": 0.2004,
        "Cr": 1.20e-05,
        "Co": 1.20e-05,
        "Cu": 1.78e-04,
        "Fe": 0.00357,
        "Pb": 2.61e-05,
        "Mg": 0.05834,
        "Mn": 9.89e-05,
        "Mo": 1.82e-05,
        "Ni": 2.23e-05,
        "K": 0.05474,
        "Na": 0.2207,
        "Sr": 0.00032409,
        "U": 3.33e-07,
        "V": 6.11e-06,
        "Zn": 1.18e-04
    }

    return slrs6, slrs6_uncertainties_ppm

def sls_sw2_10per():
# SPS-SW2 10%: Conversion from mM to ppm (without dividing by 10)
    values = {
        'Al': 0.025,
        'As': 0.005,
        'Bb': 0.025,
        'Ba': 0.025,
        'Ca': 100.0,
        'Cd': 0.00025,
        'Ce': 0.00025,
        'Co': 0.001,
        'Cr': 0.001,
        'Cs': 0.001,
        'Cu': 0.01,
        'Dy': 0.00025,
        'Er': 0.00025,
        'Eu': 0.00025,
        'Fe': 0.01,
        'Gd': 0.00025,
        'Ho': 0.00025,
        'K':  0.1,
        'La': 0.00025,
        'Lu': 0.00025,
        'Mg': 0.2,
        'Mn': 0.005,
        'Mo': 0.005,
        'Na': 1.0,
        'Nd': 0.00025,
        'Ni': 0.005,
        'P': 0.05,
        'Pb': 0.0025,
        'Pr': 0.00025,
        'Rb': 0.005,
        'S': 1.0,
        'Sc': 0.00025,
        'Se': 0.001,
        'Si': 0.5,
        'Sm': 0.00025,
        'Sr': 0.025,
        'Tb': 0.00025,
        'Th': 0.00025,
        'Tl': 0.00025,
        'Tm': 0.00025,
        'U': 0.00025,
        'V': 0.005,
        'Y': 0.00025,
        'Yb': 0.00025,
        'Zn': 0.01
    }

    errors = {
        'Al': 0.0001,
        'As': 0.00003,
        'Bb': None,
        'Ba': 0.0001,
        'Ca': 0.005,
        'Cd': 2e-06,
        'Ce': 2e-06,
        'Co': 5e-06,
        'Cr': 5e-06,
        'Cs': 5e-06,
        'Cu': 0.0001,
        'Dy': 2e-06,
        'Er': 2e-06,
        'Eu': 2e-06,
        'Fe': 0.0001,
        'Gd': 2e-06,
        'Ho': 2e-06,
        'K': 0.0005,
        'La': 2e-06,
        'Lu': 2e-06,
        'Mg': 0.001,
        'Mn': 0.00003,
        'Mo': 0.00003,
        'Na': 0.005,
        'Nd': 2e-06,
        'Ni': 0.00003,
        'P': 0.0003,
        'Pb': 0.00001,
        'Pr': 2e-06,
        'Rb': 0.00003,
        'S': 0.005,
        'Sc': 2e-06,
        'Se': 5e-06,
        'Si': 0.003,
        'Sm': 2e-06,
        'Sr': 0.0001,
        'Tb': 2e-06,
        'Th': 2e-06,
        'Tl': 2e-06,
        'Tm': 2e-06,
        'U': 2e-06,
        'V': 0.00003,
        'Y': 2e-06,
        'Yb': 2e-06,
        'Zn': 0.0002
    }

    return values, errors

def nist1640a():
    values = {
        'Al': 0.25,
        'As': 0.05,
        'Bb': 0.25,
        'Ba': 0.25,
        'Ca': 10.0,
        'Cd': 0.0025,
        'Ce': 0.0025,
        'Co': 0.01,
        'Cr': 0.01,
        'Cs': 0.01,
        'Cu': 0.1,
        'Dy': 0.0025,
        'Er': 0.0025,
        'Eu': 0.0025,
        'Fe': 0.1,
        'Gd': 0.0025,
        'Ho': 0.0025,
        'K': 1.0,
        'La': 0.0025,
        'Lu': 0.0025,
        'Mg': 2.0,
        'Mn': 0.05,
        'Mo': 0.05,
        'Na': 10.0,
        'Nd': 0.0025,
        'Ni': 0.05,
        'P': 0.5,
        'Pb': 0.025,
        'Pr': 0.0025,
        'Rb': 0.05,
        'S': 10.0,
        'Sc': 0.0025,
        'Se': 0.01,
        'Si': 5.0,
        'Sm': 0.0025,
        'Sr': 0.25,
        'Tb': 0.0025,
        'Th': 0.0025,
        'Tl': 0.0025,
        'Tm': 0.0025,
        'U': 0.0025,
        'V': 0.05,
        'Y': 0.0025,
        'Yb': 0.0025,
        'Zn': 0.1
    }

    errors = {
        'Al': 0.001,
        'As': 0.0003,
        'Bb': None,
        'Ba': 0.001,
        'Ca': 0.05,
        'Cd': 2e-05,
        'Ce': 2e-05,
        'Co': 5e-05,
        'Cr': 5e-05,
        'Cs': 5e-05,
        'Cu': 0.001,
        'Dy': 2e-05,
        'Er': 2e-05,
        'Eu': 2e-05,
        'Fe': 0.001,
        'Gd': 2e-05,
        'Ho': 2e-05,
        'K': 0.005,
        'La': 2e-05,
        'Lu': 2e-05,
        'Mg': 0.01,
        'Mn': 0.0003,
        'Mo': 0.0003,
        'Na': 0.05,
        'Nd': 2e-05,
        'Ni': 0.0003,
        'P': 0.003,
        'Pb': 0.0001,
        'Pr': 2e-05,
        'Rb': 0.0003,
        'S': 0.05,
        'Sc': 2e-05,
        'Se': 5e-05,
        'Si': 0.03,
        'Sm': 2e-05,
        'Sr': 0.001,
        'Tb': 2e-05,
        'Th': 2e-05,
        'Tl': 2e-05,
        'Tm': 2e-05,
        'U': 2e-05,
        'V': 0.0003,
        'Y': 2e-05,
        'Yb': 2e-05,
        'Zn': 0.002
    }
    return values, errors


def check_standards():
    # Concentration values in mM
    # Keys must match name in dataset
    standards = {'SLRS-6': slrs6(),
                 'SPS-SW2 10%': sls_sw2_10per(),
                 'NIST 1640a': nist1640a()}
    return standards