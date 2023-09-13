import pandas as pd

data = [
    ["Methane", 16.043, 0.012, 190.6, 45.99, 0.286, 98.6, 111.4],
    ["Ethane", 30.070, 0.100, 305.3, 48.72, 0.279, 145.5, 184.6],
    ["Propane", 44.097, 0.152, 369.8, 42.48, 0.276, 200.0, 231.1],
    ["n-Butane", 58.123, 0.200, 425.1, 37.96, 0.274, 255.0, 272.7],
    ["n-Pentane", 72.150, 0.252, 469.7, 33.70, 0.270, 313.0, 309.2],
    ["n-Hexane", 86.177, 0.301, 507.6, 30.25, 0.266, 371.0, 341.9],
    ["n-Heptane", 100.204, 0.350, 540.2, 27.40, 0.261, 428.0, 371.6],
    ["n-Octane", 114.231, 0.400, 568.7, 24.90, 0.256, 486.0, 398.8],
    ["n-Nonane", 128.258, 0.444, 594.6, 22.90, 0.252, 544.0, 424.0],
    ["n-Decane", 142.285, 0.492, 617.7, 21.10, 0.247, 600.0, 447.3],
    ["Isobutane", 58.123, 0.181, 408.1, 36.48, 0.282, 262.7, 261.4],
    ["Isooctane", 114.231, 0.302, 544.0, 25.68, 0.266, 468.0, 372.4],
    ["Cyclopentane", 70.134, 0.196, 511.8, 45.02, 0.273, 258.0, 322.4],
    ["Cyclohexane", 84.161, 0.210, 553.6, 40.73, 0.273, 308.0, 353.9],
    ["Methylcyclopentane", 84.161, 0.230, 532.8, 37.85, 0.272, 319.0, 345.0],
    ["Methylcyclohexane", 98.188, 0.235, 572.2, 34.71, 0.269, 368.0, 374.1],
    ["Ethylene", 28.054, 0.087, 282.3, 50.40, 0.281, 131.0, 169.4],
    ["Propylene", 42.081, 0.140, 365.6, 46.65, 0.289, 188.4, 225.5],
    ["1-Butene", 56.108, 0.191, 420.0, 40.43, 0.277, 239.3, 266.9],
    ["cis-2-Butene", 56.108, 0.205, 435.6, 42.43, 0.273, 233.8, 276.9],
    ["trans-2-Butene", 56.108, 0.218, 428.6, 41.00, 0.275, 237.7, 274.0],
    ["1-Hexene", 84.161, 0.280, 504.0, 31.40, 0.265, 354.0, 336.3],
    ["Isobutylene", 56.108, 0.194, 417.9, 40.00, 0.275, 238.9, 266.3],
    ["1,3-Butadiene", 54.092, 0.190, 425.2, 42.77, 0.267, 220.4, 268.7],
    ["Cyclohexene", 82.145, 0.212, 560.4, 43.50, 0.272, 291.0, 356.1],
    ["Acetylene", 26.038, 0.187, 308.3, 61.39, 0.271, 113.0, 189.4],
    ["Benzene", 78.114, 0.210, 562.2, 48.98, 0.271, 259.0, 353.2],
    ["Toluene", 92.141, 0.262, 591.8, 41.06, 0.264, 316.0, 383.8],
    ["Ethylbenzene", 106.167, 0.303, 617.2, 36.06, 0.263, 374.0, 409.4],
    ["Cumene", 120.194, 0.326, 631.1, 32.09, 0.261, 427.0, 425.6],
    ["o-Xylene", 106.167, 0.310, 630.3, 37.34, 0.263, 369.0, 417.6],
    ["m-Xylene", 106.167, 0.326, 617.1, 35.36, 0.259, 376.0, 412.3],
    ["p-Xylene", 106.167, 0.322, 616.2, 35.11, 0.260, 379.0, 411.5],
    ["Styrene", 104.152, 0.297, 636.0, 38.40, 0.256, 352.0, 418.3],
    ["Naphthalene", 128.174, 0.302, 748.4, 40.51, 0.269, 413.0, 491.2],
    ["Biphenyl", 154.211, 0.365, 789.3, 38.50, 0.295, 502.0, 528.2],
    ["Formaldehyde", 30.026, 0.282, 408.0, 65.90, 0.223, 115.0, 254.1],
    ["Acetaldehyde", 44.053, 0.291, 466.0, 55.50, 0.221, 154.0, 294.0],
    ["Methyl acetate", 74.079, 0.331, 506.6, 47.50, 0.257, 228.0, 330.1],
    ["Ethyl acetate", 88.106, 0.366, 523.3, 38.80, 0.255, 286.0, 350.2],
    ["Acetone", 58.080, 0.307, 508.2, 47.01, 0.233, 209.0, 329.4],
    ["Methyl ethyl ketone", 72.107, 0.323, 535.5, 41.50, 0.249, 267.0, 352.8],
    ["Diethyl ether", 74.123, 0.281, 466.7, 36.40, 0.263, 280.0, 307.6],
    ["Methyl t-butyl ether", 88.150, 0.266, 497.1, 34.30, 0.273, 329.0, 328.4],
    ["Methanol", 32.042, 0.564, 512.6, 80.97, 0.224, 118.0, 337.9],
    ["Ethanol", 46.069, 0.645, 513.9, 61.48, 0.240, 167.0, 351.4],
    ["1-Propanol", 60.096, 0.622, 536.8, 51.75, 0.254, 219.0, 370.4],
    ["1-Butanol", 74.123, 0.594, 563.1, 44.23, 0.260, 275.0, 390.8],
    ["1-Hexanol", 102.177, 0.579, 611.4, 35.10, 0.263, 381.0, 430.6],
    ["2-Propanol", 60.096, 0.668, 508.3, 47.62, 0.248, 220.0, 355.4],
    ["Phenol", 94.113, 0.444, 694.3, 61.30, 0.243, 229.0, 455.0],
    ["Ethylene glycol", 62.068, 0.487, 719.7, 77.00, 0.246, 191.0, 470.5],
    ["Acetic acid", 60.053, 0.467, 592.0, 57.86, 0.211, 179.7, 391.1],
    ["n-Butyric acid", 88.106, 0.681, 615.7, 40.64, 0.232, 291.7, 436.4],
    ["Benzoic acid", 122.123, 0.603, 751.0, 44.70, 0.246, 344.0, 522.4],
    ["Acetonitrile", 41.053, 0.338, 545.5, 48.30, 0.184, 173.0, 354.8],
    ["Methylamine", 31.057, 0.281, 430.1, 74.60, 0.321, 154.0, 266.8],
    ["Ethylamine", 45.084, 0.285, 456.2, 56.20, 0.307, 207.0, 289.7],
    ["Nitromethane", 61.040, 0.348, 588.2, 63.10, 0.223, 173.0, 374.4],
    ["Carbon tetrachloride", 153.822, 0.193, 556.4, 45.60, 0.272, 276.0, 349.8],
    ["Chloroform", 119.377, 0.222, 536.4, 54.72, 0.293, 239.0, 334.3],
    ["Dichloromethane", 84.932, 0.199, 510.0, 60.80, 0.265, 185.0, 312.9],
    ["Methyl chloride", 50.488, 0.153, 416.3, 66.80, 0.276, 143.0, 249.1],
    ["Ethyl chloride", 64.514, 0.190, 460.4, 52.70, 0.275, 200.0, 285.4],
    ["Chlorobenzene", 112.558, 0.250, 632.4, 45.20, 0.265, 308.0, 404.9],
    ["Tetrafluoroethane", 102.030, 0.327, 374.2, 40.60, 0.258, 198.0, 247.1],
    ["Argon", 39.948, 0.000, 150.9, 48.98, 0.291, 74.6, 87.3],
    ["Krypton", 83.800, 0.000, 209.4, 55.02, 0.288, 91.2, 119.8],
    ["Xenon", 131.30, 0.000, 289.7, 58.40, 0.286, 118.0, 165.0],
    ["Helium", 4.003, -0.390, 5.2, 2.28, 0.302, 57.3, 4.2],
    ["Hydrogen", 2.016, -0.216, 33.19, 13.13, 0.305, 64.1, 20.4],
    ["Oxygen", 31.999, 0.022, 154.6, 50.43, 0.288, 73.4, 90.2],
    ["Nitrogen", 28.014, 0.038, 126.2, 34.00, 0.289, 89.2, 77.3],
    ["Air", 28.851, 0.035, 132.2, 37.45, 0.289, 84.8, 0.00],
    ["Chlorine", 70.905, 0.069, 417.2, 77.10, 0.265, 124.0, 239.1],
    ["Carbon monoxide", 28.010, 0.048, 132.9, 34.99, 0.299, 93.4, 81.7],
    ["Carbon dioxide", 44.010, 0.224, 304.2, 73.83, 0.274, 94.0, 0.00],
    ["Carbon disulfide", 76.143, 0.111, 552.0, 79.00, 0.275, 160.0, 319.4],
    ["Hydrogen sulfide", 34.082, 0.094, 373.5, 89.63, 0.284, 98.5, 212.8],
    ["Sulfur dioxide", 64.065, 0.245, 430.8, 78.84, 0.269, 122.0, 263.1],
    ["Sulfur trioxide", 80.064, 0.424, 490.9, 82.10, 0.255, 127.0, 317.9],
    ["Nitric oxide (NO)", 30.006, 0.583, 180.2, 64.80, 0.251, 58.0, 121.4],
    ["Nitrous oxide (N2O)", 44.013, 0.141, 309.6, 72.45, 0.274, 97.4, 184.7],
    ["Hydrogen chloride", 36.461, 0.132, 324.7, 83.10, 0.249, 81.0, 188.2],
    ["Hydrogen cyanide", 27.026, 0.410, 456.7, 53.90, 0.197, 139.0, 298.9],
    ["Water", 18.015, 0.345, 647.1, 220.55, 0.229, 55.9, 373.2],
    ["Ammonia", 17.031, 0.253, 405.7, 112.80, 0.242, 72.5, 239.7],
    ["Nitric acid", 63.013, 0.714, 520.0, 68.90, 0.231, 145.0, 356.2],
    ["Sulfuric acid", 98.080, 0.00, 924.0, 64.00, 0.147, 177.0, 610.0]
]

columns = ["Molecule", "Molar mass", "ω", "Tc/K", "Pc/bar", "Zc", "cm3.mol-1", "Tn/K"]

df = pd.DataFrame(data, columns=columns)


# Assuming you have already created the DataFrame 'df' from your data

def calculate_Tc_mixing(substances, molar_fractions, df):
    """
    Calculate the critical temperature for a mixture of selected substances.
    
    Args:
        substances (list): List of substances ("Molecule") involved in the mixture.
        molar_fractions (list): List of molar fractions corresponding to the selected substances.
        df (pd.DataFrame): DataFrame containing substance properties, including 'Molecule' and 'Tc'.
        
    Returns:
        float: The calculated critical temperature for the mixture of selected substances.
    """
    if len(substances) != len(molar_fractions):
        raise ValueError("Number of substances and molar fractions must be the same.")
    
    Tc_mixing = 0.0
    
    for substance, molar_fraction in zip(substances, molar_fractions):
        Tc_substance = df[df['Molecule'] == substance]['Tc/K'].values[0]
        Tc_mixing += molar_fraction * Tc_substance
    
    return Tc_mixing


def calculate_Vc_mixing(substances, molar_fractions, df):
    """
    Calculate the critical volume for a mixture of selected substances.
    
    Args:
        substances (list): List of substances ("Molecule") involved in the mixture.
        molar_fractions (list): List of molar fractions corresponding to the selected substances.
        df (pd.DataFrame): DataFrame containing substance properties, including 'Molecule' and 'Vc'.
        
    Returns:
        float: The calculated critical volume for the mixture of selected substances.
    """
    if len(substances) != len(molar_fractions):
        raise ValueError("Number of substances and molar fractions must be the same.")
    
    Vc_mixing = 0.0
    
    for substance, molar_fraction in zip(substances, molar_fractions):
        Vc_substance = df[df['Molecule'] == substance]['cm3.mol-1'].values[0]
        Vc_mixing += molar_fraction * Vc_substance
    
    return Vc_mixing


def calculate_Zc_mixing(substances, molar_fractions, df):
    """
    Calculate the critical compressibility factor for a mixture of selected substances.
    
    Args:
        substances (list): List of substances ("Molecule") involved in the mixture.
        molar_fractions (list): List of molar fractions corresponding to the selected substances.
        df (pd.DataFrame): DataFrame containing substance properties, including 'Molecule' and 'Zc'.
        
    Returns:
        float: The calculated critical volume for the mixture of selected substances.
    """
    if len(substances) != len(molar_fractions):
        raise ValueError("Number of substances and molar fractions must be the same.")
    
    Zc_mixing = 0.0
    
    for substance, molar_fraction in zip(substances, molar_fractions):
        Zc_substance = df[df['Molecule'] == substance]['Zc'].values[0]
        Zc_mixing += molar_fraction * Zc_substance
    
    return Zc_mixing

def calculate_Omega_mixing(substances, molar_fractions, df):
    """
    Calculate the Omega for a mixture of selected substances.
    
    Args:
        substances (list): List of substances ("Molecule") involved in the mixture.
        molar_fractions (list): List of molar fractions corresponding to the selected substances.
        df (pd.DataFrame): DataFrame containing substance properties, including 'Molecule' and 'ω'.
        
    Returns:
        float: The calculated omega value for the mixture of selected substances.
    """
    if len(substances) != len(molar_fractions):
        raise ValueError("Number of substances and molar fractions must be the same.")
    
    Omega_mixing = 0.0
    
    for substance, molar_fraction in zip(substances, molar_fractions):
        Omega_substance = df[df['Molecule'] == substance]['ω'].values[0]
        Omega_mixing += molar_fraction * Omega_substance
    
    return Omega_mixing

def calculate_Pc_mixing(Zc_mixture: float, Tc_mixture: float, Vc_mixture: float, R: float = 83.14):
    '''
    Calculate the critical pressure for a mixture of selected substances.

    Args:
        Tc_mixture (float): Critical temperature of the mixture (in K).
        Zc_mixture (float): Compressibility factor of the mixture.
        R (float): Gas constant (in cm³ bar / mol K, default is 8.314).
        Vc_mixture (float): Critical volume of the mixture.
        
    Returns:
        float: The calculated critical pressure for the mixture (in bar).
    '''
    Pc_mixing = (Zc_mixture * R * Tc_mixture) / Vc_mixture

    return Pc_mixing

def calculate_mixture_properties(selected_substances, selected_molar_fractions, df=df):
    '''
    Calculate mixture properties for a given list of substances and their molar fractions.

    Args:
        selected_substances (list of str): List of selected substances.
        selected_molar_fractions (list of float): List of molar fractions corresponding to the selected substances.
        df (pd.DataFrame): DataFrame containing the properties of all substances.

    Returns:
        dict: A dictionary containing the calculated mixture properties.
    '''
    # Calculate Tc_mixture
    Tc_mixture = calculate_Tc_mixing(selected_substances, selected_molar_fractions, df)

    # Calculate Vc_mixture
    Vc_mixture = calculate_Vc_mixing(selected_substances, selected_molar_fractions, df)

    # Calculate Zc_mixture
    Zc_mixture = calculate_Zc_mixing(selected_substances, selected_molar_fractions, df)

    # Calculate Omega_mixture
    Omega_mixture = calculate_Omega_mixing(selected_substances, selected_molar_fractions, df)

    # Calculate Pc_mixture using previously calculated values
    Pc_mixture = calculate_Pc_mixing(Zc_mixture, Tc_mixture, Vc_mixture)

    # Print the calculated properties
    print(80*"-")
    print(f"Critical Temperature for the Mixture: {Tc_mixture:.4f} K")
    print(f"Critical Volume for the Mixture: {Vc_mixture:.4f} cm3.mol-1")
    print(f"Critical Z for the Mixture: {Zc_mixture:.4f}")
    print(f"Critical ω for the Mixture: {Omega_mixture}")
    print(f"Critical Pc for the Mixture: {Pc_mixture:.4f} bar")
    print(80*"-")