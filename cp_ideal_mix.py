from pandas import DataFrame
from scipy.integrate import quad

data = [
    ["Methane", 1500, 4.217, 1.702, 9.081, -2.164, 0.0],
    ["Ethane", 1500, 6.369, 1.131, 19.225, -5.561, 0.0],
    ["Propane", 1500, 9.011, 1.213, 28.785, -8.824, 0.0],
    ["n-Butane", 1500, 11.928, 1.935, 36.915, -11.402, 0.0],
    ["iso-Butane", 1500, 11.901, 1.677, 37.853, -11.945, 0.0],
    ["n-Pentane", 1500, 14.731, 2.464, 45.351, -14.111, 0.0],
    ["n-Hexane", 1500, 17.550, 3.025, 53.722, -16.791, 0.0],
    ["n-Heptane", 1500, 20.361, 3.570, 62.127, -19.486, 0.0],
    ["n-Octane", 1500, 23.174, 4.108, 70.567, -22.208, 0.0],
    ["Ethylene", 1500, 5.325, 1.424, 14.394, -4.392, 0.0],
    ["Propylene", 1500, 7.792, 1.637, 22.706, -6.915, 0.0],
    ["1-Butene", 1500, 10.520, 1.967, 31.630, -9.873, 0.0],
    ["1-Pentene", 1500, 13.437, 2.691, 39.753, -12.447, 0.0],
    ["1-Hexene", 1500, 16.240, 3.220, 48.189, -15.157, 0.0],
    ["1-Heptene", 1500, 19.053, 3.768, 56.588, -17.847, 0.0],
    ["1-Octene", 1500, 21.868, 4.324, 64.960, -20.521, 0.0],
    ["Acetaldehyde", 1000, 6.506, 1.693, 17.978, -6.158, 0.0],
    ["Acetylene", 1500, 5.253, 6.132, 1.952, 0.0, -1.299],
    ["Benzene", 1500, 10.259, -0.206, 39.064, -13.301, 0.0],
    ["1,3-Butadiene", 1500, 10.720, 2.734, 26.786, -8.882, 0.0],
    ["Cyclohexane", 1500, 13.121, -3.876, 63.249, -20.928, 0.0],
    ["Ethanol", 1500, 8.948, 3.518, 20.001, -6.002, 0.0],
    ["Ethylbenzene", 1500, 15.993, 1.124, 55.380, -18.476, 0.0],
    ["Ethylene oxide", 1000, 5.784, -0.385, 23.463, -9.296, 0.0],
    ["Formaldehyde", 1500, 4.191, 2.264, 7.022, -1.877, 0.0],
    ["Methanol", 1500, 5.547, 2.211, 12.216, -3.450, 0.0],
    ["Styrene", 1500, 15.534, 2.050, 50.192, -16.662, 0.0],
    ["Toluene", 1500, 12.922, 0.290, 47.052, -15.716, 0.0],
    ["Air", 2000, 3.509, 3.355, 0.575, 0.0, -0.016],
    ["Ammonia", 1800, 4.269, 3.578, 3.020, 0.0, -0.186],
    ["Bromine", 3000, 4.337, 4.493, 0.056, 0.0, -0.154],
    ["Carbon monoxide", 2500, 3.507, 3.376, 0.557, 0.0, -0.031],
    ["Carbon dioxide", 2000, 4.467, 5.457, 1.045, 0.0, -1.157],
    ["Carbon disulfide", 1800, 5.532, 6.311, 0.805, 0.0, -0.906],
    ["Chlorine", 3000, 4.082, 4.442, 0.089, 0.0, -0.344],
    ["Hydrogen", 3000, 3.468, 3.249, 0.422, 0.0, 0.083],
    ["Hydrogen sulfide", 2300, 4.114, 3.931, 1.490, 0.0, -0.232],
    ["Hydrogen chloride", 2000, 3.512, 3.156, 0.623, 0.0, 0.151],
    ["Hydrogen cyanide", 2500, 4.326, 4.736, 1.359, 0.0, -0.725],
    ["Nitrogen", 2000, 3.502, 3.280, 0.593, 0.0, 0.040],
    ["Nitrous oxide", 2000, 4.646, 5.328, 1.214, 0.0, -0.928],
    ["Nitric oxide", 2000, 3.590, 3.387, 0.629, 0.0, 0.014],
    ["Nitrogen dioxide", 2000, 4.447, 4.982, 1.195, 0.0, -0.792],
    ["Dinitrogen tetroxide", 2000, 9.198, 11.660, 2.257, 0.0, -2.787],
    ["Oxygen", 2000, 3.535, 3.639, 0.506, 0.0, -0.227],
    ["Sulfur dioxide", 2000, 4.796, 5.699, 0.801, 0.0, -1.015],
    ["Sulfur trioxide", 2000, 6.094, 8.060, 1.056, 0.0, -2.028],
    ["Water", 2000, 4.038, 3.470, 1.450, 0.0, 0.121]
]

columns = ["Substance", "Tmax", "CP_298/R", "A", "B", "C", "D"]

df = DataFrame(data, columns=columns)

# Convert values to B, C & D
df['B'] /= 10**3
df['C'] /= 10**6
df['D'] /= 10**-5

def calculate_H_ideal_mix(substances: list, molar_fractions: float, Tmin:float, Tmax:float, R: float = 8.314, df: DataFrame = df, current: str = "entry") -> float:
    """
    Calculate the calorific capacity for a mixture of selected substances using four individual function (V.0.1)[ignore it].
    Using the following equation: Cp_ideal_mix / R = A + B*T + C*T^2 + D*T^-2 
    {Smith; Van Ness}  
    
    Args:
        substances (list): List of substances ("Substance") involved in the mixture.
        molar_fractions (list): List of molar fractions corresponding to the selected substances.
        df (pd.DataFrame): DataFrame containing substance properties, including 'Substance' and 'A', 'B', 'C' & 'D' constants.
        
    Returns:
        float: The calculated calorific capacity at constant pressure for the mixture of selected substances.
    """
    def A_mix(df) -> float:
        """Calculate the A constant for the mixing as ideal gas

        Args:
            df (df): Admit the A constant values for a N substance in a mixing

        Returns:
            float: Return de A constant mixing value as ideal gas
        """
        A_mixing = 0.00
        for substance, molar_fraction in zip(substances, molar_fractions):
            A_substance = df[df['Substance'] == substance]['A'].values[0]
            A_mixing += molar_fraction * A_substance
        return A_mixing

    def B_mix(df) -> float:
        """Calculate the B constant for the mixing as ideal gas

        Args:
            df (df): Admit the B constant values for a N substance in a mixing

        Returns:
            float: Return de B constant mixing value as ideal gas
        """
        B_mixing = 0.00
        for substance, molar_fraction in zip(substances, molar_fractions):
            B_substance = df[df['Substance'] == substance]['B'].values[0]
            B_mixing += molar_fraction * B_substance
        return B_mixing

    def C_mix(df) -> float:
        """Calculate the C constant for the mixing as ideal gas

        Args:
            df (df): Admit the C constant values for a N substance in a mixing

        Returns:
            float: Return de C constant mixing value as ideal gas
        """
        C_mixing = 0.00
        for substance, molar_fraction in zip(substances, molar_fractions):
            C_substance = df[df['Substance'] == substance]['C'].values[0]
            C_mixing += molar_fraction * C_substance
        return C_mixing
    
    def D_mix(df) -> float:
        """Calculate the D constant for the mixing as ideal gas

        Args:
            df (df): Admit the D constant values for a N substance in a mixing

        Returns:
            float: Return de D constant mixing value as ideal gas
        """
        D_mixing = 0.00
        for substance, molar_fraction in zip(substances, molar_fractions):
            D_substance = df[df['Substance'] == substance]['D'].values[0]
            D_mixing += molar_fraction * D_substance
        return D_mixing
    
    
    A = A_mix(df)
    B = B_mix(df)
    C = C_mix(df)
    D = D_mix(df)
    
    # Check if Tmin is zero before performing division
    if Tmin == 0:
        enthalpy_value = (A * (Tmax - Tmin) + B / 2 * (Tmax**2 - Tmin**2) + C / 3 * (Tmax**3 - Tmin**3) - D * Tmax**-1) * R
    else:
        enthalpy_value = (A * (Tmax - Tmin) + B / 2 * (Tmax**2 - Tmin**2) + C / 3 * (Tmax**3 - Tmin**3) - D * (Tmax**-1 - Tmin**-1)) * R

    print(f"The ideal enthalpy mix value at {current} is: {enthalpy_value:.4f} kJ/kmol")
    print(80 * "-")
    return enthalpy_value

    