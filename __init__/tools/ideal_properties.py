from pandas import DataFrame
import numpy as np

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

def calculate_ideal_mix_properties(substances: list, molar_fractions: float, T_reference: float, T_state: float, P_reference:float, P_state: float, R: float= 8.314, df: DataFrame=df, current: str="entry") -> tuple[float, float]:
    """This function calculate the enthalpy and entropy ideal by the mix for a specific current.
    At moment to obtain values should be in this order:
    - Ideal enthalpy.
    - Ideal entropy. 

    Args:
        substances (list): Substance names in English.
        molar_fractions (float): Molar fractions by the components.
        T_reference (float): Reference temperature or min temperature.
        T_state (float): Temperature by the state (current).
        P_reference(float): Reference pressure or min pressure.
        P_state(float): Pressure by the state (current).
        R (float, optional): Gas constant in kJ/kmol*K. Defaults to 8.314.
        df (DataFrame, optional): Data information. Defaults to df.
        current (str, optional): Current status. Defaults to "entry".

    Returns:
        float: Enthalpy and entropy ideal mix.
    """
    A_values = df[df['Substance'].isin(substances)]['A'].values
    B_values = df[df['Substance'].isin(substances)]['B'].values
    C_values = df[df['Substance'].isin(substances)]['C'].values
    D_values = df[df['Substance'].isin(substances)]['D'].values

    A_mixing = np.dot(molar_fractions, A_values)
    B_mixing = np.dot(molar_fractions, B_values)
    C_mixing = np.dot(molar_fractions, C_values)
    D_mixing = np.dot(molar_fractions, D_values)

    if T_reference == 0:
        enthalpy_value = (A_mixing * (T_state - T_reference) + B_mixing / 2 * (T_state**2 - T_reference**2) + C_mixing / 3 * (T_state**3 - T_reference**3) - D_mixing * T_state**-1) * R
    else:
        enthalpy_value = (A_mixing * (T_state - T_reference) + B_mixing / 2 * (T_state**2 - T_reference**2) + C_mixing / 3 * (T_state**3 - T_reference**3) - D_mixing * (T_state**-1 - T_reference**-1)) * R

    if T_reference == 0:
        entropy_value = (A_mixing * np.log(T_state / T_reference) + B_mixing * (T_state - T_reference) + C_mixing / 2 * (T_state**2 - T_reference**2) - (D_mixing / 2) * T_state**-1) * R - np.log(P_state / P_reference)
    else:
        entropy_value = (A_mixing * np.log(T_state / T_reference) + B_mixing * (T_state - T_reference) + C_mixing / 2 * (T_state**2 - T_reference**2) - (D_mixing / 2) * (T_state**-2 - T_reference**-2)) * R - np.log(P_state / P_reference)

    return enthalpy_value, entropy_value

    