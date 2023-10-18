from tools.mixing_rules import MixRules4VdW, MixRules4SRK
from tools.critic_values import extract_critical_properties_4_VdW, extract_critical_properties_SRK
import numpy as np


def calc_VdW_volume(substances:list[str], molar_fraction:list[float], Pressure:list[float], Temperature:list[float], R=83.14):
    """
    Calculate the molar voluume for Van der Waals EoS using the quadritic equation.

    Args:
        P (float): Pressure.
        T (float): Temprature.
        R (float, optional): Gas constant. Default value is 83.14 bar·cm^3/(mol·K).

    Returns:
        float: Positive molar volume in the quadratic equation.
    """
    list_substances = substances
    Tc_values, Pc_values = extract_critical_properties_4_VdW(list_substances)
    mix_rules = MixRules4VdW(critical_pressure=Pc_values, critical_temperature=Tc_values)
    a = mix_rules.calc_a_VdW(molar_fraction)
    b = mix_rules.calc_b_VdW(molar_fraction)

     # Coefficent to quadratic equation: a = P, b = -(Pb + RT), c = Pb
    a = Pressure
    b_coef = -(Pressure * b + R * Temperature)
    c = Pressure * b

    # Calculate the vectorized solutions using the quadratic equation
    discriminante = b_coef**2 - 4 * a * c
    solucion_1 = (-b_coef + np.sqrt(discriminante)) / (2 * a)
    solucion_2 = (-b_coef - np.sqrt(discriminante)) / (2 * a)

    # Find the positive value
    volumen = np.where(solucion_1 >= 0, solucion_1, solucion_2)
    return volumen

def calc_SRK_volume(substances:list[str], molar_fraction:list[float], Pressure:list[float], Temperature:list[float], R=83.14):
    """
    Calculate the molar voluume for Van der Waals EoS using the quadritic equation.

    Args:
        P (float): Pressure.
        T (float): Temprature.
        R (float, optional): Gas constant. Default value is 83.14 bar·cm^3/(mol·K).

    Returns:
        float: Positive molar volume in the quadratic equation.
    """
    list_substances = substances
    Tc_values, Pc_values, w_values = extract_critical_properties_SRK(list_substances)
    # Calc a and b parameter to Soave-Redlich-Kwong (SRK):
    m_values = 0.480 + 1.574 * w_values - 0.176 * w_values**2
    m = np.tile(m_values, (2, 1))
    mix_rules = MixRules4SRK(critical_temperature=Tc_values, critical_pressure=Pc_values)
    a_ij = mix_rules.calc_a_ij(m=m, temperature=Temperature)

    a = mix_rules.calc_a_SRK(molar_fractions=molar_fraction, a_ij=a_ij)
    b = mix_rules.calc_b_SRK(molar_fractions=molar_fraction)

     # Coefficent to quadratic equation: a = P, b = -(Pb + RT), c = Pb
    a = Pressure
    b_coef = -(R * Temperature)
    c = Pressure * b  + a

    # Calculate the vectorized solutions using the quadratic equation
    discriminante = b_coef**2 - 4 * a * c
    solucion_1 = (-b_coef + np.sqrt(discriminante)) / (2 * a)
    solucion_2 = (-b_coef - np.sqrt(discriminante)) / (2 * a)

    # Find the positive value
    volumen = np.where(solucion_1 >= 0, solucion_1, solucion_2)
    return volumen    