from tools.mixing_rules import MixRules4VdW
from tools.critic_values import extract_critical_properties_4_VdW
import numpy as np


def calc_volume(substances:list[str], molar_fraction:list[float], Pressure:list[float], Temperature:list[float], R=83.14):
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