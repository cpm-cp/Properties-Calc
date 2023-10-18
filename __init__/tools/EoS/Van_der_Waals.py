from tools.mixing_rules import MixRules4VdW
from tools.critic_values import extract_critical_properties_4_VdW
import numpy as np


def calc_volume(substances:list[str], molar_fraction:list[float], Pressure:list[float], Temperature:list[float], R=83.14):
    """
    Calcula el valor de V en la ecuación de estado para un gas utilizando la ecuación cuadrática.

    Args:
        P (float): Presión.
        T (float): Temperatura.
        R (float, opcional): Constante de gas. Valor predeterminado es 83.14 bar·cm^3/(mol·K).

    Returns:
        float: La solución positiva para V en la ecuación cuadrática.
    """
    list_substances = substances
    Tc_values, Pc_values = extract_critical_properties_4_VdW(list_substances)
    mix_rules = MixRules4VdW(critical_pressure=Pc_values, critical_temperature=Tc_values)
    a = mix_rules.calc_a_VdW(molar_fraction)
    b = mix_rules.calc_b_VdW(molar_fraction)

     # Coeficientes de la ecuación cuadrática: a = P, b = -(Pb + RT), c = Pb
    a = Pressure
    b_coef = -(Pressure * b + R * Temperature)
    c = Pressure * b

    # Calcular las soluciones usando la fórmula cuadrática vectorizada
    discriminante = b_coef**2 - 4 * a * c
    solucion_1 = (-b_coef + np.sqrt(discriminante)) / (2 * a)
    solucion_2 = (-b_coef - np.sqrt(discriminante)) / (2 * a)

    # Encuentra las soluciones positivas para cada conjunto de temperaturas y presiones
    volumen = np.where(solucion_1 >= 0, solucion_1, solucion_2)

    return volumen