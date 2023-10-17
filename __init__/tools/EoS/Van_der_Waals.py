from mixing_rules import MixingRules
from critic_values import extract_critical_properties_4_VdW
import numpy as np


def calc_volume(substances:list[str], Pressure, Temperature, R=83.14):
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
    mix_rules = MixingRules(critical_pressure=Pc_values, critical_temperature=Tc_values)
    a = mix_rules.calc_a()
    b = mix_rules.calc_b()

    # Coeficientes de la ecuación cuadrática
    coef_a = 1
    coef_b = -(R * Temperature / Pressure + b)
    coef_c = -a / Pressure
    
    # Calcular las soluciones usando la fórmula cuadrática
    discriminante = coef_b**2 - 4 * coef_a * coef_c
    
    if discriminante >= 0:
        solucion_1 = (-coef_b + np.sqrt(discriminante)) / (2 * coef_a)
        solucion_2 = (-coef_b - np.sqrt(discriminante)) / (2 * coef_a)
        
        # Devuelve la solución positiva
        if solucion_1 >= 0:
            return solucion_1
        elif solucion_2 >= 0:
            return solucion_2
        else:
            return None  # No hay soluciones positivas
    else:
        return None  # No hay soluciones reales
