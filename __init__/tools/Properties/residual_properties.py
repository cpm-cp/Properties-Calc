import pandas as pd
import numpy as np
from tools.critic_values import extract_critical_properties_4_VdW
from tools.critic_values import extract_critical_properties_SRK
from tools.mixing_rules import MixRules4VdW, MixRules4SRK
from tools.EoS.Volumes import calc_VdW_volume, calc_SRK_volume

def residual_properties_virial_equation(P:float, P_critic:float, T:float, T_critic:float, w_value:float, R: float = 8.314) -> tuple[float, float]:
    """Calculate the residual enthalpy and entropy in function to the pressure, critic pressure, temperature, critic temperature and acentric value.

    At moment to obtain the function values is necessary obtain this specific order:
    - Residual enthalpy.
    - Residual entropy.

    Args:
        P (float): Pressure in bar.
        P_critic (float): Critic pressure in bar.
        T (float): Temperature in Kelvin.
        T_critic (float): Critic temperature in Kelvin
        w_value (float): Acentric value (dimensionaless)
        R (float, optional): Gas constant in kJ/kmol * Kelvin. Defaults to 8.314.

    Returns:
        float: Residual enthalpy and residual entropy.
    """
    b_0 = 0.083 - (0.422 / (T / T_critic)**1.6)
    
    b_1 = 0.139 - (0.172 / (T / T_critic)**4.2)

    db_0dt = 0.675 / (T / T_critic)**2.6

    db_1dt = 0.722 / (T / T_critic)**5.2
    
    # Use to calc H_residual:
    h_residual = R*T_critic * ((P / P_critic) * (b_0 - (T / T_critic) * db_0dt + w_value * (b_1 - (T / T_critic) * db_1dt)))

    # Use to cal S_esidual:
    s_residual = R * (-(P / P_critic) * (db_0dt + w_value * db_1dt))
    
    return h_residual, s_residual


def residual_properties_Van_der_Waals(substances:list[str], molar_fraction, P:list[float], T:list[float], R: float = 83.14) -> tuple[list, list]:
    list_susbstances = substances
    molar_fract = molar_fraction 
    Tc_values, Pc_values = extract_critical_properties_4_VdW(list_susbstances)
    

    # Calc a and b paramter to Van der Waals:
    mix_rules = MixRules4VdW(critical_temperature=Tc_values, critical_pressure=Pc_values)
    a = mix_rules.calc_a_VdW(molar_fractions=molar_fract)
    b = mix_rules.calc_b_VdW(molar_fractions=molar_fract)
    volume = calc_VdW_volume(substances=list_susbstances, molar_fraction=molar_fract , Pressure=P, Temperature=T)
    Z = (P * volume) / (R * T)
    B = (b * P) / (R * T)

    H_r = R*T * (Z - 1) - (a / volume)
    S_r = R * np.log(Z - B)
    return H_r, S_r

def residual_properties_SRK(substances:list[str], molar_fraction, P:list[float], T:list[float], R: float = 8.314) -> tuple[list, list]:
    list_substances = substances
    molar_fract = molar_fraction
    Tc_values, Pc_values, w_values = extract_critical_properties_SRK(list_substances)

    # Calc a and b parameter to Soave-Redlich-Kwong (SRK):
    m_values = 0.480 + 1.574 * w_values - 0.176 * w_values**2
    m = np.tile(m_values, (2, 1))
    mix_rules = MixRules4SRK(critical_temperature=Tc_values, critical_pressure=Pc_values)
    a_ij = mix_rules.calc_a_ij(m=m, temperature=T)

    a = mix_rules.calc_a_SRK(molar_fractions=molar_fract, a_ij=a_ij)
    b = mix_rules.calc_b_SRK(molar_fractions=molar_fract)
    volume = calc_SRK_volume(substances=list_substances, molar_fraction=molar_fract, Pressure=P, Temperature=T)
    Z = (P * volume) / (R * T)
    B = (b * P) / (R * T)

    H_r = R*T * (Z - 1) - (a / volume)
    S_r = R * np.log(Z - B)
    
    return H_r, S_r

