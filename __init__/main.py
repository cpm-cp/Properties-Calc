import numpy as np
from time import perf_counter
from tools.critic_values import critical_mixture_properties, molar_mass_mixture
from tools.ideal_properties import ideal_mix_properties
from tools.residual_properties import residual_properties_virial_equation, residual_properties_Van_der_Waals

# Example substances and molar fractions (summing to 1) for the mixture
selected_substances = ['Carbon dioxide', 'Carbon monoxide', 'Hydrogen', 'Nitrogen', 'Methane', 'Water']
selected_molar_fractions = np.array([0.136, 0.222, 0.1668, 0.4363, 0.0288, 0.0101])

if sum(selected_molar_fractions) != 1:
    raise Exception("The sum by molar fraction is more than 1, verify.")

# Reference values:
T_reference = 273.16 # K
P_reference = 6.117*10**-3 # bar

# States condition
Temperature = np.array([1173.15, 870.45]) # inner and outer in  Kelvin
Pressure = np.array([5.5, 1]) # innner and outer in bar

# Heat losses:
Q_loss = 50 # kW : kJ/s

# Mass flow rate:
m_point = 0.6242 # kg/s

# Save the start time:
start_time = perf_counter()

# Molecular mass by mixing:
Molecular_mass = molar_mass_mixture(substances=selected_substances, molar_fractions=selected_molar_fractions)

# Molar flow rate:
n_point = m_point / Molecular_mass

# Critical properties:
Tc_mixing, Vc_mixing, Zc_mixing, w_mixing, Pc_mixing = critical_mixture_properties(substances=selected_substances, molar_fractions=selected_molar_fractions)

# Redisual inner values at Reference values:
H_residual_reference, S_residual_reference = residual_properties_virial_equation(P=P_reference, P_critic=Pc_mixing, T=T_reference, T_critic=Tc_mixing, w_value=w_mixing)

def EoS(select_EoS:str) -> tuple:

    # Ideal inner and outer values:
    H_ideal_inner, S_ideal_inner = ideal_mix_properties(substances=selected_substances, molar_fractions=selected_molar_fractions, T_reference=T_reference, T_state=Temperature[0], P_reference=P_reference, P_state=Pressure[0], current="entry")
    H_ideal_out, S_ideal_out = ideal_mix_properties(substances=selected_substances, molar_fractions=selected_molar_fractions, T_reference=T_reference, T_state=Temperature[1], P_reference=P_reference, P_state=Pressure[1], current="out")

    if select_EoS == "Ideal gas":

        H_1, H_2 = H_ideal_inner, H_ideal_out
        S_1, S_2 = S_ideal_inner, S_ideal_out

    elif select_EoS == "Virial":
        # Residual innera nd outer values at Pressure and Temperature specific:
        # Virial equation.
        H_residual_inner_virial, S_residual_inner_virial = residual_properties_virial_equation(P=Pressure[0], P_critic=Pc_mixing, T=Temperature[0], T_critic=Tc_mixing, w_value=w_mixing)
        H_residual_out_virial, S_residual_out_virial = residual_properties_virial_equation(P=Pressure[1], P_critic=Pc_mixing, T=Temperature[1], T_critic=Tc_mixing, w_value=w_mixing)

        # Residual Inner values:
        H_residual_inner = H_residual_inner_virial - H_residual_reference
        H_residual_out = H_residual_out_virial  - H_residual_reference
        S_residual_inner = S_residual_inner_virial - S_residual_reference
        S_residual_out = S_residual_out_virial - S_residual_reference

        H_1, H_2 = (H_ideal_inner + H_residual_inner), (H_ideal_out + H_residual_out)
        S_1, S_2 = (S_ideal_inner + S_residual_inner), (S_ideal_out + S_residual_out)

    elif select_EoS == "Van der Waals":    
        # Residual innera nd outer values at Pressure and Temperature specific:
        # Virial equation.
        H_residual_VdW, S_residual_VdW = residual_properties_Van_der_Waals(substances=selected_substances, molar_fraction=selected_molar_fractions, P=Pressure, T=Temperature)

        H_1, H_2 = (H_ideal_inner + H_residual_VdW[0]), (H_ideal_out + H_residual_VdW[1]) 
        S_1, S_2 = (S_ideal_inner + S_residual_VdW[0]), (S_ideal_out + S_residual_VdW[1])
    else:
        raise Exception("EoS select not find, verify.")
    
    return H_1, H_2, S_1, S_2

H_1, H_2, S_1, S_2 = EoS(select_EoS="Van der Waals")

# Potencia de la turbina:
W_t = n_point*(H_1 - H_2 - Q_loss)

# Save the end time
end_time = perf_counter()

# Calculate the execution time
execution_time = (end_time - start_time) * 1000  # miliseconds

print(H_1)
print(H_2)

print(f"the molar flow rate is {Molecular_mass} kmol/s")
print(f"The power of turbine is: {W_t} kW")
print(f"Tiempo de ejecución: {execution_time:.5f} milisegundos")