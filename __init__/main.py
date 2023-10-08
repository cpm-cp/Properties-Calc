from time import perf_counter
from tools.critic_values import critical_mixture_properties, molar_mass_mixture
from tools.ideal_properties import ideal_mix_properties
from tools.residual_properties import residual_properties_virial_equation

# Example substances and molar fractions (summing to 1) for the mixture
selected_substances = ['Carbon dioxide', 'Carbon monoxide', 'Hydrogen', 'Nitrogen', 'Methane', 'Water']
selected_molar_fractions = [0.136, 0.222, 0.1668, 0.4363, 0.0288, 0.0101]

# Reference values:
T_reference = 273.16 # K
P_reference = 6.117*10**-3 # bar

# States condition
    # Inner
T_inner = 1173.15 # K
P_inner = 550 # kPa
P_inner /= 100  # bar

    # Outet
T_outlet = 1173.15 # K
P_oulet = 1 # bar

# Save the start time:
start_time = perf_counter()

# Molecular mass by mixing:
Molecular_mass = molar_mass_mixture(substances=selected_substances, molar_fractions=selected_molar_fractions)

# Critical properties:
Tc_mixing, Vc_mixing, Zc_mixing, w_mixing, Pc_mixing = critical_mixture_properties(substances=selected_substances, molar_fractions=selected_molar_fractions)

# Redisual inner values at Reference values:
H_residual_inner_reference, S_residual_inner_reference = residual_properties_virial_equation(P=P_reference, P_critic=Pc_mixing, T=T_reference, T_critic=Tc_mixing, w_value=w_mixing)


# Ideal inner values:
H_ideal_inner, S_ideal_inner = ideal_mix_properties(substances=selected_substances, molar_fractions=selected_molar_fractions, T_reference=T_reference, T_state=T_inner, P_reference=P_reference, P_state=P_inner, current="entry")

# Residual inner values at Pressure and Temperature specific:
H_residual_inner_specific, S_residual_inner_specific = residual_properties_virial_equation(P=P_inner, P_critic=Pc_mixing, T=T_inner, T_critic=Tc_mixing, w_value=w_mixing)

# Residual Inner values:
    # Residual H inner:
H_residual_inner = H_residual_inner_specific - H_residual_inner_reference
    # Residual S inner
S_residual_inner = S_residual_inner_specific - S_residual_inner_reference


# Save the end time
end_time = perf_counter()

# Calculate the execution time
execution_time = end_time - start_time

# Print the time
print(f"Execution time: {execution_time:.4f} seconds")