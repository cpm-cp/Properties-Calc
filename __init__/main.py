import time
from tools.critic_values import calculate_critical_mixture_properties
from tools.ideal_properties import calculate_H_ideal_mix, calculate_S_ideal_mix

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

# Save the start time
start_time = time.perf_counter()

Tc_mixing, Vc_mixing, Zc_mixing, w_mixing, Pc_mixing = calculate_critical_mixture_properties(substances=selected_substances, molar_fractions=selected_molar_fractions)
calculate_H_ideal_mix(substances=selected_substances, molar_fractions=selected_molar_fractions, T_reference=T_reference, T_state=T_inner, current="entry")
calculate_S_ideal_mix(substances=selected_substances, molar_fractions=selected_molar_fractions, T_reference=T_reference, T_state=T_inner, P_reference=P_reference, P_state=P_inner, current="entry")

# Save the end time
end_time = time.perf_counter()

# Calculate the execution time
execution_time = end_time - start_time

# Print the time
print(f"Execution time: {execution_time:.4f} seconds")