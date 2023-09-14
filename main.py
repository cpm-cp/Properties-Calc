from critic_values import calculate_mixture_properties
from cp_ideal_mix import calculate_H_ideal_mix

# Example substances and molar fractions (summing to 1) for the mixture
selected_substances = ['Nitrogen', 'Carbon dioxide', 'Water']
selected_molar_fractions = [0.7, 0.1, 0.2]
Tmin = 273.16 # K
Tmax = 500 # K

calculate_mixture_properties(selected_substances, selected_molar_fractions)
calculate_H_ideal_mix(selected_substances, selected_molar_fractions, Tmin, Tmax)