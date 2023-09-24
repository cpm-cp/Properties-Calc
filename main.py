from critic_values import calculate_mixture_properties
from cp_ideal_mix import calculate_H_ideal_mix

# Example substances and molar fractions (summing to 1) for the mixture
selected_substances = ['Carbon dioxide', 'Carbon monoxide', 'Hydrogen', 'Nitrogen', 'Methane', 'Water']
selected_molar_fractions = [0.136, 0.222, 0.1668, 0.4363, 0.0288, 0.0102]
Tmin = 273.16 # K
Tmax = 1173.15 # K

calculate_mixture_properties(selected_substances, selected_molar_fractions)
calculate_H_ideal_mix(selected_substances, selected_molar_fractions, Tmin, Tmax)