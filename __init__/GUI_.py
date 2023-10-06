import tkinter as tk
from tkinter import ttk
from customtkinter import CTk 
from time import perf_counter
from tools.critic_values import calculate_critical_mixture_properties
from tools.ideal_properties import calculate_ideal_mix_properties
from tools.residual_properties import residual_properties

root = CTk()

substances = ['Carbon dioxide', 'Carbon monoxide', 'Hydrogen', 'Nitrogen', 'Methane', 'Water']

def calculate_properties():
    selected_substances = [substance_comboboxes[i].get() for i in range(num_substances)]
    selected_molar_fractions = [float(fraction_entries[i].get()) for i in range(num_substances)]

    T_reference = float(T_reference_entry.get())
    P_reference = float(P_reference_entry.get())
    T_inner = float(T_inner_entry.get())
    P_inner = float(P_inner_entry.get()) / 100  # Converting kPa to bar

    start_time = perf_counter()

    Tc_mixing, Vc_mixing, Zc_mixing, w_mixing, Pc_mixing = calculate_critical_mixture_properties(substances=selected_substances, molar_fractions=selected_molar_fractions)
    H_residual_inner_reference, S_residual_inner_reference = residual_properties(P=P_reference, P_critic=Pc_mixing, T=T_reference, T_critic=Tc_mixing, w_value=w_mixing)
    H_ideal_inner, S_ideal_inner = calculate_ideal_mix_properties(substances=selected_substances, molar_fractions=selected_molar_fractions, T_reference=T_reference, T_state=T_inner, P_reference=P_reference, P_state=P_inner, current="entry")
    H_residual_inner_specific, S_residual_inner_specific = residual_properties(P=P_inner, P_critic=Pc_mixing, T=T_inner, T_critic=Tc_mixing, w_value=w_mixing)
    H_residual_inner = H_residual_inner_specific - H_residual_inner_reference
    S_residual_inner = S_residual_inner_specific - S_residual_inner_reference
    H_inner = H_ideal_inner + H_residual_inner
    S_inner = S_ideal_inner + S_residual_inner

    result_label.config(text=f"The value for the Enthalpy is: {H_inner} kJ/mol\nThe value for the Entropy is: {S_inner} kJ/mol*K")

    end_time = perf_counter()
    execution_time = end_time - start_time
    time_label.config(text=f"Execution time: {execution_time:.4f} seconds")

# Creating GUI elements
num_substances_label = tk.Label(root, text="Number of Substances:")
num_substances_label.pack()
num_substances_var = tk.StringVar()
num_substances_entry = tk.Entry(root, textvariable=num_substances_var)
num_substances_entry.pack()

confirm_button = tk.Button(root, text="Confirm", command=lambda: create_substance_inputs(int(num_substances_var.get())))
confirm_button.pack()

substance_comboboxes = []
fraction_entries = []
num_substances = 0

def create_substance_inputs(num):
    global num_substances
    num_substances = num
    for widget in root.winfo_children():
        if widget not in [num_substances_label, num_substances_entry, confirm_button]:
            widget.destroy()
    global substance_comboboxes
    global fraction_entries
    substance_comboboxes = []
    fraction_entries = []
    for i in range(num):
        label = tk.Label(root, text=f"Substance {i + 1}:")
        label.pack()
        substance_var = tk.StringVar()
        substance_combobox = ttk.Combobox(root, textvariable=substance_var, values=substances)
        substance_combobox.pack()
        substance_comboboxes.append(substance_combobox)

        fraction_label = tk.Label(root, text=f"Molar Fraction {i + 1}:")
        fraction_label.pack()
        fraction_var = tk.StringVar()
        fraction_entry = tk.Entry(root, textvariable=fraction_var)
        fraction_entry.pack()
        fraction_entries.append(fraction_var)

    T_reference_label = tk.Label(root, text="Reference Temperature (K):")
    T_reference_label.pack()
    global T_reference_entry
    T_reference_entry = tk.Entry(root)
    T_reference_entry.pack()

    P_reference_label = tk.Label(root, text="Reference Pressure (bar):")
    P_reference_label.pack()
    global P_reference_entry
    P_reference_entry = tk.Entry(root)
    P_reference_entry.pack()

    T_inner_label = tk.Label(root, text="Inner Temperature (K):")
    T_inner_label.pack()
    global T_inner_entry
    T_inner_entry = tk.Entry(root)
    T_inner_entry.pack()

    P_inner_label = tk.Label(root, text="Inner Pressure (kPa):")
    P_inner_label.pack()
    global P_inner_entry
    P_inner_entry = tk.Entry(root)
    P_inner_entry.pack()

    calculate_button = tk.Button(root, text="Calculate", command=calculate_properties)
    calculate_button.pack()

    global result_label
    result_label = tk.Label(root, text="")
    result_label.pack()

    global time_label
    time_label = tk.Label(root, text="")
    time_label.pack()

root.mainloop()
