#!/usr/bin/env python
import numpy as np
import json



with open('data_SC_ION_G.json', 'r') as file:
    data_SC_ION = json.load(file)





def rmsd(data_SC_ION, ref_key1, ref_key2):
    ref1_values = [values for values in data_SC_ION[ref_key1]]
    ref2_values = [values for values in data_SC_ION[ref_key2]]

    if len(ref1_values) != len(ref2_values):
        raise ValueError("Mismatched number of values between ref_key1 and ref_key2")

    ref1_values = np.array(ref1_values)
    ref2_values = np.array(ref2_values)

    return np.sqrt(np.mean((ref1_values - ref2_values)**2))

rmsd_value1 = rmsd(data_SC_ION, 'E$_{elec}$ (kJ/mol) GAFF_RESP', 'E$_{elec}$ (kJ/mol) SAPT')
rmsd_value2 = rmsd(data_SC_ION, 'E$_{elec}$ (kJ/mol) ACT$_{S}$', 'E$_{elec}$ (kJ/mol) SAPT')
rmsd_value3 = rmsd(data_SC_ION, 'E$_{elec}$ (kJ/mol) GAFF_BCC', 'E$_{elec}$ (kJ/mol) SAPT')
rmsd_value22 = rmsd(data_SC_ION, 'E$_{elec}$ (kJ/mol) ACT$_{GC}$', 'E$_{elec}$ (kJ/mol) SAPT')

def mse(data_SC_ION, ref_key1, ref_key2):
    ref1_values = [values for values in data_SC_ION[ref_key1]]
    ref2_values = [values for values in data_SC_ION[ref_key2]]

    MSE = np.mean(np.subtract(ref1_values,ref2_values))
    return MSE

mse_value1 = mse(data_SC_ION, 'E$_{elec}$ (kJ/mol) GAFF_RESP', 'E$_{elec}$ (kJ/mol) SAPT')
mse_value2 = mse(data_SC_ION, 'E$_{elec}$ (kJ/mol) ACT$_{S}$', 'E$_{elec}$ (kJ/mol) SAPT')
mse_value3 = mse(data_SC_ION, 'E$_{elec}$ (kJ/mol) GAFF_BCC', 'E$_{elec}$ (kJ/mol) SAPT')
mse_value22 = mse(data_SC_ION, 'E$_{elec}$ (kJ/mol) ACT$_{GC}$', 'E$_{elec}$ (kJ/mol) SAPT')

file_path = "AC-MA-IONS-GAFF.tex"

with open(file_path, "w") as file:

    file.write("\\begin{table}[ht]\n")
    file.write("\\centering\n")
    file.write("\\caption{Electrostatic energy (kJ/mol) between alkali halides or water (oxygen) and amino acid side chain analogs, formate (oxygen), acetate (oxygen), methylammonium (nitrogen), ethylammonium (nitrogen) from SAPT2+(CCD)$\\delta$MP2/aug-cc-pVTZ, and charges determined using either RESP~\\cite{Bayly1993a} or BCC~\\cite{Jakalian2000a}, as well as two models generated using the ACT.}\n")
    file.write("\label{tab:ac_ma_ions}\n")
    file.write("\\begin{tabular}{lcccccc} \n")
    file.write("\\hline \n")
    file.write("Ion & r$_{min}$ & SAPT & RESP & BCC & GC+PGV & PC+GVS \\\\\n")
    file.write("\\hline \n")


    for i in range(len(data_SC_ION['Ion'])):
        file.write(f"{data_SC_ION['Ion'][i]} & {data_SC_ION['r$_{min}$'][i]} & {data_SC_ION['E$_{elec}$ (kJ/mol) SAPT'][i]:.1f} & {data_SC_ION['E$_{elec}$ (kJ/mol) GAFF_RESP'][i]:.1f} & {data_SC_ION['E$_{elec}$ (kJ/mol) GAFF_BCC'][i]:.1f} & {data_SC_ION['E$_{elec}$ (kJ/mol) ACT$_{GC}$'][i]:.1f} & {data_SC_ION['E$_{elec}$ (kJ/mol) ACT$_{S}$'][i]:.1f} \\\\\n")
    file.write("\\hline\n")
    file.write(f"RMSD & & & {rmsd_value1:.1f} & {rmsd_value3:.1f} & {rmsd_value22:.1f} & {rmsd_value2:.1f} \\\\\n")
    file.write(f"MSE & & & {mse_value1:.1f} & {mse_value3:.1f} & {mse_value22:.1f} & {mse_value2:.1f} \\\\\n")

    file.write("\\hline \n")
    file.write("\\end{tabular} \n")
    file.write("\n")
    file.write("\\end{table}")

print("Please check %s" % file_path)
