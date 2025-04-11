#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import json
import math
import os
from potential_elec_functions import *


T=100 #delta z= 0.01

def main():
   with open(f'../AnalyticalFitting/params_4_{T}.json', 'r') as json_f:
       params = json.load(json_f)


   cations = ["Li", "Na", "K"]
   anions = ["F", "Cl", "Br"]
   dir = 'SAPT_Alkali_Halides'

   files = [
       os.path.join(dir, f'distances_Electrostatics-{cation.lower()}{anion.lower()}.txt')
       for cation in cations
       for anion in anions
   ]

   data2 = [np.loadtxt(file) for file in files]
   distances = np.arange(1, 4.6, 0.1)
   K=1389

   function_values = {}

   func_index_to_name = {
       0: "P+G",
       1: "P+1S",
       2: "P+2Gs",
       3: "P+1S+2S"
   }

   func_index_to_function = {
       0: Point_core_gaussian_shell,
       1: Point_core_1slater_shell,
       2: Point_core_2gaussian_shell,
       3: Point_core_1slater_2slater_shell,
   }


   for file, dataset in zip(files, data2):
       fig, axes = plt.subplots(4, 1, figsize=(3, 9))

       for i in range(4):
           x, y = dataset[:, 0], dataset[:, 1]
           sorted_indices = np.argsort(x)
           x1 = x[sorted_indices]
           y1 = y[sorted_indices]

           for cation in cations:
               for anion in anions:
                   if file == f'SAPT_Alkali_Halides/distances_Electrostatics-{cation.lower()}{anion.lower()}.txt':
                       func_index = i
                       function_name = f"{func_index_to_name[func_index]}_{cation}{anion}"
                       function_values[function_name] = None




                       if func_index == 0 or func_index ==1:
                           if  func_index == 0:
                                 CM= Point_core_gaussian_shell
                           if  func_index == 1:
                                 CM= Point_core_1slater_shell

                           function_values[function_name] = CM(
                               distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"],
                               params[f"{cation}"][f"q_s_{cation}_{func_index}"],
                               params[f"{anion}"][f"q_c_{anion}_{func_index}"],
                               params[f"{anion}"][f"q_s_{anion}_{func_index}"],
                               params[f"{cation}"][f"z2_{cation}_{func_index}"],
                               params[f"{anion}"][f"z2_{anion}_{func_index}"]
                           )


                       elif func_index == 2 or func_index == 3:
                           if  func_index == 2:
                                 CM= Point_core_2gaussian_shell
                           if  func_index == 3:
                                 CM= Point_core_1slater_2slater_shell

                           function_values[function_name] = CM(
                               distances, params[f"{cation}"][f"q_c_{cation}_{func_index}"],
                               params[f"{cation}"][f"q_s1_{cation}_{func_index}"],
                               params[f"{cation}"][f"q_s2_{cation}_{func_index}"],
                               params[f"{anion}"][f"q_c_{anion}_{func_index}"],
                               params[f"{anion}"][f"q_s1_{anion}_{func_index}"],
                               params[f"{anion}"][f"q_s2_{anion}_{func_index}"],
                               params[f"{cation}"][f"z1_{cation}_{func_index}"],
                               params[f"{anion}"][f"z1_{anion}_{func_index}"],
                               params[f"{cation}"][f"z2_{cation}_{func_index}"],
                               params[f"{anion}"][f"z2_{anion}_{func_index}"]
                           )



                       axes[i].plot(x1, y1, label=f"SAPT_{cation}{anion}", color='r')
                       axes[i].plot(distances, function_values[function_name], color='black', label=function_name)
                       axes[i].set_xlabel('Distance ($\AA$)', fontsize=12)
                       axes[i].set_ylabel('Electrostatic energies (kJ/mol)', fontsize=8)
                       axes[i].tick_params(axis='x', labelsize=8)
                       axes[i].tick_params(axis='y', labelsize=8)
                       axes[i].legend(fontsize=12)
                       plt.tight_layout()
                       plt.savefig(f'SAPT_{anion}_{cation}_{T}.pdf')



       # for i in range(4):
       #     axes[i].set_xlabel('Distance ($\AA$)', fontsize=12)
       #     axes[i].set_ylabel('Electrostatic energies (kJ/mol)', fontsize=8)
       #     axes[i].tick_params(axis='x', labelsize=8)
       #     axes[i].tick_params(axis='y', labelsize=8)
       #     axes[i].legend(fontsize=12)
       #     plt.tight_layout()
       #     plt.savefig(f'SAPT_{anion}_{cation}.pdf')


       # plt.tight_layout()
       # plt.savefig(f'SAPT_{anion}_{cation}_0.1.pdf')
       plt.show()

if __name__ == "__main__":
    main()
