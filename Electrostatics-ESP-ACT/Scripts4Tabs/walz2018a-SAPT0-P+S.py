#!/usr/bin/env python3

import json, math
import numpy as np
from Electrostatic_4models import Point_core_1slater_2slater_shell

from Electrostatic_4models import Point_core_1slater_shell

one_4pi_eps0 = 1389

sapt = {
    "Li-F":  { "rmin": 1.564, "eelec":  -867.15},
    "Li-Cl": { "rmin": 2.021, "eelec":  -661.43},
    "Li-Br": { "rmin": 2.17, "eelec":   -612.35},
    "Na-F":  { "rmin":1.926, "eelec":   -748.97},
    "Na-Cl": { "rmin": 2.361, "eelec":  -607.68},
    "Na-Br": { "rmin": 2.502, "eelec":  -572.36},
    "K-F":   { "rmin":2.171, "eelec":   -704.64},
    "K-Cl":  { "rmin": 2.667, "eelec":  -566.07},
    "K-Br":  { "rmin": 2.821, "eelec":  -535.34}, }

walz = {
    "Li": { "zeta": 1.3772, "qcore": 3.391 },
    "Na": { "zeta": 1.2609, "qcore": 4.007 },
    "K":  { "zeta": 1.2412, "qcore": 5.464 },
    "F":  { "zeta": 0.9139, "qcore": 5.590 },
    "Cl": { "zeta": 0.6967, "qcore": 5.713 },
    "Br": { "zeta": 0.6465, "qcore": 5.778 } }

def gauss(distance:float, z1:float, z2:float)->float:
    zeta = z1*z1/math.sqrt(z1*z1+z2*z2)
    return one_4pi_eps0*math.erf(zeta*distance)/distance


def rmsd(ref_values, values):
    return np.sqrt(np.mean((np.array(ref_values) - np.array(values))**2))


def mse(ref_values, values):
    return np.mean(np.subtract(values,ref_values))


def print_walz(outf):


    with open('params_4_10.json', 'r') as json_f:
        params = json.load(json_f)


    with open('params_4_100.json', 'r') as json_f0:
        params0 = json.load(json_f0)


    rmsd_epc_values = []
    rmsd_eesp_values = []
    rmsd_eesp0_values = []

    mse_epc_values = []
    mse_eesp_values = []
    mse_eesp0_values = []

    for rm in sapt.keys():
        ions = rm.split("-")
        a1   = ions[0]
        a2   = ions[1]
        rrr  = sapt[rm]["rmin"]
        qc1  = walz[a1]["qcore"]
        qs1  = 1-qc1
        qc2  = walz[a2]["qcore"]
        qs2  = -1-qc2
        eelec = -gauss(rrr, walz[a1]["zeta"], walz[a2]["zeta"])
        # q_c_na, q_s1_na, q_s2_na, q_c_cl, q_s1_cl, q_s2_cl, z_c_na, z_c_cl, z_s_na, z_s_cl):

        #(distances, q_c_na, q_s_na, q_c_cl, q_s_cl, z_na, z_cl)
        # [eesp]  = Point_core_1slater_shell([rrr],
        #                                            params[a1]["q_c_"+a1+"_1"], params[a1]["q_s_"+a1+"_1"],
        #                                            params[a2]["q_c_"+a2+"_1"], params[a2]["q_s_"+a2+"_1"],
        #                                            params[a1]["z2_"+a1+"_1"],
        #                                            params[a2]["z2_"+a2+"_1"])
        #
        # [eesp0]  = Point_core_1slater_shell([rrr],
        #                                            params0[a1]["q_c_"+a1+"_1"], params0[a1]["q_s_"+a1+"_1"],
        #                                            params0[a2]["q_c_"+a2+"_1"], params0[a2]["q_s_"+a2+"_1"],
        #                                            params0[a1]["z2_"+a1+"_1"],
        #                                            params0[a2]["z2_"+a2+"_1"])

        [eesp]  = Point_core_1slater_2slater_shell([rrr],
                                                   params[a1]["q_c_"+a1+"_3"], params[a1]["q_s1_"+a1+"_3"], params[a1]["q_s2_"+a1+"_3"],
                                                   params[a2]["q_c_"+a2+"_3"], params[a2]["q_s1_"+a2+"_3"], params[a2]["q_s2_"+a2+"_3"],
                                                   params[a1]["z1_"+a1+"_3"], params[a1]["z2_"+a1+"_3"],
                                                   params[a2]["z1_"+a2+"_3"], params[a2]["z2_"+a2+"_3"])

        [eesp0]  = Point_core_1slater_2slater_shell([rrr],
                                                   params0[a1]["q_c_"+a1+"_3"], params0[a1]["q_s1_"+a1+"_3"], params0[a1]["q_s2_"+a1+"_3"],
                                                   params0[a2]["q_c_"+a2+"_3"], params0[a2]["q_s1_"+a2+"_3"], params0[a2]["q_s2_"+a2+"_3"],
                                                   params0[a1]["z1_"+a1+"_3"], params0[a1]["z2_"+a1+"_3"],
                                                   params0[a2]["z1_"+a2+"_3"], params0[a2]["z2_"+a2+"_3"])

        epc   = -one_4pi_eps0/rrr
        outf.write("%s & %g & %.0f & %.0f & %.0f & %.0f \\\\\n" % ( rm, rrr, sapt[rm]["eelec"], epc, eesp, eesp0) )
        rmsd_epc = rmsd([sapt[rm]["eelec"]], [epc])
        rmsd_eesp = rmsd([sapt[rm]["eelec"]], [eesp])
        rmsd_eesp0 = rmsd([sapt[rm]["eelec"]], [eesp0])

        mse_epc = mse([sapt[rm]["eelec"]], [epc])
        mse_eesp = mse([sapt[rm]["eelec"]], [eesp])
        mse_eesp0 = mse([sapt[rm]["eelec"]], [eesp0])

        rmsd_epc_values.append(rmsd_epc)
        rmsd_eesp_values.append(rmsd_eesp)
        rmsd_eesp0_values.append(rmsd_eesp0)

        mse_epc_values.append(mse_epc)
        mse_eesp_values.append(mse_eesp)
        mse_eesp0_values.append(mse_eesp0)

    mean_rmsd_epc = np.mean(rmsd_epc_values)
    mean_rmsd_eesp = np.mean(rmsd_eesp_values)
    mean_rmsd_eesp0 = np.mean(rmsd_eesp0_values)

    mean_mse_epc = np.mean(mse_epc_values)
    mean_mse_eesp = np.mean(mse_eesp_values)
    mean_mse_eesp0 = np.mean(mse_eesp0_values)

    outf.write("%s & & & %.1f & %.1f & %.1f \\\\\n" % ( "RMSD", mean_rmsd_epc, mean_rmsd_eesp, mean_rmsd_eesp0) )
    outf.write("%s & & & %.1f & %.1f & %.1f \\\\\n" % ( "MSE", mean_mse_epc, mean_mse_eesp, mean_mse_eesp0) )

if __name__ == "__main__":
    filename = "electable.tex"
    with open(filename, "w") as outf:
        outf.write("""\\begin{table}
        \\centering
        \\begin{tabular}{lccccc}
        \\hline
        Ion pair & r$_{min}$ & \\multicolumn{4}{c}{E$_{elec}$ (kJ/mol)}\\\\
        & ({\\AA})&                 SAPT & PC & ESP$_{vdw}$ & ESP$_0$\\\\
        \\hline
        """)
        print_walz(outf)
        outf.write("""\\hline
    \\end{tabular}
\\caption{Electrostatic energies from SAPT at the experimental minimum energy distance~\cite{NISTa}  and corresponding energy of two point charges. The
SAPT0~\cite{Sherill2013a} level of theory was applied and the total Electrostatics energy is reported. Point charge (PC) energy follows
from Coulomb's law. ESP indicates model consisting of a point charge, a 1S Slater distribution, and a 2S Slater distribution, fitted to the electrostatic potential from 0.0 or 2.0 to 4.5 {\AA}.}
    \\label{tab:sapt}
    \\end{table}
    """)

    print("Output written to %s" % filename)
