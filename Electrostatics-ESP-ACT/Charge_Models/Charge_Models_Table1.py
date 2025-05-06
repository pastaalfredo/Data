#!/usr/bin/env python3

import os
import sys

def run_command(command):
    os.system(command)

def run_one(qtype: str, suffix: str = "") -> dict:
    molprops = "../AlexandriaFF/merged6.xml" #sapt-0.015.xml"
    log_filename = f"{qtype}{suffix}.log"
    base_command = f"alexandria train_ff -nooptimize -g {log_filename}"
    mps = { "15": "../AlexandriaFF/merged8.xml",
            "16": "../AlexandriaFF/merged8.xml" }
    suffix_commands = {
        "1": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/coul-p.xml -fc_elec 1",
        "2": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/coul-p.xml -fc_elec 1",
        "3": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/all-p.xml -fc_elec 1",
        "4": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/all-p.xml -fc_elec 1",
        "5": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/coul-g.xml -fc_elec 1",
        "6": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/coul-g.xml -fc_elec 1",
        "7": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/all-g.xml -fc_elec 1",
        "8": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/all-g.xml -fc_elec 1",
        "9": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/coul-gv.xml -fc_elec 1",
        "10": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/coul-gv.xml -fc_elec 1",
        "11": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/all-gv.xml -fc_elec 1",
        "12": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/all-gv.xml -fc_elec 1",
        "13": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/all-pg.xml -fc_elec 1",
        "14": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/all-pg.xml -fc_elec 1",
        "15": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/esp-gv.xml -fc_elec 1",
        "16": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/esp-gv.xml -fc_elec 1",
    }

    suffix_commands2 = {
        "1": " -sel ../Selection/ac-train.dat ",
        "2": " -sel ../Selection/ac-test.dat "
    }

    if qtype == "qACM" and suffix in suffix_commands:
        print(f"Running command for {qtype}{suffix} - COUL")
        mycmd = base_command + suffix_commands[suffix]
        if suffix in mps:
            mycmd += ( " -mp %s " % mps[suffix] )
        else:
            mycmd += ( " -mp %s " % molprops )
        run_command(mycmd)
    elif qtype != "qACM" and suffix in suffix_commands2:
        run_command(base_command + suffix_commands2[suffix]+ f"-charges ../AlexandriaFF/esp-paper-gaussian.xml -fc_elec 1  -ff ../ForceFields/GAFF.xml -qtype {qtype} -mp {molprops} ")

    print(f"Reading log file {log_filename} for COUL")
    mydict = {"COUL": {}, "ALLELEC": {}}
    with open(log_filename, "r") as inf:
        for line in inf:
            if "COULOMB (kJ" in line:
                words = line.strip().split()
                try:
                    mydict["COUL"]["N"] = int(words[2])
                    mydict["COUL"]["RMSD"] = round(float(words[5]),1)
                    mydict["COUL"]["MSE"] = round(float(words[6]),1)
                except ValueError:
                    sys.exit(f"Strange line {line.strip()}")

    if qtype == "qACM" and suffix in suffix_commands:
        print(f"Running command for {qtype}{suffix} - ALLELEC")
        mycmd = base_command + suffix_commands[suffix].replace("-fc_elec", "-fc_allelec")
        if suffix in mps:
            mycmd += ( " -mp %s " % mps[suffix] )
        else:
            mycmd += ( " -mp %s " % molprops )
        run_command(mycmd)
    elif qtype != "qACM" and suffix in suffix_commands2:
        run_command(base_command + suffix_commands2[suffix]+ f"-charges ../AlexandriaFF/esp-paper-gaussian.xml -fc_allelec 1  -ff ../ForceFields/GAFF.xml -qtype {qtype} -mp {molprops}")

    print(f"Reading log file {log_filename} for ALLELEC")
    with open(log_filename, "r") as inf:
        for line in inf:
            if "ALLELEC (kJ" in line:
                words = line.strip().split()
                try:
                    mydict["ALLELEC"]["N"] = int(words[2])
                    mydict["ALLELEC"]["RMSD"] = round(float(words[5]),1)
                    mydict["ALLELEC"]["MSE"] = round(float(words[6]),1)
                except ValueError:
                    sys.exit(f"Strange line {line.strip()}")

    return mydict

myqt = {
    "Mulliken": "Mulliken1955a",
    "Hirshfeld": "Hirshfeld1977a",
    "ESP": "Besler1990a",
    "CM5": "Marenich2012a",
    "BCC": "Jakalian2000a",
    "RESP": "Bayly1993a",
    "ACM": "ACT"
}

mytable = {}
couls = ""
allelecs = ""
labels = ""

charge_models = [
    ("header", "Existing charge models" ),
    ("Mulliken", "1"), ("Mulliken", "2"),
    ("Hirshfeld", "1"), ("Hirshfeld", "2"),
    ("ESP", "1"), ("ESP", "2"),
    ("CM5", "1"), ("CM5", "2"),
    ("BCC", "1"), ("BCC","2"),
    ("RESP","1"), ("RESP","2"),
    ("header", "Non-polarizable ACT models" ),
    ("ACM", "1"), ("ACM", "2"), ("ACM", "3"),
    ("ACM", "4"), ("ACM", "5"), ("ACM", "6"),
    ("ACM", "7"), ("ACM", "8"), ("ACM", "9"),
    ("ACM", "10"), ("ACM", "11"), ("ACM", "12"), 
    ("ACM", "15"), ("ACM", "16"),
    ("header", "Polarizable ACT model" ),
    ("ACM", "13"), ("ACM", "14")
]

nparams = { "1": 32, "3": 32, "5": 48, "7": 48, "9": 55, "11": 55, "13": 123 }

for qt, suffix in charge_models:
    if qt == "header":
        continue
    mytable[qt + suffix] = run_one("q" + qt, suffix)
    newcoul = f"COULOMB-{qt}{suffix}.xvg"
    run_command(f"mv COULOMB.xvg {newcoul}")
    couls += f" {newcoul}"
    newallelec = f"ALLELEC-{qt}{suffix}.xvg"
    run_command(f"mv ALLELEC.xvg {newallelec}")
    allelecs += f" {newallelec}"
    labels += f" {qt}{suffix}"
    for fn in [ "EXCHANGE.xvg", "DISPERSION.xvg", "INDUCTIONCORRECTION.xvg", "EPOT.xvg", "INDUCTION.xvg" ]:
        if os.path.exists(fn):
            os.unlink(fn)

run_command(f"viewxvg -f {couls} -label {labels} -ls None -mk -res -noshow -pdf legacy_coul.pdf")
run_command(f"viewxvg -f {allelecs} -label {labels} -ls None -mk -res -noshow -pdf legacy_allelec.pdf")

with open("legacy.tex", "w") as outf:
    outf.write("\\begin{table}[htb]\n")
    outf.write("\\centering\n")
    
    outf.write("\\caption{Root mean square deviation (RMSE) and mean signed error (MSE) of electrostatic energies (Elec, kJ/mol) and the sum of electrostatics and induction (Elec+Induc, kJ/mol) for popular charge models compared to SAPT2+(CCD)$\\delta$MP2 with the aug-cc-pVTZ basis set. The dataset consisted of 77 dimers (Table S5), and the number of data points (energies) is indicated as N. For the point charge (PC) and Gaussian charge (GC) models, the training targets are indicated, and RMSD and MSE values corresponding to the training set are indicated in bold. A non-polarizable model with virtual sites with a Gaussian distributed charge (on anions and potassium ion only) is labeled as GC+PGV. The polarizable point charge + Gaussian virtual site and shell (PC+GVS) model was trained on electrostatic and induction energies in one step. The \"test\" and \"train\" labels in the table for Mulliken, Hirshfeld, ESP, CM5, BCC, and RESP do not imply that we trained them. Instead, they indicate that we evaluated the models using compounds in the test or training sets (Table S5). The RMSD and MSE were calculated with respect to the SAPT electrostatic energy. \#P indicates the number of parameters in the model.}\n")
    
    
    outf.write("\\label{legacy}\n")
    outf.write("\\begin{tabular}{lccccccc}\n")
    outf.write("\\hline\n")
    outf.write(" & Target & \# P & N & \\multicolumn{2}{c}{Elec}  & \\multicolumn{2}{c}{Elec+Induc}\\\\\n")
    outf.write("Model & & & & RMSD & MSE & RMSD & MSE \\\\\n")

    label_map = {
        "ACM1": ("PC (Train)", "Elec"),
        "ACM2": ("PC (Test)", ""),
        "ACM3": ("PC (Train)", "Elec+Induc"),
        "ACM4": ("PC (Test)", ""),
        "ACM5": ("GC (Train)", "Elec"),
        "ACM6": ("GC (Test)", ""),
        "ACM7": ("GC (Train)", "Elec+Induc"),
        "ACM8": ("GC (Test)", ""),
        "ACM9": ("GC+PGV (Train)", "Elec"),
        "ACM10": ("GC+PGV (Test)", ""),
        "ACM11": ("GC+PGV (Train)", "Elec+Induc"),
        "ACM12": ("GC+PGV (Test)", ""),
        "ACM15": ("GC+PGV (Train)", "ESP"),
        "ACM16": ("GC+PGV (Test)", ""),
        "ACM13": ("PC+GVS (Train)", "Elec,Induc"),
        "ACM14": ("PC+GVS (Test)", "")
    }

    for qt, suffix in charge_models:
        if qt == "header":
            outf.write("\\hline\n")
            outf.write("\\multicolumn{8}{c}{\\bf %s}\\\\\n" % suffix)
            continue
        qt_with_suffix = qt + suffix
        cite = f"~\\cite{{{myqt[qt]}}}" if qt != "ACM" else ""
        label = label_map.get(qt_with_suffix, (qt, ""))[0]
        star = label_map.get(qt_with_suffix, (qt, ""))[1]

        coul_data = mytable[qt_with_suffix]["COUL"]
        elec_data = mytable[qt_with_suffix]["ALLELEC"]

        rmsd = 'RMSD'
        mse  = 'MSE'
        na   = 'N/A'
        if rmsd in coul_data:
            rmsd_coul_str = f"\\textbf{{{coul_data['RMSD']}}}" if "Train" in label else f"{coul_data['RMSD']}"
        else:
            rmsd_coul_str = na
        if mse in coul_data:
            mse_coul_str = f"\\textbf{{{coul_data['MSE']}}}" if "Train" in label else f"{coul_data['MSE']}"
        else:
            mse_coul_str = na
        if rmsd in elec_data:
            rmsd_elec_str = f"\\textbf{{{elec_data['RMSD']}}}" if "Train" in label else f"{elec_data['RMSD']}"
        else:
            rmsd_elec_str = na
        if mse in elec_data:
            mse_elec_str = f"\\textbf{{{elec_data['MSE']}}}" if "Train" in label else f"{elec_data['MSE']}"
        else:
            mse_elec_str = na
        np = ""
        if qt == "ACM" and suffix in nparams:
            np = nparams[suffix]
        outf.write(f"{label}{cite} & {star} & {coul_data['N']} & {np} & {rmsd_coul_str} & {mse_coul_str} & {rmsd_elec_str} & {mse_elec_str} \\\\\n")

    outf.write("\\hline\n")
    outf.write("\\end{tabular}\n")
    outf.write("\\end{table}\n")
