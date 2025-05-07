#!/usr/bin/env python3

import os, glob, sys

def run_command(command):
    os.system(command)

merged8 = "../AlexandriaFF/merged8.xml"
train   = "Train"
test    = "Test"
acmparm = {
    "21": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/coul-p.xml -fc_elec 1", "nparm": 32, "label": ("PC", "Elec") },
    "22": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/coul-p.xml -fc_elec 1", "label": ("PC", "") },
    "23": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/all-p.xml -fc_elec 1", "nparm": 32, "label": ("PC", "Elec+Induc") },
    "24": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/all-p.xml -fc_elec 1", "label": ("PC", "") },
    "5": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/coul-g.xml -fc_elec 1", "nparm": 48, "label": ("GC", "Elec") },
    "6": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/coul-g.xml -fc_elec 1", "label": ("GC", "") },
    "7": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/all-g.xml -fc_elec 1", "nparm": 48, "label": ("GC", "Elec+Induc") },
    "8": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/all-g.xml -fc_elec 1", "label": ("GC", "") },
    "9": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/coul-gv.xml -fc_elec 1", "nparm": 54, "label": ("GC+PGV", "Elec") },
    "10": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/coul-gv.xml -fc_elec 1", "label": ("GC+PGV", "") },
    "11": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/all-gv.xml -fc_elec 1", "nparm": 54, "label": ("GC+PGV", "Elec+Induc") },
    "12": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/all-gv.xml -fc_elec 1", "label": ("GC+PGV", "") },
    "13": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/all-pg.xml -fc_elec 1", "nparm": 123, "label": ("PC+GVS", "Elec,Induc") },
    "14": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/all-pg.xml -fc_elec 1", "label": ("PC+GVS", "") },
    "15": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/esp-g.xml -fc_elec 1", "mp": merged8, "nparm": 48, "label": ("GC", "ESP") },
    "16": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/esp-g.xml -fc_elec 1", "mp": merged8, "label": ("GC", "") },
    "17": { "cmd": " -sel ../Selection/ac-train.dat -ff ../AlexandriaFF/esp-gv.xml -fc_elec 1", "mp": merged8, "nparm": 54, "label": ("GC+PGV", "ESP") },
    "18": { "cmd": " -sel ../Selection/ac-test.dat -ff ../AlexandriaFF/esp-gv.xml -fc_elec 1", "mp": merged8, "label": ("GC+PGV", "") }
    }

def run_one(qtype: str, suffix: str = "") -> dict:
    molprops = "../AlexandriaFF/merged6.xml" #sapt-0.015.xml"
    log_filename = f"{qtype}{suffix}.log"
    base_command = f"alexandria train_ff -nooptimize -g {log_filename}"

    suffix_commands2 = {
        "1": " -sel ../Selection/ac-train.dat ",
        "2": " -sel ../Selection/ac-test.dat "
    }

    if qtype == "qACM" and suffix in acmparm:
        print(f"Running command for {qtype}{suffix} - COUL")
        mycmd = base_command + acmparm[suffix]["cmd"]
        if "mp" in acmparm[suffix]:
            mycmd += ( " -mp %s " % acmparm[suffix]["mp"] )
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

    if qtype == "qACM" and suffix in acmparm:
        print(f"Running command for {qtype}{suffix} - ALLELEC")
        mycmd = base_command + acmparm[suffix]["cmd"].replace("-fc_elec", "-fc_allelec")
        if "mp" in acmparm[suffix]:
            mycmd += ( " -mp %s " % acmparm[suffix]["mp"] )
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
    ("header", "Non-polarizable ESP-based ACT models" ),
    ("ACM", "15"), ("ACM", "16"),
    ("ACM", "17"), ("ACM", "18"),
    ("header", "Non-polarizable SAPT-based ACT models" ),
    ("ACM", "21"), ("ACM", "22"), ("ACM", "23"), ("ACM", "24"),
    ("ACM", "5"), ("ACM", "6"), ("ACM", "7"), ("ACM", "8"), 
    ("ACM", "9"), ("ACM", "10"), ("ACM", "11"), ("ACM", "12"), 
    ("header", "Polarizable SAPT-based ACT model" ),
    ("ACM", "13"), ("ACM", "14")
]


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
    for fn in [ "EXCHANGE.xvg", "EXCHIND.xvg", "DISPERSION.xvg", "INDUCTIONCORRECTION.xvg", "EPOT.xvg", "INDUCTION.xvg" ] + glob.glob("#*#"):
        if os.path.exists(fn):
            os.unlink(fn)

run_command(f"viewxvg -f {couls} -label {labels} -ls None -mk -res -noshow -pdf legacy_coul.pdf")
run_command(f"viewxvg -f {allelecs} -label {labels} -ls None -mk -res -noshow -pdf legacy_allelec.pdf")

with open("legacy.tex", "w") as outf:
    outf.write("\\begin{table}[htb]\n")
    outf.write("\\centering\n")
    
    outf.write("\\caption{Root mean square deviation (RMSE) and mean signed error (MSE) of electrostatic energies (Elec, kJ/mol) and the sum of electrostatics and induction (Elec+Induc, kJ/mol) for popular charge models compared to SAPT2+(CCD)$\\delta$MP2 with the aug-cc-pVTZ basis set. The dataset consisted of 77 dimers (Table S5), with 8753 data points for training and 2729 for testing (for ESP training see Methods). \\#P indicates the number of parameters in the model. The training targets are indicated and RMSD and MSE values that correspond to the training set are indicated in {\\bf bold font}. A non-polarizable model with virtual sites with a Gaussian distributed charge (on anions and potassium ion only) is labeled as GC+PGV. The polarizable point charge + Gaussian virtual site and shell (PC+GVS) model was trained on electrostatic and induction energies in one step. The \"test\" and \"train\" labels in the table for Mulliken, Hirshfeld, ESP, CM5, BCC, and RESP do not imply that we trained them. Instead, they indicate that we evaluated the models using compounds in the test or training sets (Table S5). The RMSD and MSE were calculated with respect to the SAPT electrostatic energy.}\n")
    
    
    outf.write("\\label{legacy}\n")
    outf.write("\\begin{tabular}{lcccccccc}\n")
    outf.write("\\hline\n")
    outf.write(" & Dataset & Training & \\#P & \\multicolumn{2}{c}{Elec}  & \\multicolumn{2}{c}{Elec+Induc}\\\\\n")
    outf.write("Model & & target & & RMSD & MSE & RMSD & MSE \\\\\n")


    for qt, suffix in charge_models:
        if qt == "header":
            outf.write("\\hline\n")
            outf.write("\\multicolumn{8}{c}{\\bf %s}\\\\\n" % suffix)
            continue
        qt_with_suffix = qt + suffix
        cite = f"~\\cite{{{myqt[qt]}}}" if qt != "ACM" else ""
        label = qt
        star  = None
        if suffix in acmparm:
            label = acmparm[suffix]["label"][0]
            star  = acmparm[suffix]["label"][1]

        rmsd = 'RMSD'
        mse  = 'MSE'
        na   = 'N/A'
        rmsd_str = {}
        mse_str  = {}
        for mydata in [ "COUL", "ALLELEC" ]:
            ttable = mytable[qt_with_suffix][mydata]
            if rmsd in ttable and mse in ttable:
                bold    = False
                dataset = train
                if int(suffix) % 2 == 0:
                    dataset = test
                if star and "nparm" in acmparm[suffix]:
                    if "Elec,Induc" in star:
                        bold = True
                    elif "Elec+Induc" in star:
                        bold = mydata == "ALLELEC"
                    else:
                        bold = mydata == "COUL"
                        
                if  bold:
                    rmsd_str[mydata] = f"\\textbf{{{ttable[rmsd]}}}"
                    mse_str[mydata]  = f"\\textbf{{{ttable[mse]}}}"
                else:
                    rmsd_str[mydata] = f"{ttable[rmsd]}"
                    mse_str[mydata]  = f"{ttable[mse]}"
            else:
                print("Something wrong with table for %s" % qt_with_suffix)
                sys.exit(ttable)
        np = ""
        if qt == "ACM" and "nparm" in acmparm[suffix]:
            np = acmparm[suffix]["nparm"]
        N = mytable[qt_with_suffix][mydata]["N"]
        target = ""
        if star:
            target = star
        outf.write(f"{label}{cite} & {dataset} &{target} & {np} & {rmsd_str['COUL']} & {mse_str['COUL']} & {rmsd_str['ALLELEC']} & {mse_str['ALLELEC']} \\\\\n")

    outf.write("\\hline\n")
    outf.write("\\end{tabular}\n")
    outf.write("\\end{table}\n")
