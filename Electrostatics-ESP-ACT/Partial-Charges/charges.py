#!/usr/bin/env python
import os

compounds_of_interest = [
    "ethylammonium#chloride",
    "ammonium#potassium",
    "methylammonium#potassium",
    "acetate#potassium",
    "formate#potassium",
    "water#potassium"
]

def extract_v3bw_value(log_file):
    with open(log_file, 'r') as file:
        for line in file:
            if 'v3bw' in line and line.strip().startswith('|'):
                parts = [col.strip() for col in line.strip().split('|') if col.strip()]
                if len(parts) >= 3:
                    try:
                        return float(parts[2])
                    except ValueError:
                        continue
    return None

v3bw = extract_v3bw_value('tune_ff-elec-h.log')
v3bw_allG = extract_v3bw_value('tune_ff-allelec-n.log')
v3bw_elG = extract_v3bw_value('tune_ff-elec-v.log')

def extract_data_from_log(log_files):
    data = {
        "ESP": {compound: [] for compound in compounds_of_interest},
        "ACM-P": {compound: [] for compound in compounds_of_interest},
        "ACM-I": {compound: [] for compound in compounds_of_interest},
        "ACM-G": {compound: [] for compound in compounds_of_interest},
        "ACM-S": {compound: [] for compound in compounds_of_interest},
        "ACM-V": {compound: [] for compound in compounds_of_interest},
        "ACM-N": {compound: [] for compound in compounds_of_interest},
        "ACM-H": {compound: [] for compound in compounds_of_interest}
    }

    for log_file in log_files:
        if not os.path.isfile(log_file):
            print(f"File does not exist: {log_file}")
            continue

        if 'ESP' in log_file:
            file_type = "ESP"
        elif 'elec-p' in log_file:
            file_type = "ACM-P"
        elif 'allelec-i' in log_file:
            file_type = "ACM-I"
        elif 'elec-g' in log_file:
            file_type = "ACM-G"
        elif 'allelec-s' in log_file:
            file_type = "ACM-S"
        elif 'elec-v' in log_file:
            file_type = "ACM-V"
        elif 'allelec-n' in log_file:
            file_type = "ACM-N"
        elif 'elec-h' in log_file:
            file_type = "ACM-H"
        else:
            print(f"Unknown log file type: {log_file}")
            continue

        print(f"Processing {file_type} file: {log_file}")

        try:
            with open(log_file, 'r') as file:
                read_data = False
                current_compound = None
                core_shell_dict = {}

                for line in file:
                    line = line.strip()

                    if "Name:" in line:
                        for compound in compounds_of_interest:
                            if f"Name: {compound}" in line:
                                current_compound = compound
                                read_data = True
                                break
                        if current_compound:
                            continue  

                    if "EPOT" in line:
                        read_data = False
                    elif read_data:
                        columns = line.split()
                        if len(columns) > 4 and columns[0].isdigit():
                            atom_type = columns[2]
                            acm_value = columns[3]
                            if file_type in ["ESP", "ACM-P", "ACM-I", "ACM-G", "ACM-S", "ACM-V", "ACM-N"]:
                                data[file_type][current_compound].append((atom_type, acm_value, "", "", ""))
                            elif file_type == "ACM-H":
                                if "_s" in atom_type:
                                    atom_core_type = atom_type.split("_")[0]
                                    if atom_core_type in core_shell_dict:
                                        core_value = core_shell_dict[atom_core_type]
                                        try:
                                            core_value = float(core_value)
                                            shell_value = float(acm_value)
                                            total_value = core_value + shell_value
                                        except ValueError:
                                            core_value = "N/A"
                                            shell_value = "N/A"
                                            total_value = "N/A"
                                        data[file_type][current_compound].append((atom_core_type, core_value, shell_value, total_value))
                                        del core_shell_dict[atom_core_type] 
                                    else:
                                        print(f"Warning: Found shell value before core value for atom: {atom_type}")
                                else:
                                    core_shell_dict[atom_type] = acm_value

        except Exception as e:
            print(f"Error processing file {log_file}: {e}")

    return data

def save_data_as_latex(data, output_dir):
    combined_output_file = os.path.join(output_dir, 'combined_data.tex')

    with open(combined_output_file, 'w') as file:
        for compound in compounds_of_interest:
            compound2 = compound.split('#')[0]
            file.write(r"\begin{sidewaystable}")
            file.write("\n")
            file.write(r"\caption{Partial charges for " + compound2 + " from ESP and from ACT models, point charge (PC), Gaussian charge (GC), point core+Gaussian vsite (GC+PGV), and point charge + Gaussian vsite and shell (PC+GVS).  Partial charges for the PC, GC, and GC+PGV models trained on either electrostatic energy (e) or the sum of the electrostatic and induction energy (ei) from the SAPT2+(CCD)-$\\delta$MP2 method with the aug-cc-pVTZ basis set are reported. Partial charges for the PC+GVS model, trained on the electrostatic and induction energies are also provided.}")
            file.write("\n")
            file.write("\\hspace{-1cm}\n")
            file.write("\\begin{tabular}{lcccccccccccccccc}\n")
            file.write("\\hline\n")
            file.write(fr" Atom type & ESP & PC$_{{e}}$ & PC$_{{ei}}$ & GC$_{{e}}$ & GC$_{{ei}}$ & GC+PGV$_{{e}}$ & GC+PGV$_{{ei}}$ & \multicolumn{{3}}{{c}}{{PC+GVS}} \\\\")
            file.write("\n")
            file.write(" & & & & & & & & core & shell & total \\\\")
            file.write("\n")
            file.write("\\hline\n")

            esp_data = []
            for value in data["ESP"].get(compound, []):
                if isinstance(value, tuple):
                    esp_data.append(tuple(round(v, 2) if isinstance(v, (int, float)) else v for v in value))
                else:
                    esp_data.append(round(value, 2) if isinstance(value, (int, float)) else value)

            acm_p_data = []
            for value in data["ACM-P"].get(compound, []):
                if isinstance(value, tuple):
                    acm_p_data.append(tuple(round(v, 2) if isinstance(v, (int, float)) else v for v in value))
                else:
                    acm_p_data.append(round(value, 2) if isinstance(value, (int, float)) else value)

            acm_p2_data = []
            for value in data["ACM-I"].get(compound, []):
                if isinstance(value, tuple):
                    acm_p2_data.append(tuple(round(v, 2) if isinstance(v, (int, float)) else v for v in value))
                else:
                    acm_p2_data.append(round(value, 2) if isinstance(value, (int, float)) else value)

            acm_g_data = []
            for value in data["ACM-G"].get(compound, []):
                if isinstance(value, tuple):
                    acm_g_data.append(tuple(round(v, 2) if isinstance(v, (int, float)) else v for v in value))
                else:
                    acm_g_data.append(round(value, 2) if isinstance(value, (int, float)) else value)

            acm_g2_data = []
            for value in data["ACM-S"].get(compound, []):
                if isinstance(value, tuple):
                    acm_g2_data.append(tuple(round(v, 2) if isinstance(v, (int, float)) else v for v in value))
                else:
                    acm_g2_data.append(round(value, 2) if isinstance(value, (int, float)) else value)

            acm_v_data = []
            for value in data["ACM-V"].get(compound, []):
                if isinstance(value, tuple):
                    acm_v_data.append(tuple(round(v, 2) if isinstance(v, (int, float)) else v for v in value))
                else:
                    acm_v_data.append(round(value, 2) if isinstance(value, (int, float)) else value)

            acm_n_data = []
            for value in data["ACM-N"].get(compound, []):
                if isinstance(value, tuple):
                    acm_n_data.append(tuple(round(v, 2) if isinstance(v, (int, float)) else v for v in value))
                else:
                    acm_n_data.append(round(value, 2) if isinstance(value, (int, float)) else value)

            acm_h_data = []
            for value in data["ACM-H"].get(compound, []):
                if isinstance(value, tuple):
                    acm_h_data.append(tuple(round(v, 6) if isinstance(v, (int, float)) else v for v in value))
                else:
                    acm_h_data.append(round(value, 6) if isinstance(value, (int, float)) else value)



                    
            for item in esp_data:
                atom_type, esp_value, _, _, _ = item
                acm_p_value = next((val for val in acm_p_data if val[0] == atom_type), ("", "", ""))
                acm_p2_value = next((val for val in acm_p2_data if val[0] == atom_type), ("", "", ""))
                acm_g_value = next((val for val in acm_g_data if val[0] == atom_type), ("", "", ""))
                acm_g2_value = next((val for val in acm_g2_data if val[0] == atom_type), ("", "", ""))
                acm_v_value = next((val for val in acm_v_data if val[0] == atom_type), ("", "", ""))
                acm_n_value = next((val for val in acm_n_data if val[0] == atom_type), ("", "", ""))
                acm_h_value = next((val for val in acm_h_data if val[0] == atom_type), ("", "", "", ""))
                core_value = acm_h_value[1] if acm_h_value[1] != "N/A" else "N/A"
                shell_value = acm_h_value[2] if acm_h_value[2] != "N/A" else "N/A"
                total_value = acm_h_value[3] if acm_h_value[3] != "N/A" else "N/A"
                if atom_type in ["Cl-", "K+"] :
                    continue
                file.write(f"{atom_type} & {esp_value} & {acm_p_value[1]} & {acm_p2_value[1]} & {acm_g_value[1]} & {acm_g2_value[1]} & {acm_v_value[1]} & {acm_n_value[1]} & {core_value} & {shell_value} & {total_value} \\\\")
                file.write("\n")
            if compound == "water#potassium" and acm_h_data:
                 file.write(f" v3bw & 0 & 0 & 0 & 0 & 0 & {v3bw_allG} & {v3bw_elG} & {v3bw} & 0 & {v3bw} \\\\\n")    
            file.write("\\hline\n")
            file.write("\\end{tabular}\n")
            file.write("\\end{sidewaystable}\n")

    print("LaTeX file generated successfully!")

log_files = [
    "ESP.log", "tune_ff-elec-p.log", "tune_ff-allelec-i.log", "tune_ff-elec-g.log", 
    "tune_ff-allelec-s.log", "tune_ff-elec-v.log", "tune_ff-allelec-n.log", "tune_ff-elec-h.log"
]

output_directory = "./latex_output" 
os.makedirs(output_directory, exist_ok=True)

data = extract_data_from_log(log_files)
save_data_as_latex(data, output_directory)
