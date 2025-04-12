#!/usr/bin/env python3
import os

output_dir = "PARAMETERS"
os.makedirs(output_dir, exist_ok=True)

latex_labels = {
    "a1dexp": "Induction correction A",
    "bdexp": "Induction correction b",
    "alpha": "Polarizability $\\alpha$",
    "charge": "Charge",
    "chi": "Electronegativity $\\chi$",
    "delta_chi": "Bond electronegativity $\\Delta\\chi$",
    "delta_eta": "Bond hardness $\\Delta\\eta$",
    "eta": "Hardness $\\eta$",
    "zeta": "Gaussian distribution width $\\zeta$",
    "vs3sa": "Virtual site position $a$"
}

def extract_table_from_log(log_file):
    table = []
    with open(log_file, 'r') as file:
        start_reading = False
        for line in file:
            if "Here are the best parameters I found, together with some summary statistics of the last population:" in line:
                start_reading = True
                next(file)  
                continue
            elif "Final best genome." in line:
                break
            elif start_reading:
                if not line.strip():
                    continue
                if line.startswith('|'):
                    row = line.strip().split('|')
                    row = [cell.strip() for cell in row]
                    table.append(row[1:4])
    return table

def convert_labels(table):
    new_table = []
    for row in table:
        label = latex_labels.get(row[0], row[0])
        new_table.append([label, f"\\verb^{row[1]}^", row[2]])
    return new_table

def save_table_as_latex(table, output_file, log_path):
    caption_label = os.path.basename(log_path).replace('_', ' ').replace('.log', '')
    with open(output_file, 'w') as file:
        file.write("\\begin{table}[ht]\n")
        file.write(f"\\caption{{Parameters after training. For details please see the ACT manual.}}\n")
        file.write("\\begin{tabular}{lcc}\n")
        file.write("\\hline\n")
        file.write("Parameter & Atom type & Value \\\\ \n")
        file.write("\\hline\n")
        for row in table:
            file.write(" & ".join(row) + " \\\\ \n")
        file.write("\\hline\n")
        file.write("\\end{tabular}\n")
        file.write("\\end{table}\n")

for subdir in ["Elec", "AllElec"]:
    if not os.path.isdir(subdir):
        continue
    for root, dirs, files in os.walk(subdir):
        for file in files:
            if file.endswith('.log'):
                log_path = os.path.join(root, file)
                table = extract_table_from_log(log_path)
                if table:
                    table = convert_labels(table)
                    out_name = f"{os.path.basename(root)}_{file.replace('.log', '')}_{subdir}.tex"
                    output_file = os.path.join(output_dir, out_name)
                    save_table_as_latex(table, output_file, log_path)
