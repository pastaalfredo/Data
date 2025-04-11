#!/usr/bin/env python3
import re
import json

log_file_path = {
    "log/walz.log": "eelec-Walz2018a",
    "log/pg.log":   "eelec-ACT-S",
    "log/gc.log":   "eelec-ACT-G",
    "log/jc.log":   "eelec-CT"
}

# Molecules + QM values
molecules = ['lithium#fluoride', 'lithium#chloride', 'lithium#bromide', 'sodium#fluoride', 'sodium#chloride', 
             'sodium#bromide', 'potassium#fluoride', 'potassium#chloride', 'potassium#bromide']

QM_values = [
    -850.53,
    -656.85,
    -609.54,
    -750.62,
    -607.90,
    -572.68,
    -719.80,
    -570.19,
    -538.38
]

molecule_data = {
    "lithium#fluoride": {"rmin": 1.564, "eelec-SAPT":  -850.53},
    "lithium#chloride": {"rmin": 2.021, "eelec-SAPT": -656.85},
    "lithium#bromide": {"rmin": 2.17,  "eelec-SAPT": -609.54},
    "sodium#fluoride": {"rmin": 1.926, "eelec-SAPT":  -750.62},
    "sodium#chloride": {"rmin": 2.361, "eelec-SAPT": -607.90},
    "sodium#bromide": {"rmin": 2.502, "eelec-SAPT": -572.68},
    "potassium#fluoride": {"rmin": 2.171, "eelec-SAPT":   -719.80},
    "potassium#chloride": {"rmin": 2.667, "eelec-SAPT":  -570.19},
    "potassium#bromide": {"rmin": 2.821, "eelec-SAPT":  -538.38}
}
results = {}

def extract_log_values(log_file, molecule):
    try:
        with open(log_file, 'r') as file:
            log_content = file.read()

        molecule_match = re.search(fr"Name: {re.escape(molecule)}.*", log_content, re.DOTALL)
        if molecule_match:
            coulomb_section = molecule_match.group(0)
            return coulomb_section
    except Exception as e:
        print(f"Error reading log file {log_file}: {e}")
    return None

for molecule, qm_value in zip(molecules, QM_values):
    molecule_name = molecule.replace('#', '#')
    
    data = molecule_data.get(molecule_name, {}).copy()
    data["QM"] = qm_value
    
    for log_file, field_name in log_file_path.items():
        log_values = extract_log_values(log_file, molecule)
        if log_values:
            coulomb_section = log_values.split("\n")
            
            for line in coulomb_section:
                if str(qm_value) in line:
                    line_values = line.split()
                    qm_index = line_values.index(str(qm_value)) if str(qm_value) in line_values else -1
                    if qm_index != -1 and qm_index + 1 < len(line_values):
                        try:
                            data[field_name] = float(line_values[qm_index + 1])
                        except ValueError:
                            data[field_name] = None

    results[molecule_name] = data

with open('sapt2_data.json', 'w') as json_file:
    json.dump(results, json_file, indent=4)

print("Data saved to sapt2_data.json")
