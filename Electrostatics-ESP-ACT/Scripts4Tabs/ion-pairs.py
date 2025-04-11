#!/usr/bin/env python
import os
import re
import json

log_files = ["log/pg.log", "log/jc.log", "log/gc.log", "log/walz.log"]

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

QM_values = [f'{value:.2f}' for value in QM_values]
coulomb_value_pattern = re.compile(r"^\s*(?:\d+\s+)?(?:-?\d+\.\d+\s+){3}(-?\d+\.\d+)\s+", re.MULTILINE)

for log_file in log_files:
    results = {}

    if os.path.exists(log_file):
        with open(log_file, 'r') as file:
            log_content = file.read()

        for molecule, qm_value in zip(molecules, QM_values):
            molecule_match = re.search(fr"Name: {re.escape(molecule)}.*", log_content, re.DOTALL)
            if molecule_match:
                coulomb_section = molecule_match.group(0)
                qm_index = coulomb_section.split().index(str(qm_value))
                act_value = float(coulomb_section.split()[qm_index + 1])
            else:
                act_value = None

            results[molecule] = {'QM': qm_value, 'ACT': act_value}

        output_filename = f"results-ions-{os.path.basename(log_file)}.json"
        with open(output_filename, 'w') as json_file:
            json.dump(results, json_file, indent=4)

        print(f"Data saved to {output_filename}")
    else:
        print(f"Log file {log_file} not found.")
