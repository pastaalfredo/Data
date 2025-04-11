#!/usr/bin/env python
import os
import re
import json

log_file_path = "train-Inter6.log"

# Molecules + QM values
molecules = ['formate#lithium', 'formate#sodium', 'formate#potassium', 'formate#water',
             'acetate#lithium', 'acetate#sodium', 'acetate#potassium', 'acetate#water',
             'methylammonium#fluoride', 'methylammonium#chloride', 'methylammonium#bromide', 'methylammonium#water',
             'ethylammonium#fluoride', 'ethylammonium#chloride', 'ethylammonium#bromide', 'ethylammonium#water']

QM_values=[
        -734.50,
        -673.98,
        -604.73,
        -136.61,
        -624.68,
        -646.24,
        -574.45,
        -142.25,
        -458.91,
        -452.92,
        -434.94,
        -100.75,
        -440.31,
        -399.66,
        -385.98,
        -99.48]

QM_values = [f'{value:.2f}' for value in QM_values]

coulomb_value_pattern = re.compile(r"^\s*(?:\d+\s+)?(?:-?\d+\.\d+\s+){3}(-?\d+\.\d+)\s+", re.MULTILINE)

results = {}

if os.path.exists('results.json'):
    with open('results.json', 'r') as json_file:
        results = json.load(json_file)

with open(log_file_path, 'r') as file:
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

with open('results.json', 'w') as json_file:
    json.dump(results, json_file, indent=4)

print("Data saved to results.json")
