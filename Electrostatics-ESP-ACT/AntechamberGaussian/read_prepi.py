#!/usr/bin/env python3

input_files = [
    'acetate.prepi',
    'ammonium.prepi',
    'formate.prepi',
    'ethylammonium.prepi',
    'methylammonium.prepi'
]
output_file = 'combined_charges.txt'

def extract_columns(input_file, outfile, skip_lines=7):
    compound_name = input_file.split('.')[0]
    outfile.write(f"{compound_name}\n")
    with open(input_file, 'r') as infile:
        for _ in range(skip_lines):
            next(infile)
        
        for line in infile:
            if line.strip() == 'LOOP':
                break
            
            parts = line.split()
            if len(parts) >= 2 and parts[1] != 'DUMM':  
                second_column = parts[1]  
                last_column = parts[-1] 
                outfile.write(f"{second_column} {last_column}\n")
    outfile.write("\n")  


with open(output_file, 'w') as outfile:
    for input_file in input_files:
        extract_columns(input_file, outfile)
