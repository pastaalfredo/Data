#!/usr/bin/env python
import os

a="SC-Water-Ion"
b="PC+GV"

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



def save_table_as_latex(table, output_file):
    with open(output_file, 'w') as file:
        file.write("\\begin{table}[ht]\n")
        file.write(f"\\caption{{Parameters in {b} trained on the electrostatic energy.}}\n")
        file.write("\\begin{tabular}{|c|c|c|}\n")
        file.write("\\hline\n")
        file.write("CLASS & NAME & BEST (Train) \\\\ \n")
        file.write("\\hline\n")
        for row in table:
            file.write(" & ".join(row) + " \\\\ \n")
        file.write("\\hline\n")
        file.write("\\end{tabular}\n")
        file.write("\\end{table}")

log_file_path = os.path.join(os.getcwd(), 'train-Inter6.log')
output_file_path = os.path.join(os.getcwd(), f'output_table_{a}_{b}.tex')

table = extract_table_from_log(log_file_path)

save_table_as_latex(table, output_file_path)
