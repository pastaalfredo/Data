#!/usr/bin/env python

import numpy as np

#make grids for ESP calcs

delta=0.01 
x = 0
y = 0

z_values = np.arange(0, 4.51, delta)

with open("grid.dat", "w") as file:
    for z in z_values:
        file.write(f"{x} {y} {z:.2f}\n")
