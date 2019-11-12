#!/usr/bin/env python3
# This file is in the public domain.

import numpy as np
import matplotlib.pyplot as plt

# http://matplotlib.org/examples/color/colormaps_reference.html
cmaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis']

gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

def gen_table(cmap_name):
    # matplotlib.colors.ListedColormap
    cmap=plt.get_cmap(cmap_name)
    #print("// Size: {}".format(cmap.N))
    #print("// Size: {}".format(len(cmap.colors)))
    print("float {}Colormap[{}][3] = {{".format(cmap_name, cmap.N))
    for colors in cmap.colors:
        print("\t{{{}, {}, {}}},".format(colors[0], colors[1], colors[2]))
    print("};\n")

for cmap_name in cmaps:
    gen_table(cmap_name)
