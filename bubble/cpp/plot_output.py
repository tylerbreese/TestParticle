import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data
# Using pandas to read the CSV (equivalent to readtable)
A = pd.read_csv("output/sim_data_2026-05-03_12-13-07.csv")
sim = A.values # equivalent to table2array
B = pd.read_csv("output/sim_data_2026-05-03_12-13-07.csv")
split = B.values

# Extract and scale columns (Note: Python uses 0-based indexing)
x0 = sim[:, 1] / 1e5
xf = sim[:, 2] / 1e5
xs = split[:, 3] / 1e5
En0 = sim[:, 3]
En1 = sim[:, 4]
# Ens = split[:, 5]

# Histogram Figure
plt.figure()
# MATLAB's histcounts(data, 25) returns 25 bins (26 edges)
n0, bins0 = np.histogram(En0, bins=25)
n1, bins1 = np.histogram(En1, bins=25)

# Normalizing and plotting (scatter uses the left edge of bins here)
# Added 0.0 to match MATLAB's manual padding if needed for alignment
plt.scatter(bins0[:-1], n0 / np.sum(n0), label='En0', facecolors='none', edgecolors='b')
plt.scatter(bins1[:-1], n1 / np.sum(n1), label='En1', facecolors='none', edgecolors='r')

plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.title("Normalized Histograms")

# Scatter Plot Figure
plt.figure()
plt.scatter(x0, En0, label='Initial', alpha=0.6)
plt.scatter(xf, En1, label='Final', alpha=0.6)
# plt.scatter(xs, Ens, label='Split', alpha=0.6)

plt.yscale('log')
plt.grid(True)
plt.legend()
plt.xlabel('Position (scaled)')
plt.ylabel('Energy')

plt.show()
