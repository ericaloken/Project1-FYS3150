"""
Problem 7b:
Plotting x, v and u for different values of n
"""

import numpy as np
import matplotlib.pyplot as plt
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

# extracting and plotting data from analytical solution
data_exact = np.loadtxt("../results/output_problem2.txt")
x_exact = data_exact[:,0]
u_exact = data_exact[:,1]
plt.figure(figsize=(12,10))
plt.plot(x_exact, u_exact, label="analytic u(x)")


# plotting data from general algorithm
n_list = [10, 100, 1000, 100000, 10000000]

for n in n_list:
    data_algorithm = np.loadtxt(f"../results/output_problem9_n{n}.txt")
    x = data_algorithm[:,0]
    v = data_algorithm[:,1]
    plt.plot(x, v, label=f"n = {n}")

plt.title("Analytic solution u(x) and v(x) from special algorithm", fontsize=16)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$u(x)$", fontsize=14)
plt.legend()
plt.grid()
plt.show()