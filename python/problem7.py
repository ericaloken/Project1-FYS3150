"""
Problem 7b:
Plotting x, v(x) and u(x) for different values of n
"""

import numpy as np
import matplotlib.pyplot as plt
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

# extracting data from analytical solution
data_exact = np.loadtxt("../results/output_problem2.txt")
x_exact = data_exact[:,0]
u_exact = data_exact[:,1]

# plotting
n_list = [10, 100, 1000, 100000, 10000000]
plt.figure(figsize=(12,10))

for n in n_list:
    data_algorithm = np.loadtxt(f"../results/output_problem7_n{n}.txt")
    x = data_algorithm[:,0]
    v = data_algorithm[:,1]
    exp = int(np.log10(n))
    plt.plot(x, v, label=fr"n = $10^{exp}$")

plt.plot(x_exact, u_exact, color='mediumvioletred', linestyle='--', label=r"analytic $u(x)$")
plt.title(r"Analytic solution $u(x)$ and numerical solution $v(x)$ using general algorithm", fontsize=16)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$u(x)$", fontsize=14)
plt.legend()
plt.grid()
plt.show()