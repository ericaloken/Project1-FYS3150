import numpy as np
import matplotlib.pyplot as plt
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

# opening txt file from problem2.cpp and extracting data
data = np.loadtxt("../results/output_problem2.txt")
x = data[:,0]
u = data[:,1]

plt.figure(figsize=(12,8))
plt.plot(x, u, color='mediumvioletred')
plt.title(r"Analytic solution $u(x)$", fontsize=16)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$u(x)$", fontsize=14)
plt.grid()
plt.show()