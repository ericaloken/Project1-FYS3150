"""
Problem 8a
Making a plot of the absoute error log10(Δ) = log10(|u_i-v_i|)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

# analytic solution
def u_exact(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10*x)

# calculating and plotting absolute error
n_list = [10, 100, 1000, 100000, 10000000]
plt.figure(figsize=(12,8))

for n in n_list:
    data = np.loadtxt(f"../results/output_problem7_n{n}.txt")
    x = data[:,0]
    v = data[:,1]
    exp = int(np.log10(n))

    # computing analytical solution
    u = u_exact(x)

    # calculation absolute error
    delta = np.abs(u - v)
    
    # plotting absolute error (excluding boundary points with zero error)
    plt.plot(x[1:-1], np.log10(delta[1:-1]), label=fr"n = $10^{exp}$")
    
plt.title(r"Absolute error $\log_{10}(\Delta) = \log_{10}(|u_i - v_i|)$", fontsize=16)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$\log_{10}\!\left(\left|\frac{u_i - v_i}{u_i}\right|\right)$", fontsize=14)
plt.legend()
plt.grid()
plt.show()

"""
Problem 8b
Plotting the relative error log10(ε) = log10(|(u_i-v_i) / u_i|)
"""

plt.figure(figsize=(12,8))

for n in n_list:
    data = np.loadtxt(f"../results/output_problem7_n{n}.txt")
    x = data[:,0]
    v = data[:,1]
    exp = int(np.log10(n))

    # computing analytical solution
    u = u_exact(x)

    # calculation relative error
    epsilon = np.zeros_like(u)
    for i in range(n):
        if abs(u[i]) > 1e-15:
            epsilon[i] = np.abs((u[i] - v[i]) / u[i])

        else:
            epsilon[i] = 0.0


    # printing maximum relative error
    print(f'n = {n}         max err_rel = {np.max(epsilon)}')
    
    # plotting relative error (excluding boundary points with zero error)
    plt.plot(x[1:-1], np.log10(epsilon[1:-1]), label=fr"n = $10^{exp}$")

plt.title(r"Relative error $\log_{10}(\epsilon) = \log_{10}\left(\left|\frac{u - v}{u}\right|\right)$", fontsize=16)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$\log_{10}\left(\left|\frac{u - v}{u}\right| \right)$", fontsize=14)
plt.legend()
plt.grid()
plt.show()