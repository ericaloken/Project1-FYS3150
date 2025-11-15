"""
Problem 10
Plotting the timings of the special and general algorithm
"""

import matplotlib.pyplot as plt
import numpy as np
import os 

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

n_list = [10, 100, 1000, 100000, 1000000]
runs = list(range(1, 101))                 # run algorithm 100 times per n for reliability

# ------------ plotting run nr. vs time ------------
fig, ax = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

# General algorithm
for n in n_list:
    t_general = np.loadtxt(f"../results/output_problem10_general_n{n}.txt")
    t_general = np.clip(t_general, 0, None)    # removing negative time
    exp = int(np.log10(n))
    ax[0].plot(runs, t_general*1e3, marker=".", markersize=4, label=fr"n = $10^{exp}$")

# Special algorithm
for n in n_list:
    t_special = np.loadtxt(f"../results/output_problem10_special_n{n}.txt")
    exp = int(np.log10(n))
    ax[1].plot(runs, t_special*1e3, marker=".", markersize=4, label=fr"n = $10^{exp}$")

ax[0].set_xlabel("Run nr.", fontsize=14)
ax[0].set_ylabel("Time [ms]", fontsize=14)
ax[0].set_title("General algorithm")
ax[0].legend()

ax[1].set_xlabel("Run nr.", fontsize=14)
ax[1].set_title("Special algorithm", fontsize=14)
ax[1].legend()

plt.suptitle("Time stability across 100 runs", fontsize=16)
plt.tight_layout()
plt.show()



# ------------ plotting n_steps vs time with errorbars------------

means_general = []
stds_general = []

for n in n_list:
    t_general = np.loadtxt(f"../results/output_problem10_general_n{n}.txt")
    t_general = np.clip(t_general, 0, None)  # removing negative time
    
    # calculating means and std for each n
    means_general.append(t_general.mean())
    stds_general.append(t_general.std())

means_special, stds_special = [], []
for n in n_list:
    t_special = np.loadtxt(f"../results/output_problem10_special_n{n}.txt")
    t_special = np.clip(t_special, 0, None)
    means_special.append(t_special.mean())
    stds_special.append(t_special.std())

plt.figure(figsize=(7,6))
plt.errorbar(n_list, means_general, yerr=stds_general, marker="o", capsize=3, label="General", color='darkmagenta')
plt.errorbar(n_list, means_special, yerr=stds_special, marker="s", capsize=3, label="Special", color='hotpink')

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$n_{\mathrm{steps}}$", fontsize=14)
plt.ylabel("Time [s]", fontsize=14)
plt.title(r"Computation time as function of $n_{steps}$", fontsize=16)
plt.legend()
plt.tight_layout()
plt.show()