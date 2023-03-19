import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 8
plt.rcParams['figure.autolayout'] = True

FLOPS = 8 * 1.8 * 1e9 * 16
f = lambda n : (2*n**3 / 3) / FLOPS
g = lambda n : (2*n**2) / FLOPS
fband = lambda n : 13.52*n**2 / FLOPS
gband = lambda n : 10.4*n**(3/2) / FLOPS

def plot(n, lu, solve, lu_band, solve_band, cholesky, solve_cholesky):
    
    _, ax = plt.subplots(nrows=1, ncols=2, figsize=(7, 9), dpi=350, num="Analyse Numérique: Devoir 2")
    sz = 8

    # LU decomposition function measurements
    ax[0].loglog(n, lu, label="lu", color='slateblue')
    ax[0].loglog(n, lu_band, label="lu_band", color='tomato')
    ax[0].loglog(n, cholesky, label="cholesky", color='lightsalmon')
    ax[0].set_xlabel("$m$", fontsize=sz)
    ax[0].set_ylabel("Temps d'exécution " + "$[s]$", fontsize=sz)
    ax[0].grid(True, which='both', lw=0.5)
    # Theorical results
    ax[0].loglog(n, [f(x) for x in n], label="lu théorique", color='seagreen')
    ax[0].loglog(n, [fband(x) for x in n], label="lu_band théorique", color='darkorange')
    ax[0].legend(loc='best', fontsize=sz-2)

    # Solve function measurements
    ax[1].loglog(n, solve, label="solve", color='slateblue')
    ax[1].loglog(n, solve_band, label="solve_band", color='tomato')
    ax[1].set_xlabel("$m$", fontsize=sz)
    
    ax[1].grid(True, which='both', lw=0.5, color='lightgrey')
    # Theorical results
    ax[1].loglog(n, [g(x) for x in n], label="solve théorique", color='seagreen')
    ax[1].loglog(n, [gband(x) for x in n], label="solve_band théorique", color='darkorange')

    ax[1].legend(loc='best', fontsize=sz-2)
    
    plt.show()


if __name__ == "__main__":

    # parsing command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_file', type=str, help='file.csv to plot')
    args = parser.parse_args()

    # parsing file
    data = np.loadtxt(args.csv_file, delimiter=',')

    # plot data
    plot(data[:, 0], data[:, 1], data[:, 2], data[:, 5], data[:, 6], data[:, 7], data[:, 8])