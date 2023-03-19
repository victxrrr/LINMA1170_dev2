import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 6
plt.rcParams['figure.autolayout'] = True
plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = False

def plot(X, Y, Z, W, V, U):
    _, ax = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(10, 6), dpi=400, num="Analyse Numérique: Devoir 2")
    
    # plot data
    ax[0][0].spy(X)
    ax[0][0].set_ylabel("$K$")
    ax[1][0].spy(Y)
    ax[1][0].set_ylabel("$K=L*U$")
    ax[0][1].spy(Z)
    ax[0][1].set_ylabel(r'$K_{\mathrm{perm}}$')
    ax[1][1].spy(W)
    ax[1][1].set_ylabel(r'$K_{\mathrm{perm}} = L*U$')
    ax[0][2].spy(V)
    ax[0][2].set_ylabel("$K_{\mathrm{sym}}$")
    ax[1][2].spy(U)
    ax[1][2].set_ylabel("$K_{\mathrm{sym}} = L*L^T$")

    plt.show()

if __name__ == "__main__":

    # parsing command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_file1', type=str, help='file1.csv to plot')
    parser.add_argument('csv_file2', type=str, help='file2.csv to plot')
    parser.add_argument('csv_file3', type=str, help='file3.csv to plot')
    parser.add_argument('csv_file4', type=str, help='file4.csv to plot')
    parser.add_argument('csv_file5', type=str, help='file5.csv to plot')
    parser.add_argument('csv_file6', type=str, help='file6.csv to plot')
    args = parser.parse_args()

    # parsing file
    A = np.loadtxt(args.csv_file1, delimiter=',')
    B = np.loadtxt(args.csv_file2, delimiter=',')
    C = np.loadtxt(args.csv_file3, delimiter=',')
    D = np.loadtxt(args.csv_file4, delimiter=',')
    E = np.loadtxt(args.csv_file5, delimiter=',')
    F = np.loadtxt(args.csv_file6, delimiter=',')

    # plot data
    plot(A, B, C, D, E, F)