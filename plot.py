import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 8

def plot(X, Y):
    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(5, 3), dpi=200, num="Analyse Numérique: Devoir 2")
    #fig.suptitle("Mesure des performances des algorithmes définis dans lu.h")
    fig.tight_layout()

    # plot data
    ax[0].spy(X)
    ax[0].text(0.5, -0.1, "Taux de remplissage: {:.2f}%".format(np.count_nonzero(X)/X.size), horizontalalignment='center', verticalalignment='center', transform=ax[0].transAxes)
    ax[1].spy(Y)
    ax[1].text(0.5, -0.1, "Taux de remplissage: {:.2f}%".format(np.count_nonzero(Y)/Y.size), horizontalalignment='center', verticalalignment='center', transform=ax[1].transAxes)

    plt.show()


if __name__ == "__main__":

    # parsing command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_file1', type=str, help='file1.csv to plot')
    parser.add_argument('csv_file2', type=str, help='file2.csv to plot')
    args = parser.parse_args()

    # parsing file
    A = np.loadtxt(args.csv_file1, delimiter=',')
    B = np.loadtxt(args.csv_file2, delimiter=',')

    # plot data
    plot(A, B)