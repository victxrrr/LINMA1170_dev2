import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 8
plt.rcParams['figure.figsize'] = [5, 3]
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200

def plot(n, x, y, z):
    
    plt.plot(n, x)
    plt.plot(n, y)
    plt.plot(n, z)
    plt.show()


if __name__ == "__main__":

    # parsing command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('csv_file', type=str, help='file.csv to plot')
    args = parser.parse_args()

    # parsing file
    data = np.loadtxt(args.csv_file, delimiter=',')

    # plot data
    plot(data[:, 0], data[:, 1], data[:, 2], data[:, 3])