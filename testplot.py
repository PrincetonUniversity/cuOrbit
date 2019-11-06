#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')

def main(fn):
    if not os.path.exists(fn):
        raise ValueError("Input file {} does not exist!".format(fn))

    with open(fn, 'r') as f:
        lines = [l.strip() for l in f.readlines()]

    header = lines[:6]
    nz_lines = lines[6:]
    NNZ = len(nz_lines)

    #E, Pz, Mu, DE, DP = [np.array(l.split()) for l in header]
    # feature
    DE, DP, E, Pz, Mu  = [np.array(l.split()) for l in header[1:]]

    X = np.zeros(NNZ)
    Y = np.zeros(NNZ)
    W = np.zeros(NNZ)
    U = np.zeros(NNZ)
    V = np.zeros(NNZ)

    F = np.zeros(NNZ)

    # I'm going to guess we want to plot x and DE and y and DP based on Mario
    Z = np.zeros((DE.shape[0], DP.shape[0]))

    for n, line in enumerate(nz_lines):
        ls = line.split()
        val = np.float64(ls[-1])
        F[n] = val
        # these are indices
        e, pz, mu, de, dp = [int(item) for item in ls[:-1]]
        # we'll map to values, maybe wont use, 
        X[n] = E[e]
        Y[n] = Pz[pz]
        W[n] = Mu[mu]
        U[n] = DE[de]
        V[n] = DP[dp]

        # for now, just agg them up and lets take a peek
        Z[de,dp] += val

    plt.contour(DE, DP, Z, linewidths=2, cmap='RdGy')
    # We can still add a colorbar for the image, too.
    plt.colorbar(orientation='horizontal', shrink=0.8)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv)!=2:
        raise RuntimeError("Usage is ./{0} sparse_input_file.nz".format(sys.argv[0]))

    fn = sys.argv[1]
    main(fn)
