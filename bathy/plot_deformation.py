#!/usr/bin/env python

from __future__ import division

import numpy
import matplotlib.pyplot as plt


def read_deformation_file(path):
    # Load data
    data = numpy.loadtxt(path)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2:]

    # Find number of points in each direction
    for i in xrange(1, x.shape[0]):
        if x[i] < x[i-1]:
            break
    nx = i - 1
    ny = x.shape[0] / nx

    import pdb; pdb.set_trace()

def plot_deformation_file(X, Y, Z):

    fig = plt.figure()

    return figure

if __name__ == "__main__":
    read_deformation_file("gapThi.xyzt")