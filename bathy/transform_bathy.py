#!/usr/bin/env python

import os
import sys
import glob

import numpy
import scipy.interpolate as interp
import matplotlib.pyplot as plt

import clawpack.geoclaw.topotools as topotools

def transform_coords(x, y):
    # Rotate
    longitude = (x + 180.0) * numpy.cos(22.0 * numpy.pi / 180.0) \
              + (y + 60.00) * numpy.sin(22.0 * numpy.pi / 180.0)
    latitude  = -(x + 180.0) * numpy.sin(22.0 * numpy.pi / 180.0) \
              +  (y + 60.00) * numpy.cos(22.0 * numpy.pi / 180.0)

    # Transform to lat-long
    longitude = (longitude / 106.0) - 102.3
    latitude = (latitude / 111.13) + 16.2

    return longitude, latitude


def write_deformation_to_file(t, X, Y, Z, outfile, topo_type=1):
    r"""Write out a dtopo file"""

    if topo_type == 1:
        Y_flipped = numpy.flipud(Y)
        Z_flipped = numpy.flipud(Z)
        for j in xrange(Y.shape[0]):
            for i in xrange(X.shape[1]):
                outfile.write("%s %s %s %s\n" % (t, X[j,i], Y_flipped[j,i], Z_flipped[j,i]))
    elif topo_type == 2 or topo_type == 3:
        raise NotImplementedError("Topography types 2 and 3 are not yet supported.")
    else:
        raise ValueError("Only topography types 1, 2, and 3 are supported.")


def transform_deformation_file(path, out_path='./', t_start=0.0, dt=5.0, 
                                      scaling=1.0, debug=False, debug_frame=25):

    num_dim = 2

    # Read in data from original deformation file
    print "Loading data from %s" % path
    data = numpy.loadtxt(path)

    # Extract data arrays
    x = data[:,0]
    y = data[:,1]
    z = data[:,num_dim:]

    # Transform (x,y) -> (long,lat)
    longlat_coords = transform_coords(x, y)

    # Find number of points in each direction of original data
    for i in xrange(1, data.shape[0]):
        if data[i,0] == data[0,0]:
            break
    num_cells = [0, 0]
    num_cells[0] = i
    num_cells[1] = int(data.shape[0] / num_cells[0])
    if num_cells[0] * num_cells[1] != data.shape[0]:
        raise ValueError("Error in calculating extents.")

    if debug:
        # Rearrange input data to sanity
        X = longlat_coords[0].reshape((num_cells[1], num_cells[0]))
        Y = longlat_coords[1].reshape((num_cells[1], num_cells[0]))
        Z = numpy.zeros((num_cells[1], num_cells[0], data.shape[1] - 2))
        for n in xrange(data.shape[1] - num_dim):
            Z[:,:,n] = z[:,n].reshape((num_cells[1], num_cells[0]))

        # Plot original data in lat-long coordinates
        fig = plt.figure()
        axes = fig.add_subplot(121)
        axes.pcolor(X, Y, Z[:,:,debug_frame])
        axes.set_aspect('equal')

    print "Done reading in data."

    # Construct new grid
    new_num_cells = []
    for dim in xrange(2):
        new_num_cells.append(int(scaling * num_cells[dim]))
    x_new = numpy.linspace(numpy.min(longlat_coords[0]), numpy.max(longlat_coords[0]), new_num_cells[0])
    y_new = numpy.linspace(numpy.min(longlat_coords[1]), numpy.max(longlat_coords[1]), new_num_cells[1])
    X_new, Y_new = numpy.meshgrid(x_new, y_new)
    Z_new = numpy.zeros((new_num_cells[1], new_num_cells[0], data.shape[1] - 2))

    # Construct each interpolating function and evaluate at new grid
    print "Writing out new deformation file rot_%s" % os.path.basename(path)
    output_file = os.path.join(out_path, "rot_%s" % os.path.basename(path))
    try:
        file_handle = open(output_file, 'w')
        for n in xrange(data.shape[1] - 2):
            # Project onto rotated grid
            Z_new[:,:,n] = interp.griddata(longlat_coords, z[:,n], (X_new, Y_new), 
                                                    method='linear', fill_value=0.0)
            Z_new[:,:,n] *= 0.1

            # Write out new gridded file
            t = t_start + n * dt
            write_deformation_to_file(t, X_new, Y_new, Z_new[:,:,n], file_handle)
        file_handle.close()
    except IOError as e:
        raise e
    finally:
        file_handle.close()

    # Plot rotation if requested
    if debug:
        axes = fig.add_subplot(122)
        plot = axes.pcolor(X_new, Y_new, Z_new[:,:,debug_frame])
        axes.set_aspect('equal')
        fig.colorbar(plot)

        plt.show()
    print "Done writing out deformation file."


if __name__ == "__main__":
    file_list = glob.glob('./gap*.xyzt')
    if len(sys.argv) > 1:
        file_list = sys.argv[1:]
    for deformation_file in file_list:
        transform_deformation_file(deformation_file, scaling=4.0)