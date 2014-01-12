#!/usr/bin/env python

import os
import sys
import subprocess
import glob
import re

import numpy
import scipy.interpolate as interp
import matplotlib.pyplot as plt

import clawpack.visclaw.colormaps as colormaps

def transform_coords(x, y, strike=22.0):
    r"""Transform coordinates from local fault coordinates to long-lat

    TODO: Add more parameters so that this function can be generalized to any
    fault setup.  Also want to make this an input to the *read_compsys_file*
    function.
    """
    # Rotate
    longitude = (x + 180.0) * numpy.cos(strike * numpy.pi / 180.0) \
              + (y + 60.00) * numpy.sin(strike * numpy.pi / 180.0)
    latitude  = -(x + 180.0) * numpy.sin(strike * numpy.pi / 180.0) \
              +  (y + 60.00) * numpy.cos(strike * numpy.pi / 180.0)

    # Transform to lat-long
    longitude = (longitude / 106.0) - 102.3
    latitude = (latitude / 111.13) + 16.2

    return longitude, latitude


def read_compsys_file(path, transform_coords, flip=False):
    r"""Read a COMPSYS generated deformation file

    Returns array of long lat tuples and Z array.  If X, Y arrays for plotting 
    are desired than use
    ```
    X = longlat_coords[0].reshape((num_cells[1], num_cells[0]))
    Y = longlat_coords[1].reshape((num_cells[1], num_cells[0]))
    Z = numpy.zeros((num_cells[1], num_cells[0], num_cells[2]))
    for n in xrange(num_cells[2]):
        Z[:,:,n] = z[:,n].reshape((num_cells[1], num_cells[0]))
    ```
    to get the correct shapes.

    input
    -----
     - *path* (path) - Path to COMPSYS deformation file.
     - *transform_coords (func) - Function that will perform the coordinate
       transformation from local fault coordinates (meters) to lat-long.
     - *flip* (bool) - Whether to flip the fault rupture propagation, default is
       *False*.

    output
    ------
     - *num_cells* (list) - Dimensions of the data in the x, y, and t dimensions
     - *longlat_coords* (list) - List of numpy.ndarrays containing the new
       transformed coordinates.
     - *z* (numpy.ndarray[:,:]) - Flattened array containing deformation data.
       The first indice is spatial and the second temporal, i.e. 
       *z.shape = (num_cells[0] * num_cells[1], num_cells[2])*

    """

    # Read in data from COMPSYS deformation file
    print "Loading data from %s" % path
    data = numpy.loadtxt(path)

    # Extract data arrays
    x = data[:,0]
    y = data[:,1]
    z = data[:,2:]

    # Flip fault rupture into other direction (E-W flip)
    if flip:
        x = numpy.flipud(x)

    # Transform (x,y) -> (long,lat)
    longlat_coords = transform_coords(x, y)

    # Find number of points in each direction of original data
    for i in xrange(1, data.shape[0]):
        if data[i,0] == data[0,0]:
            break
    num_cells = [i, int(data.shape[0] / i), data.shape[1] - 2]
    if num_cells[0] * num_cells[1] != data.shape[0]:
        raise ValueError("Error in calculating extents.")

    print "Done reading in data."

    # Convert from centimeters to meters
    z *= 0.01

    # Add zero-deformation time
    num_cells[2] += 1
    x = numpy.concatenate((x[0:num_cells[0]*num_cells[1]],x))
    y = numpy.concatenate((y[0:num_cells[0]*num_cells[1]],y))
    z_new = numpy.zeros((num_cells[0] * num_cells[1], num_cells[2]))
    z_new[:,1:] = z
    del z

    return num_cells, longlat_coords, z_new


def project_deformation(num_cells, longlat_coords, z, t_start=0.0, dt=5.0, 
                        scaling=1.0):
    r"""Project the input deformation onto a grid aligned with lat-long

    Given an array of coordinates *longlat_coords*, deformations *z* and the 
    number of cells in longitude, latitude and time project onto a grid that
    with coordinates aligned with longitude and latitude.  The time vector *t*
    that is returned is calculated via *t_start* and *dt*.

    input
    -----
     - *num_cells*
     - *longlat_coords*
     - *z*
     - *t_start*
     - *dt*
     - *scaling*

    output
    ------
     - *t* (numpy.ndarray[num_cells[2]]) - Array containing times which the
       deformation array *Z* has.
     - *X* (numpy.ndarray[:,:]) - Array containing x (longitude) coordinates as 
       would be given back by *numpy.meshgrid*.
     - *Y* (numpy.ndarray[:,:]) - Array containing y (latitude) coordinates as 
       would be given back by *numpy.meshgrid*.
     - *Z* (numpy.ndarray[:,:,:]) - Array containing deformation information at
       each time in *t*.

    """

    # Construct new grid
    new_num_cells = []
    new_num_cells.append(int(scaling * num_cells[1]))
    new_num_cells.append(int(scaling * num_cells[0]))
    x_new = numpy.linspace(numpy.min(longlat_coords[0]), numpy.max(longlat_coords[0]), new_num_cells[0])
    y_new = numpy.linspace(numpy.min(longlat_coords[1]), numpy.max(longlat_coords[1]), new_num_cells[1])
    X_new, Y_new = numpy.meshgrid(x_new, y_new)
    Z_new = numpy.zeros((new_num_cells[1], new_num_cells[0], num_cells[2]))
    t = numpy.arange(t_start, num_cells[2] * dt, dt)

    # Construct each interpolating function and evaluate at new grid
    for n in xrange(num_cells[2]):
        # Project onto rotated grid
        Z_new[:,:,n] = interp.griddata(longlat_coords, z[:,n], (X_new, Y_new), 
                                                method='linear', fill_value=0.0)

    return t, X_new, Y_new, Z_new


def write_deformation_to_file(t, X, Y, Z, output_file, topo_type=1):
    r"""Write out a dtopo file to *output_file*

    input
    -----
     - *t* (numpy.ndarray[:]) - Array containing time points, note that 
       *t.shape[0] == Z.shape[2]*.
     - *X* (numpy.ndarray[:,:]) - Array containing x-coodinates (longitude), 
       should be in the form given by *numpy.meshgrid*.
     - *Y* (numpy.ndarray[:,:]) - Array containing y-coordinates (latitude),
       should be in the form given by *numpy.meshgrid*.
     - *Z* (numpy.ndarray[:,:,:]) - Array containing deformation from original
       bathymetry.
     - *output_file* (path) - Path to the output file to written to.
     - *topo_type* (int) - Type of topography file to write out.  Default is 1.

    """

    # Temporary catch for non implemented topo_type input
    if topo_type != 1:
        raise NotImplementedError("Topography types 2 and 3 are not yet supported.")

    # Construct each interpolating function and evaluate at new grid
    try:
        outfile = open(output_file, 'w')

        if topo_type == 1:
            # Topography file with 4 columns, t, x, y, dz written from the upper
            # left corner of the region
            Y_flipped = numpy.flipud(Y)
            for n in xrange(t.shape[0]):
                Z_flipped = numpy.flipud(Z[:,:,n])
                for j in xrange(Y.shape[0]):
                    for i in xrange(X.shape[1]):
                        outfile.write("%s %s %s %s\n" % (t[n], X[j,i], Y_flipped[j,i], Z_flipped[j,i]))
    
        elif topo_type == 2 or topo_type == 3:
            raise NotImplementedError("Topography types 2 and 3 are not yet supported.")
        else:
            raise ValueError("Only topography types 1, 2, and 3 are supported.")

    except IOError as e:
        raise e
    finally:
        outfile.close()


def read_dtopo_file(path, topo_type=1):
    r""""""

    if topo_type != 1:
        raise ValueError("Topography type != 1 is not implemented.")

    # Load raw data
    data = numpy.loadtxt(path)

    # Parse data
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]

    # Initialize extents
    t0 = t[0]
    lower = [x[0], y[0]]
    upper = [None, None]
    num_cells = [0,0]

    # Count total x-values
    for row in xrange(1,data.shape[0]):
        if x[row] == x[0]:
            num_cells[0] = row
            break

    # Count total y-values
    for row in xrange(num_cells[0], data.shape[0]):
        if t[row] != t0:
            num_cells[1] = row / num_cells[0]
            num_times = data.shape[0] / row
            break

    # Check extents
    assert(t[0] != t[num_cells[0] * num_cells[1] + 1])
    # assert(t[0] == t[num_cells[0] * num_cells[1]])

    # Fill in rest of pertinent data
    t = data[::num_cells[0] * num_cells[1], 0]
    x = data[:num_cells[0], 1]
    y = data[:num_cells[0] * num_cells[1]:num_cells[0], 2]
    upper = [x[-1], y[-1]]
    X, Y = numpy.meshgrid(x, y)
    Z = numpy.empty( (num_times, num_cells[0], num_cells[1]) )

    for (n,time) in enumerate(t):
        Z[n,:,:] = data[num_cells[0] * num_cells[1] * n:
                        num_cells[0] * num_cells[1] * (n+1), 3].reshape(
                                (num_cells[0], num_cells[1]))

    return t, x, y, Z


def plot_deformations(plot, frames, make_movies=True):
    r""""""
    if frames is None or frames == 'all':
        frames = range(num_cells[2])

    if not os.path.exists(plot_output):
        os.makedirs(plot_output)

    # Choose which plots to make
    plot_objects = {}
    if plot == 'original' or plot == 'all':
        plot_objects['Original'] = [X_orig, Y_orig, Z_orig]
    if plot == 'projected' or plot == 'all':
        plot_objects['Projected'] = [X_new, Y_new, Z_new]
        if rotate:
            plot_objects['Projected Rotated'] = [X_new_rotated, 
                                                 Y_new_rotated, 
                                                 Z_new_rotated]

    # Find what colorbar limits would be best
    colorbar_limits = [numpy.infty, -numpy.infty]
    for frame in frames:
        colorbar_limits[0] = min(colorbar_limits[0], numpy.min(Z_orig[:,:,frame]))
        colorbar_limits[0] = min(colorbar_limits[0], numpy.min(Z_new[:,:,frame]))
        colorbar_limits[1] = max(colorbar_limits[1], numpy.max(Z_orig[:,:,frame]))
        colorbar_limits[1] = max(colorbar_limits[1], numpy.max(Z_new[:,:,frame]))

    # Equalize limits and choose type of colorbar
    if colorbar_limits[0] > 0.0 and colorbar_limits[1] > 0.0:
        # Use sequential colormap
        colorbar_limits[0] = colormaps.make_colormap({1.0:'r',0.0:'w'})
        cmap = plt.get_cmap('PuBu')
    elif colorbar_limits[0] <= 0.0 and colorbar_limits[1] <= 0.0:
        colorbar_limits[1] = 0.0
        cmap = colormaps.make_colormap({1.0:'w',0.0:'b'})
    elif colorbar_limits[0] <= 0.0 and colorbar_limits[1] >= 0.0:
        colorbar_limits[0] = -max(-colorbar_limits[0], colorbar_limits[1])
        colorbar_limits[1] = -colorbar_limits[0]
        cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})

    # Plot original and transformed next to each other
    for (n,frame) in enumerate(frames):
        for (plot_type, arrays) in plot_objects.iteritems():
            fig = plt.figure()
            axes = fig.add_subplot(111)

            im = axes.pcolormesh(arrays[0], arrays[1], arrays[2][:,:,frame],
                            vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                            cmap=cmap)


            axes.set_aspect('equal')
            axes.set_title('%s %s - Time t = %s seconds' %  \
                                    (plot_type, deformation_file, t[frame]))

            # Add colorbar
            colorbar = fig.colorbar(im, ax=axes)
            colorbar.set_label("Deformation (m)")

            # Save this figure
            plot_name = "%s_%s_%s.png" % (os.path.splitext(deformation_file)[0], 
                                          "_".join(plot_type.lower().split()), 
                                          str(frame).zfill(2))
            fig.savefig(os.path.join(plot_output, plot_name))
            plt.close(fig)

    if make_movies:
        print "Making movies..."
        for plot_type in plot_objects.iterkeys():
            loop = 0
            delay = 30
            name = "%s_%s" % (os.path.splitext(deformation_file)[0], 
                              "_".join(plot_type.lower().split()))
            png_files = [os.path.join(plot_output,png_file) 
                            for png_file in os.listdir(plot_output) 
                                if re.search(r'%s_[0-9].*\.png' % name, png_file)]
            cmd = 'convert -delay %s -loop %s %s %s.gif' % (delay, loop, 
                                                  " ".join(png_files), name)
            subprocess.Popen(cmd, shell=True)
        print "Done making movies."



if __name__ == "__main__":
    # Default parameters
    file_list = glob.glob('./gap*.xyzt')
    rotate = False
    frames = range(1,50,2)
    # frames = 'all'
    scaling = 5.0
    plot = 'all' #  == 'original' and 'projected'
    plot_output = "_plots"

    # Simple command line parsing
    if len(sys.argv) > 1:
        rotate = bool(sys.argv[1])
        if len(sys.argv) > 2:
            file_list = sys.argv[2:]

    # Loop over each file performing a transform, projection, plotting and writing
    for deformation_file in file_list:
        # =====================================
        #  Read and transform deformation file
        # =====================================
        num_cells, longlat_coords, z = read_compsys_file(deformation_file,
                                                         transform_coords)

        # TEMPORARY:  Remove last 2 invalid time points
        z = z[:,:-2]
        num_cells[2] -= 2

        # Put data into correct shapes for later
        X_orig = longlat_coords[0].reshape((num_cells[1], num_cells[0]))
        Y_orig = longlat_coords[1].reshape((num_cells[1], num_cells[0]))
        Z_orig = numpy.zeros((num_cells[1], num_cells[0], num_cells[2]))
        for n in xrange(num_cells[2]):
            Z_orig[:,:,n] = z[:,n].reshape((num_cells[1], num_cells[0]))

        # Project deformation into lat-long coordinates
        print "Projecting deformation..."
        t, X_new, Y_new, Z_new = project_deformation(num_cells, longlat_coords, 
                                                        z, scaling=scaling)
        print "done."
        
        # ===========================
        #  Write out projected dtopo
        # ===========================
        prefix = 'rot'
        out_path = os.getcwd()
        print "Writing out new deformation file %s_%s" % (prefix, os.path.basename(deformation_file))
        output_file = os.path.join(out_path, "%s_%s" % (prefix, os.path.basename(deformation_file)))
        write_deformation_to_file(t, X_new, Y_new, Z_new, output_file)
        print "Done writing out deformation file."

        # ====================================
        #  East-West Flip variation for dtopo
        # ====================================
        if rotate:
            num_cells, longlat_coords_flipped, z = \
                read_compsys_file(deformation_file, transform_coords, flip=True)

            t, X_new_rotated, Y_new_rotated, Z_new_rotated = \
                project_deformation(num_cells, longlat_coords_flipped, z,
                    scaling=scaling)

            # TEMPORARY:  Remove last 2 invalid time points
            z = z[:,:-2]
            num_cells[2] -= 2

            prefix = 'rot_west'
            out_path = os.getcwd()
            print "Writing out new deformation file %s_%s" % (prefix, os.path.basename(deformation_file))
            output_file = os.path.join(out_path, "%s_%s" % (prefix, os.path.basename(deformation_file)))
            write_deformation_to_file(t, X_new_rotated, Y_new_rotated, 
                                                     Z_new_rotated, output_file)
            print "Done writing out deformation file."

        # ======================================
        #  Make plots and movies of deformation
        # ======================================
        plot_deformations(plot, frames, make_movies=True)

        print