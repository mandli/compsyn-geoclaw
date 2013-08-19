#!/usr/bin/env python

import os
import sys
import glob

import numpy
import scipy.interpolate as interp
import matplotlib.pyplot as plt

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


def read_compsys_file(path, transform_coords):
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

    return num_cells, longlat_coords, z


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


def plot_deformation_comparison(path, transform_coords, frames='all'):
    r"""Plot deformation comparisons between original COMPSYS file and rotated

    input
    -----
     - *path* (path) - Path to COMPSYS deformation file
     - *frames* (list or string) - List of time frames to plot or keyword 'all'
       to plot all available frames.  Default is 'all'.

    output
    ------
     - (matplotlib.pyplot.figure) - Figure containing all plots.
    """

    # Read in data from original deformation file
    num_cells, longlat_coords, z = read_compsys_file(path, transform_coords)

    # Put data into correct shapes
    X = longlat_coords[0].reshape((num_cells[1], num_cells[0]))
    Y = longlat_coords[1].reshape((num_cells[1], num_cells[0]))
    Z = numpy.zeros((num_cells[1], num_cells[0], num_cells[2]))
    for n in xrange(num_cells[2]):
        Z[:,:,n] = z[:,n].reshape((num_cells[1], num_cells[0]))

    # Transform to new arrays
    t, X_new, Y_new, Z_new =        \
                        project_deformation(num_cells, longlat_coords, z,
                                                    scaling=4.0)

    # Find what colorbar limits would be best
    colorbar_limits = [numpy.infty, -numpy.infty]
    for frame in frames:
        colorbar_limits[0] = min(colorbar_limits[0], numpy.min(Z[:,:,frame]))
        colorbar_limits[0] = min(colorbar_limits[0], numpy.min(Z_new[:,:,frame]))
        colorbar_limits[1] = max(colorbar_limits[1], numpy.max(Z[:,:,frame]))
        colorbar_limits[1] = max(colorbar_limits[1], numpy.max(Z_new[:,:,frame]))

    # Equalize limits?


    # Create frames list if not given or 'all' is used
    if frames is None or frames == 'all':
        frames = range(num_cells[2])

    # Plot original and transformed next to each other
    fig, axes = plt.subplots(2, len(frames))
    fig.suptitle('Deformation for %s' % os.path.basename(path))
    for (n,frame) in enumerate(frames):
        # Plot original data in lat-long coordinates
        axes[0,n].pcolormesh(X, Y, Z[:,:,frame], vmin=colorbar_limits[0], 
                                                      vmax=colorbar_limits[1])
        axes[0,n].set_aspect('equal')
        axes[0,n].set_title('Time t = %s seconds' % t[frame])

        im = axes[1,n].pcolormesh(X_new, Y_new, Z_new[:,:,frame], 
                                                   vmin=colorbar_limits[0], 
                                                   vmax=colorbar_limits[1])
        axes[1,n].set_aspect('equal')

    # Add colorbar
    colorbar = fig.colorbar(im, ax=list(axes.flatten()))
    colorbar.set_label("Deformation (m)")
    # colorbar.set_clim(colorbar_limits)

    return fig


if __name__ == "__main__":
    # Simple command line parsing
    command = 'transform'
    file_list = 'all'
    if len(sys.argv) > 1:
        command = sys.argv[1]
        if len(sys.argv) > 2:
            file_list = sys.argv[2:]
    if file_list == 'all':
        file_list = glob.glob('./gap*.xyzt')

    # Execute approrpiate command
    if command == 'plot':
        for deformation_file in file_list:
            fig = plot_deformation_comparison(deformation_file, transform_coords, 
                                          frames=[0,25,49])
        plt.show()
    elif command == 'transform':
        for deformation_file in file_list:
            # Read and transform deformation file
            num_cells, longlat_coords, z = read_compsys_file(deformation_file,
                                                             transform_coords)

            # Project deformation into lat-long coordinates
            t, X, Y, Z = project_deformation(num_cells, longlat_coords, z,
                                             scaling=4.0)
            
            # Write out new grid
            prefix = 'rot'
            out_path = os.getcwd()
            print "Writing out new deformation file %s_%s" % (prefix, os.path.basename(deformation_file))
            output_file = os.path.join(out_path, "%s_%s" % (prefix, os.path.basename(deformation_file)))
            write_deformation_to_file(t, X, Y, Z, output_file)
            print "Done writing out deformation file."
    else:
        raise ValueError("Invalid command, use 'plot' or 'transform'")