#!/usr/bin/env python

import numpy
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import mpl_toolkits.basemap.pyproj as pyproj

import clawpack.visclaw.colormaps as colormaps

import bathy

def UTM2longlat(x,y,zone='14'):
    r"""
    Convert from UTM to longitude-latitude coordinates
    """
    projector = pyproj.Proj(proj='utm', zone=zone)
    longitude, latitude = projector(x, y, inverse=True)

    return longitude, latitude


def project(x, y, z, X_fill, Y_fill, Z_fill, extent=None,
                                                    method='nearest', 
                                                    delta_limit=20.0, 
                                                    no_data_value=-9999,
                                                    buffer_length=100.0,
                                                    proximity_radius=100.0):
        
    # Calculate new grid coordinates
    if extent is None:
        buffer_degrees = bathy.meters2deg(buffer_length, 0.5 * (y[0] + y[-1]))
        extent = [numpy.min(x) - buffer_degrees, numpy.max(x) + buffer_degrees, 
                  numpy.min(y) - buffer_degrees, numpy.max(y) + buffer_degrees]
    delta = max( min(numpy.min(numpy.abs(x[1:] - x[:-1])), 
                     numpy.min(numpy.abs(y[1:] - y[:-1])) ), 
                 bathy.meters2deg(delta_limit, 0.5 * (y[0] + y[-1])) )
    N = ( numpy.ceil((extent[1] - extent[0]) / delta),
          numpy.ceil((extent[3] - extent[2]) / delta) )
    assert numpy.all(N[:] < numpy.ones((2)) * resolution_limit), \
           ValueError("Calculated resolution too high, N=%s!" % str(N))
    X, Y = numpy.meshgrid( numpy.linspace(extent[0], extent[1], N[0]),
                           numpy.linspace(extent[2], extent[3], N[1]) )
    Z = numpy.empty(N)

    # Create extent mask
    extent_mask = extent[0] > X_fill
    extent_mask = numpy.logical_or(extent_mask,extent[1] < X_fill)
    extent_mask = numpy.logical_or(extent_mask,extent[2] > Y_fill)
    extent_mask = numpy.logical_or(extent_mask,extent[3] < Y_fill)
    
    # Create fill no-data value mask
    no_data_mask = numpy.logical_or(extent_mask, Z_fill == no_data_value)

    all_mask = numpy.logical_or(extent_mask, no_data_mask)

    # Create proximity mask
    proximity_mask = numpy.ndarray(X_fill.shape, dtype=bool)
    proximity_radius_deg = bathy.meters2deg(proximity_radius, 0.5 * (y[0] + y[-1]))
    
    indices = (~all_mask).nonzero()
    count = 0
    for n in xrange(indices[0].shape[0]):
        i = indices[0][n]
        j = indices[1][n]
        all_mask[i,j] = numpy.any(numpy.sqrt( (x - X_fill[i,j])**2 + (y - Y_fill[i,j])**2 ) < proximity_radius_deg)
        count += int(all_mask[i,j])
    print "Masked %s/%s points." % (count, indices[0].shape[0])

    X_fill_masked = numpy.ma.masked_where(all_mask, X_fill)
    Y_fill_masked = numpy.ma.masked_where(all_mask, Y_fill)
    Z_fill_masked = numpy.ma.masked_where(all_mask, Z_fill)    

    # Stick both the input data and fill data into arrays
    fill_points = numpy.column_stack((X_fill_masked.compressed(),
                                   Y_fill_masked.compressed()))
    points = numpy.concatenate((numpy.array([x,y]).transpose(), fill_points))
    values = numpy.concatenate((z, Z_fill_masked.compressed()))

    # Use nearest-neighbor interpolation
    Z = interpolate.griddata(points, values, (X,Y), method=method)

    return X, Y, Z, delta, extent



def in_poly(X, Y, polygon):
    r"""
    Masks points (x,y) that are not in the specified polygon

    Uses simple ray casting algorithm for speed so beware of corner cases!

    Input
    -----
     - *X* (numpy.ndarray) Coordinates in x direction in a meshgrid type of
       configuration.
     - *Y* (numpy.ndarray) Coordinates in y direction in a meshgrid type of
       configuration.
     - *polygon* (list) List of points that comprise the polygon.  Note that
       order of the points will effect if this works (positive versus negative
       winding order).  Points should be in counter-clockwise arrangement.

    Returns
    -------
     - *X_mask* (numpy.ma.MaskedArray) Masked array of X coordinates where those
       points outside of the polygon have been masked.
     - *Y* (numpy.ndarray) Coordinates in y direction in a meshgrid type of
       configuration.
    """

    TOLERANCE = 1e-6

    # Flatten the input arrays to make this a bit easier
    x = X.flatten()
    y = Y.flatten()

    # Construct edges
    edges = []
    for edge in xrange(len(polygon) - 1):
        edges.append([polygon[edge], polygon[edge+1]])
    edges.append([polygon[-1], polygon[0]])

    # Check for intersections
    num_intersections = numpy.zeros(x.shape[0])

    for edge in edges:
        # Check for a vertical line
        if numpy.abs(edge[0][0] - edge[1][0]) < TOLERANCE:
            x_intersect = edge[0][0]        
        else:
            edge_slope = (edge[0][1] - edge[1][1]) / (edge[0][0] - edge[1][0])
            x_intersect = (y - edge[0][1]) / edge_slope + edge[0][0]

        num_intersections += (min(edge[0][1], edge[1][1]) <= y) * \
                             (max(edge[0][1], edge[1][1]) >= y) * \
                             (x_intersect <= x)
                             

    # General intersection of two lines
    intersect = (numpy.mod(num_intersections, numpy.ones(x.shape) * 2) != 1)

    # Return masked arrays that are reshaped back to the input shapes
    return numpy.ma.masked_where(intersect, x, copy=False).reshape(X.shape), \
           numpy.ma.masked_where(intersect, y, copy=False).reshape(Y.shape)
           

if __name__ == "__main__":

    # Parameters
    # Output paths
    orig_bathy = "acapulco30m.xyz"
    out_bathy = "acapulco_projected_30m.tt3"
    fill_bathy = 'srtm_17_09.tt3'

    # Algorithm parameters
    delta_limit = 20.0
    resolution_limit = 1000
    no_data_value = -9999
    bathy_limits = [-5.0,0.0]

    # Read in transform data
    orig_data = numpy.loadtxt(orig_bathy)
    x, y = UTM2longlat(orig_data[:,0], orig_data[:,1])
    z = -orig_data[:,2]

    # Read in fill data
    X_fill, Y_fill, Z_fill = bathy.read(fill_bathy)

    # Project data into new grid and write it out
    X, Y, Z, delta, extent = project(x, y, z, X_fill, Y_fill, Z_fill)
    bathy.write(out_bathy, Z, (extent[0], extent[2]), delta)

    # ===============
    #  Plot new data
    # ===============
    plt.ticklabel_format(format="plain", useOffset=False)
    
    fig = plt.figure()
    axes = fig.add_subplot(1,1,1)
    mean_lat = 0.5 * (extent[3] - extent[2])
    axes.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))    
    axes.set_autoscale_on(True)
    axes.set_title("Acapulco Bathymetry")
    axes.set_xlabel("Longitude (degrees)")
    axes.set_ylabel("Latitude (degrees)")

    # Create colorbar
    cmap = plt.get_cmap("terrain")
    bathy_limits = (numpy.min(Z), numpy.max(Z))
    color_norm = colors.Normalize(bathy_limits[0], bathy_limits[1], clip=True)
    
    axes.contourf(X_fill, Y_fill, Z_fill, vmin=bathy_limits[0], 
                                          vmax=bathy_limits[1],
                                          cmap=cmap, 
                                          norm=color_norm)
    plot = axes.contourf(X, Y, Z, vmin=bathy_limits[0], vmax=bathy_limits[1],
                                  cmap=cmap, norm=color_norm)

    axes.set_xlim(extent[0:2])
    axes.set_ylim(extent[2:])

    # Add colorbar
    cbar = fig.colorbar(plot, ax=axes)
    cbar.set_label("Depth (m)")

    plt.show()

    # # Plot overlay of old bathymetry extent and SRTM data for new bathy
    # axes = bathy.plot(new_bathy, region_extent=rect, limits=[-20,20], coastlines=False)
    # # axes.scatter(longitude, latitude, c=z, cmap=plt.get_cmap('terrain'), alpha=1.0)

    # X_fill, Y_fill, Z_fill = bathy.read(base_bathy_file)
    # extent_mask = rect[0] > X_fill
    # extent_mask = numpy.logical_or(extent_mask,rect[1] < X_fill)
    # extent_mask = numpy.logical_or(extent_mask,rect[2] > Y_fill)
    # extent_mask = numpy.logical_or(extent_mask,rect[3] < Y_fill)
    # all_mask = numpy.logical_or(extent_mask, Z_fill == -9999)

    # X_fill_mask = numpy.ma.masked_where(all_mask,X_fill)
    # Y_fill_mask = numpy.ma.masked_where(all_mask,Y_fill)
    # Z_fill_mask = numpy.ma.masked_where(all_mask,Z_fill,-9999)

    # axes.pcolor(X_fill_mask, Y_fill_mask, Z_fill_mask)

    # theta = numpy.linspace(0,2 * numpy.pi,360)
    # axes.plot(radius * numpy.cos(theta) + center[0], 
    #           radius * numpy.sin(theta) + center[1],'y:')

    # plt.show()

    # Filter before projection
    # for i in xrange(z.shape[0]):
    #     if radius >= numpy.sqrt((longitude[i] - center[0])**2 + (latitude[i] - center[1])**2):
    #         running_ave_depth = z[i]
    #         num_points_found = 1
    #         for j in xrange(z.shape[0]):
    #             if filter_radius >= numpy.sqrt((longitude[j] - longitude[i])**2 + (latitude[j] - latitude[i])**2):
    #                 running_ave_depth += z[j]
    #                 num_points_found += 1
    #         # z[i] = running_ave_depth / num_points_found
    #         z[i] = -20.0
    # print "  ...done."


    # Smooth out region near shore
    # Use previous axes with data to show where smoothing will be
    # print "Running filter..."
    # center = (-99.9029,16.840)
    # radius = bathy.meters2deg(200.0, center[1])
    # filter_radius = 2
    # import pdb; pdb.set_trace()
    # for (i,longitude) in enumerate(X[0,:]):
    #     for (j,latitude) in enumerate(Y[:,0]):
    #         if radius >= numpy.sqrt((longitude - center[0])**2 + (latitude - center[1])**2):
    #             running_ave_depth = Z[i,j]
    #             num_points = 1
    #             for n in xrange(1,filter_radius):
    #                 for m in xrange(1,filter_radius):
    #                     running_ave_depth += Z[i - n, j - m]
    #                     running_ave_depth += Z[i + n, j - m]
    #                     running_ave_depth += Z[i - n, j + m]
    #                     running_ave_depth += Z[i + n, j + m]
    #                     num_points += 4

    #             print "Changing Z[%s,%s] from %s to %s." % (i,j,Z[i,j],running_ave_depth / num_points)
    #             Z[i,j] = running_ave_depth / num_points

    # print "   ...done."