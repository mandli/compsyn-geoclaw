#!/usr/bin/env python

import numpy
import matplotlib.pyplot as plt
import mpl_toolkits.basemap.pyproj as pyproj

import bathy

def convert_UTM2longlat(file_path, zone=14, out_path=None):
    data = numpy.loadtxt(file_path)
    projector = pyproj.Proj(proj='utm',zone='14')
    longitude, latitude = projector(data[:,0], data[:,1], inverse=True)
    z = -data[:,2]

    if isinstance(out_path, basestring):
        bathy.write_unstructured(out_path, longitude, latitude, z)

    return longitude, latitude, z


if __name__ == "__main__":

    # base_bathy_file = "mexican_coast_pacific.tt3"
    base_bathy_file = "srtm_17_09.tt3"
    region_bathy = "acapulco30m.xyz"
    conv_region_bathy = "acapulco_converted_30m.xyz"
    new_bathy = "acapulco_projected_30m.tt3"

    # Convert from UTM and write out original data for use later in extraction
    print "Converting UTM to long-lat..."
    longitude, latitude, z = convert_UTM2longlat(region_bathy, 
                                                     out_path=conv_region_bathy)
    print "  ...done."

    # Smooth out region near shore
    # Use previous axes with data to show where smoothing will be
    print "Running filter..."
    center = (-99.9029,16.840)
    radius = bathy.meters2deg(200.0, center[1])
    filter_radius = bathy.meters2deg(30.0, center[1])
    for i in xrange(z.shape[0]):
        if radius >= numpy.sqrt((longitude[i] - center[0])**2 + (latitude[i] - center[1])**2):
            running_ave_depth = z[i]
            num_points_found = 1
            for j in xrange(z.shape[0]):
                if filter_radius >= numpy.sqrt((longitude[j] - longitude[i])**2 + (latitude[j] - latitude[i])**2):
                    running_ave_depth += z[j]
                    num_points_found += 1
            z[i] = running_ave_depth / num_points_found
    print "  ...done."

    # Find extent of acapulco data, buffered a bit
    dx = bathy.meters2deg(250.0, 0.5 * (latitude[-1] + latitude[0]))
    rect = [numpy.min(longitude) - dx, numpy.max(longitude) + dx, 
            numpy.min(latitude) - dx, numpy.max(latitude) + dx]

    # Extract and project original data
    Z, delta = bathy.extract(conv_region_bathy, base_bathy_file, 
                                    extent=rect, no_data_value=-9999)


    


    # Write out data
    bathy.write(new_bathy, Z, (rect[0], rect[2]), delta, 
                                        no_data_value=-9999)

    # Plot overlay of old bathymetry extent and SRTM data for new bathy
    axes = bathy.plot(new_bathy, region_extent=rect, limits=[-20,20])
    # axes.scatter(longitude, latitude, c=z, cmap=plt.get_cmap('terrain'), alpha=1.0)
    theta = numpy.linspace(0,2 * numpy.pi,360)
    axes.plot(radius * numpy.cos(theta) + center[0], 
              radius * numpy.sin(theta) + center[1],'y-',linewidth=10)

    plt.show()
