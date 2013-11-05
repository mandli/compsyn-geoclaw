#!/usr/bin/env python

import numpy

import matplotlib.pyplot as plt
import mpl_toolkits.basemap.pyproj as pyproj

import extract_bathy
import bathy

def convert_UTM2longlat(file_path, zone=14, out_path=None):
    data = numpy.loadtxt(file_path)
    projector = pyproj.Proj(proj='utm',zone='14')
    longitude, latitude = projector(data[:,0], data[:,1], inverse=True)
    z = -data[:,2]

    if isinstance(out_path, basestring):
        write_bathy(out_path, longitude, latitude, z)

    return longitude, latitude, z


def write_bathy(file_path, longitude, latitude, z, topotype=1):
    bathy_file = open(file_path, 'w')
    for (i, depth) in enumerate(z):
        bathy_file.write("%s %s %s\n" % (longitude[i], latitude[i], depth))

    bathy_file.close()

if __name__ == "__main__":

    # base_bathy_file = "mexican_coast_pacific.tt3"
    base_bathy_file = "srtm_17_09.tt3"
    region_bathy = "acapulco30m.xyz"
    conv_region_bathy = "acapulco_converted_30m.xyz"
    new_bathy = "acapulco_projected_30m.tt3"

    # Convert and write out UTM data
    acapulco = convert_UTM2longlat(region_bathy, out_path=conv_region_bathy)

    # Find extent of acapulco data
    rect = [numpy.min(acapulco[0]), numpy.max(acapulco[0]), 
            numpy.min(acapulco[1]), numpy.max(acapulco[1])]

    # Extract and project data
    Z, delta = extract_bathy.extract(conv_region_bathy, base_bathy_file, 
                                        extent=rect, no_data_value=-9999)
    extract_bathy.write_bathy(new_bathy, Z, (rect[0], rect[2]), delta, 
                                            no_data_value=-9999)
    # write_bathy(new_bathy, )
    # Z, delta = extract_bathy.extract(new_bathy, "mexican_coast_pacific.tt3",
                                        # extent=rect)

    # Plot new data
    bathy.plot(new_bathy, coastlines=True)
    plt.show()