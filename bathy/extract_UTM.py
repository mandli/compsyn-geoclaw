#!/usr/bin/env python

import sys

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

    tasks = ["extract","plot"]
    if len(sys.argv) > 1:
        tasks = sys.argv[1:]

    for task in tasks:
        if task == "extract":
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
        elif task == "plot":
            # Plot new data
            axes = [None, None, None]
            axes[0] = bathy.plot(new_bathy, contours=[0.0,1.0,2.0,3.0,4.0,5.0], coastlines=False)
            axes[1] = bathy.plot(base_bathy_file, contours=numpy.linspace(-10,5,16), coastlines=False)
            axes[2] = bathy.plot('./mexican_coast_pacific.tt3', contours=[0.0,1.0,2.0,3.0,4.0,5.0], coastlines=False)

            # Plot acapulco gauge location
            for axis in axes:
                axis.plot([-99.90333333], [16.83833333], 'ro')

            plt.show()
