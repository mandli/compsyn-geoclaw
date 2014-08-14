#!/usr/bin/env python

r"""Script extracts topography information from raw data files."""

import sys
import subprocess

import numpy
import matplotlib.pyplot as plt
import mpl_toolkits.basemap.pyproj as pyproj
import matplotlib.mlab as mlab

import clawpack.geoclaw.topotools as topotools

def UTM2longlat(x, y, zone='14'):
    r"""Convert from UTM to longitude-latitude coordinates."""
    projector = pyproj.Proj(proj='utm', zone=zone)
    longitude, latitude = projector(x, y, inverse=True)

    return longitude, latitude


def extract_bathy(fill_paths, region, topo_path='./bathymix.xyz', 
                                      out_bathy="new_acapulco_bathy.tt2"):

    # Read in concatenated data and transform it
    print "Reading fine data..."
    orig_data = topotools.Topography(topo_path, unstructured=True)
    orig_data.x, orig_data.y = UTM2longlat(orig_data.x, orig_data.y)
    orig_data.z = -orig_data.z
    print "Done."

    # Project data into new grid and write it out
    # Specify custom region for projection
    filter_region = list(orig_data.extent)
    buffer_degrees = topotools.dist_meters2latlong(1e3, 0.0, numpy.mean(orig_data.y))[0]
    filter_region[0] -= buffer_degrees
    filter_region[1] += buffer_degrees
    # filter_region[2] = 16.8128
    filter_region[2] -= buffer_degrees
    filter_region[3] += buffer_degrees

    # Create fill topography objects using filters and transforming coordinates
    print "Reading fill topography..."
    fig = plt.figure()
    axes = fig.add_subplot(1,1,1)
    fill_topo = [topotools.Topography(path, unstructured=True) for path in fill_paths]
    for (n, topo) in enumerate(fill_topo):
        topo.no_data_value = 0.0
        topo.read(mask=True)
        topo.x, topo.y = UTM2longlat(topo.x, topo.y)
        indices = mlab.find((numpy.abs(topo.z) > 1e-2))
        print "Throwing out %s/%s points." % (indices.shape[0],topo.x.shape[0])
        topo.x = topo.x.take(indices)
        topo.y = topo.y.take(indices)
        topo.z = topo.z.take(indices)
        if n == 0:
            topo.plot(axes=axes, region_extent=region, add_colorbar=True)
        else:
            topo.plot(axes=axes, region_extent=region, add_colorbar=False)
    plt.savefig("fill_topo.png")
    print "Done reading fill topography..."

    print "Interpolating unstructured data..."
    orig_data.interp_unstructured(fill_topo, extent=region, 
                                              delta_limit=20.0)
    orig_data.write(out_bathy, topo_type=2)
    print "Done projecting data."

    # Plot final topography
    print "Plotting final topography."
    acapulco_final_topo = topotools.Topography(out_bathy)
    acapulco_final_topo.plot()

    plt.savefig("final_topo.png")
    # subprocess.Popen('touch final_topo.png', shell=True)


if __name__ == '__main__':

    if len(sys.argv) > 1:
        pass

    fill_paths = ["./702825794163_xyz/MT_XYZ/E14C57D2_MT.xyz", 
                  "./702825794200_xyz/MT_XYZ/E14C57E1_MT.xyz", 
                  "./702825794286_xyz/MT_XYZ/E14C57E3_MT.xyz", 
                  "./702825794323_xyz/MT_XYZ/E14C57E4_MT.xyz",
                  "./reskeleton/E14C57E2_MT.xyz", 
                  "./reskeleton/E14C57F1_MT.xyz",
                  "./reskeleton/E14C57F3_MT.xyz",
                  "./reskeleton/E14C67C1_MT.xyz",
                  "./reskeleton/E14C67C2_MT.xyz"]
                  # "../mexican_coast_pacific.tt3"]
    # fill_paths = ["../srtm_17_09.tt3"]
    extract_bathy(fill_paths, (-99.930021, -99.730477, 16.710640, 16.870122))
    # extract_bathy(fill_paths, (-99.930021, -99.830477, 16.780640, 16.870122))