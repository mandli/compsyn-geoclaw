#!/usr/bin/env python

r"""Script extracts topography information from raw data files."""

import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mpl_toolkits.basemap.pyproj as pyproj

import clawpack.geoclaw.topo as topo
import clawpack.geoclaw.topotools as topotools


def UTM2longlat(x, y, zone='14'):
    r"""Convert from UTM to longitude-latitude coordinates."""
    projector = pyproj.Proj(proj='utm', zone=zone)
    longitude, latitude = projector(x, y, inverse=True)

    return longitude, latitude

# Parameters
# Output paths
orig_bathy = "acapulco30m.xyz"
out_bathy = "acapulco_projected_30m.tt3"
fill_bathy = 'srtm_17_09.tt3'

# Read in transform data
print "Reading fine data..."
orig_data = topo.Topography(orig_bathy, unstructured=True)
orig_data.x, orig_data.y = UTM2longlat(orig_data.x, orig_data.y)

orig_data.z = - orig_data.z
print "Done."

# Read in fill data
print "Reading fill data..."
fill_data = topo.Topography(fill_bathy)
# Explicitly read in data to add filter
fill_data.read(mask=True, filter_region=(-100.00041620659999, 
                                         -98.788687206402713,
                                         16.462378038046992,
                                         17.158344215509878))
print "Done."

# Project data into new grid and write it out
print "Projecting unstructured data..."
# Specify custom region for projection
filter_region = list(orig_data.extent)
buffer_degrees = topotools.dist_meters2latlong(1000.0, 0.0, numpy.mean(orig_data.y))[0]
filter_region[0] -= buffer_degrees
filter_region[1] += buffer_degrees
filter_region[2] = 16.8128
filter_region[3] += buffer_degrees
orig_data.project_unstructured(fill_data.X, fill_data.Y, fill_data.Z, 
                                                           extent=filter_region)
print "Done."

print "Writing out data..."
orig_data.write(out_bathy)
fill_data.write('./srtm_subsection.tt3')
print "Done."

#  Plot new data
fig = plt.figure(figsize=[8.0*2, 6.0])
axes = fig.add_subplot(1,2,1)
mean_lat = 0.5 * (orig_data.extent[3] - orig_data.extent[2])
axes.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))    
axes.set_autoscale_on(True)
axes.set_title("Acapulco Bathymetry")
axes.set_xlabel("Longitude (degrees)")
axes.set_ylabel("Latitude (degrees)")
axes.ticklabel_format(format="plain", useOffset=False)

# Create colorbar
cmap = plt.get_cmap("terrain")
bathy_limits = (numpy.min(orig_data.Z), numpy.max(orig_data.Z))
color_norm = colors.Normalize(bathy_limits[0], bathy_limits[1], clip=True)
plot = axes.contourf(orig_data.X, orig_data.Y, orig_data.Z, 
                     vmin=bathy_limits[0], vmax=bathy_limits[1],
                     cmap=cmap, norm=color_norm)
axes.scatter(orig_data.x, orig_data.y, alpha=0.25)
axes.set_xlim(orig_data.extent[0:2])
axes.set_ylim(orig_data.extent[2:])
cbar = fig.colorbar(plot, ax=axes)
cbar.set_label("Depth (m)")

# Plot SRTM subsection
axes = fig.add_subplot(1,2,2)
mean_lat = 0.5 * (orig_data.extent[3] - orig_data.extent[2])
axes.set_aspect(1.0 / numpy.cos(numpy.pi / 180.0 * mean_lat))    
axes.set_autoscale_on(True)
axes.set_title("SRTM Topography")
axes.set_xlabel("Longitude (degrees)")
axes.set_ylabel("Latitude (degrees)")
axes.ticklabel_format(format="plain", useOffset=False)
bathy_limits = (numpy.min(fill_data.Z), numpy.max(fill_data.Z))
color_norm = colors.Normalize(bathy_limits[0], bathy_limits[1], clip=True)
axes.contourf(fill_data.X, fill_data.Y, fill_data.Z, vmin=bathy_limits[0], 
                                      vmax=bathy_limits[1],
                                      cmap=cmap, 
                                      norm=color_norm)
axes.set_xlim(fill_data.extent[0:2])
axes.set_ylim(fill_data.extent[2:])
cbar = fig.colorbar(plot, ax=axes)
cbar.set_label("Depth (m)")

plt.show()