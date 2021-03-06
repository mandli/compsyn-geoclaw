#!/usr/bin/env python
"""
Plot fgmax output from GeoClaw runs, assuming points are on a rectangular
grid.

"""

import os

import numpy

import matplotlib.pyplot as plt

import clawpack.clawutil.data as clawdata
import clawpack.geoclaw.data as geodata

def make_fgmax_plots(output_path="_output", data_path="_output", 
                     plot_path="_plots", plot_zeta=True, plot_zeta_times=True, 
                     plot_arrival_times=True):
    r""""""

    # Load appropriate data files
    claw_data = clawdata.ClawpackInputData()
    claw_data.read(os.path.join(data_path, 'claw.data'))
    geo_data = geodata.GeoclawInputData()
    geo_data.read(os.path.join(data_path, 'geo.data'))

    # Parse fgmax input data file
    fgmax_input = open(os.path.join(data_path, 'fgmax_data.txt'), 'r')
    num_cells = [int(value) for value in fgmax_input.readline().split()]
    min_level = 2

    # Parameters
    zeta_contours = []
    t_contours = []
    topo_contours = []

    # Load fgmax output data
    fg_data = numpy.loadtxt(os.path.join(output_path, 'fort.FG1.valuemax'))
    
    x = numpy.reshape(fg_data[:,0], (mx,my), order='F')
    y = numpy.reshape(fg_data[:,1], (mx,my), order='F')
    y0 = 0.5 * (y.min() + y.max())   # mid-latitude for scaling plots
    eta_tilde = numpy.reshape(fg_data[:,3], (mx,my), order='F')
    level = numpy.reshape(fg_data[:,2].astype('int'), (mx,my), order='F')

    # Load fgmax bathymetry data
    fg_aux_data = numpy.loadtxt(os.path.join(output_path, 'fort.FG1.aux1'))
    for n in xrange(min_level)



def make_plots(outdir='_output', plotdir='_plots'):

    # Some things that might need to change...
    fgmax_input_file = 'fgmax_grid.txt'

    plot_zeta = True
    plot_zeta_times = True
    plot_arrival_times = True

    sea_level = 0.

    #clines_zeta = None  # can set to desired contours of zeta 
    clines_zeta = [0.01] + list(linspace(0.05,0.3,6)) + [0.5,1.0,10.0]

    #clines_t = None  # can set to desired contours of arrival time or zeta time
    clines_t = linspace(0,8,17)  # hours
    clines_t_label = clines_t[::2]  # which ones to label 
    clines_t_colors = [.5,.5,.5]
    clines_topo = [0]

    if not os.path.isdir(outdir):
        raise Exception("Missing directory: %s" % outdir)

    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)

    print outdir
    print fgmax_input_file

    # read mx and my from the input file:
    try:
        fid = open(fgmax_input_file)
    except:
        raise Exception("cannot open %s" % fgmax_input_file)

    # skip some lines:
    for i in range(5):
        line = fid.readline()

    line = fid.readline().split()
    fid.close()
    mx = int(line[1])
    my = int(line[2])


    fname = outdir + '/fort.FG1.valuemax' 
    print "Reading %s ..." % fname
    try:
        d = loadtxt(fname)
    except:
        raise Exception("*** Cannot read file: %s" % fname)

    x = reshape(d[:,0],(mx,my),order='F')
    y = reshape(d[:,1],(mx,my),order='F')
    y0 = 0.5*(y.min() + y.max())   # mid-latitude for scaling plots
    eta_tilde = reshape(d[:,3],(mx,my),order='F')

    # AMR level used for each zeta value:
    level = reshape(d[:,2].astype('int'),(mx,my),order='F')
    
    # Determine topo B at each point from the same level of AMR:
    fname = outdir + '/fort.FG1.aux1' 
    print "Reading %s ..." % fname
    daux = loadtxt(fname)
    topo = []
    nlevels = daux.shape[1]
    for i in range(2,nlevels):
        topoi = reshape(daux[:,i],(mx,my),order='F')
        topoi = ma.masked_where(topoi < -1e50, topoi)
        topo.append(topoi)

    B = ma.masked_where(level==0, topo[0])  # level==0 ==> never updated
    levelmax = level.max()
    for i in range(levelmax):
        B = where(level==i+1, topo[i], B)

    h = where(eta_tilde > B, eta_tilde - B, 0.)

    # zeta = max h on land or max eta offshore:
    zeta = where(B>sea_level, h, eta_tilde)
    zeta = where(zeta > -1e20, zeta, sea_level)

    tzeta = reshape(d[:,4],(mx,my),order='F')  # Time maximum h recorded
    tzeta = ma.masked_where(tzeta < -1e50, tzeta)      
    tzeta = ma.masked_where(zeta == 0., tzeta) / 3600.  # hours 

    inundated = logical_and((B>0), (h>0))

    atimes = reshape(d[:,5],(mx,my),order='F')
    atimes = ma.masked_where(atimes < -1e50, atimes)  
    atimes = ma.masked_where(zeta == 0., atimes) / 3600.  # hours 

    if plot_zeta:

        # Plot h or eta along with contours of topo:
        figure(101)
        clf()
        zeta = ma.masked_where(zeta==0.,zeta)
        if clines_zeta is None:
            cmax = zeta.max()
            cmin = zeta.min()
            clines_zeta = linspace(cmin,cmax,10)
        colors = discrete_cmap(clines_zeta)
        contourf(x,y,zeta,clines_zeta,colors=colors)

        cbar = colorbar()
        cbar.set_ticks(clines_zeta)
        cbar.set_label('meters', fontsize=15)


        # Contours of topo:
        contour(x,y,B,clines_topo,colors='g',linestyles='-')

        if plot_arrival_times:
            # Contours of arrival time
            cs = contour(x,y,atimes,clines_t,colors=clines_t_colors)
            clabel(cs,clines_t_label)


        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))

        title("Zeta Maximum",fontsize=20)
        if plot_arrival_times:
            title("Zeta Maximum and arrival times",fontsize=15)
        
        fname = plotdir + '/zeta.png' 
        savefig(fname)
        print "Created ",fname


    if plot_zeta_times:

        # Plot time max h recorded:
        figure(102)
        clf()

        if clines_t is None:
            clines_t = linspace(tzeta.min(), tzeta.max(), 10)
            

        colors = discrete_cmap_times(clines_t)
        contourf(x,y,tzeta,clines_t,colors=colors)
        cbar = colorbar()
        cbar.set_ticks(clines_t)
        cbar.set_label('hours',fontsize=15)

        # Contours of topo:
        contour(x,y,B,clines_topo,colors='k',linestyles='-')

        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))
        title('Time of max zeta', fontsize=20)
        
        fname = plotdir + '/zetatimes.png' 
        savefig(fname)
        print "Created ",fname


    if plot_arrival_times:

        # Plot time max h recorded:
        figure(103)
        clf()
        if clines_t is None:
            clines_t = linspace(atimes.min(), atimes.max(), 10)
        colors = discrete_cmap_times(clines_t)
        contourf(x,y,atimes,clines_t,colors=colors)
        cbar = colorbar()
        cbar.set_ticks(clines_t)
        cbar.set_label('hours',fontsize=15)

        # Contours of topo:
        contour(x,y,B,clines_topo,colors='k',linestyles='-')
        #
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))
        title('Arrival time', fontsize=20)
        fname = plotdir + '/arrival_times.png' 
        savefig(fname)
        print "Created ",fname



def discrete_cmap(clines):
    """
    Construct a discrete color map for the regions between the contour lines
    given in clines. Colors go from turqouise through yellow to red.
    """
    nlines = len(clines)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
    Red = hstack([linspace(0,0.8,n1), ones(n2)])
    Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
    colors = zip(Red,Green,Blue)
    return colors

def discrete_cmap_times(clines):
    """
    Construct a discrete color map for the regions between the contour lines
    given in clines. For arrival times, colors go from red to turquoise.
    """
    nlines = len(clines)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = flipud(hstack([linspace(1,1,n1),linspace(1,0,n2)]))
    Red = flipud(hstack([linspace(0,0.8,n1), ones(n2)]))
    Blue = flipud(hstack([linspace(1,0.2,n1), zeros(n2)]))
    colors = zip(Red,Green,Blue)
    return colors

if __name__ == "__main__":
    make_plots()
