
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import sys

import numpy as np

# Plot customization
import matplotlib

# Markers and line widths
matplotlib.rcParams['lines.linewidth'] = 2.0
matplotlib.rcParams['lines.markersize'] = 6
matplotlib.rcParams['lines.markersize'] = 8

# Font Sizes
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['axes.labelsize'] = 15
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12

# DPI of output images
matplotlib.rcParams['savefig.dpi'] = 100

# Change backend
# This does not work as pylab is imported before we are called
# matplotlib.use('Agg') 

import matplotlib.pyplot as plt

from clawpack.visclaw import geoplot
import clawpack.clawutil.data
import clawpack.amrclaw.data
import clawpack.geoclaw.dtopotools as dtopo

if not (sys.platform == 'darwin'):
    print "Switching to", plt.switch_backend('Agg')

try:
    from setplotfg import setplotfg
except:
    print "Did not find setplotfg.py"
    setplotfg = None

# Try and load Acapulco gauge data for 1957 quake
try:
    gauge_path = "./gauge_obs/fig2.txt"
    gauge_data = np.loadtxt(gauge_path, skiprows=1, delimiter=',')
except IOError:
    pass

# Create subfault for plotting
try:    
    subfault = dtopo.SubFault(units={"slip":"cm", "dimensions":"km", "depth":"km"})
    subfault.coordinates = [-99.25, 16.6]
    subfault.coordinate_specification = "top center"
    subfault.slip = 200
    subfault.rake = 90.0
    subfault.strike = 296
    subfault.dip = 25.0
    subfault.depth = 12.0
    subfault.dimensions = (90.0, 70.0)
except:
    subfault = None


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'

    data = clawpack.clawutil.data.ClawInputData(2)
    data.read('./claw.data')
    regiondata = clawpack.amrclaw.data.RegionData()
    regiondata.read('./regions.data')
    xlimits = (data.lower[0], data.upper[0])
    ylimits = (data.lower[1], data.upper[1])
    surface_limit = 1.0
    speed_limit = 1.0

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    
    # ========================================================================
    #  Water helper functions
    # ========================================================================
    def b(cd):
        return cd.q[3,:,:] - cd.q[0,:,:]
        
    def extract_eta(h,eta,DRY_TOL=10**-3):
        index = np.nonzero((np.abs(h) < DRY_TOL) + (h == np.nan))
        eta[index[0],index[1]] = np.nan
        return eta
    
    def extract_velocity(h,hu,DRY_TOL=10**-8):
        u = np.zeros(hu.shape)
        index = np.nonzero((np.abs(h) > DRY_TOL) * (h != np.nan))
        u[index[0],index[1]] = hu[index[0],index[1]] / h[index[0],index[1]]
        return u
    
    def eta(cd):
        return extract_eta(cd.q[0,:,:],cd.q[3,:,:])
        
    def water_u(cd):
        return extract_velocity(cd.q[0,:,:],cd.q[1,:,:])
        
    def water_v(cd):
        return extract_velocity(cd.q[0,:,:],cd.q[2,:,:])
        
    def water_speed(current_data):
        u = water_u(current_data)
        v = water_v(current_data)
            
        return np.sqrt(u**2+v**2)


    #-----------------------------------------
    # Figure for surface
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Full Surface')

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Surface at %4.2f hours' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)

        # Plot fault plane
        subfault.plot_fault_rect(pylab.gca())
        subfault.plot_rake(pylab.gca())

    plotaxes.afteraxes = fixup

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -surface_limit
    plotitem.pcolor_cmax = surface_limit
    plotitem.add_colorbar = True
    plotitem.colorbar_label = "Surface (m)"
    plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
    plotitem.patchedges_show = 0
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = np.linspace(-3000,-3000,1)
    plotitem.amr_contour_colors = ['y']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = 0
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    #-----------------------------------------
    # Zoom in of shore
    #-----------------------------------------
    # Note that for the region names to come out correctly they either have 
    # to be in the html versions or latex versions as these are just titles
    # region_names = [r'Acapulco', r"Acapulco Zoom"]
    #regiondata.regions.append([1,10,0.0,1e10,-99.9,-99.8516,16.7728,16.8589])
    regiondata.regions.append([1, 10, 0.0, 1e10, 
                      -99.9364923773, -99.80, 16.7714776209, 16.86007248])
    regiondata.regions.append([1, 10, 0.0, 1e10, 
                      -99.9797, -99.8335, 16.72, 16.876])
    regiondata.regions.append([1, 10, 0.0, 1e10, 
                      -99.9119, -99.8961, 16.8366, 16.8511])
    for (n,region) in enumerate(regiondata.regions):
        plotfigure = plotdata.new_plotfigure()#name='Surface Zoom of %s' % region_names[n])

        # Set up for axes in this figure:
        plotaxes = plotfigure.new_plotaxes('pcolor')
        plotaxes.title = 'Surface'
        plotaxes.scaled = True

        def fixup(current_data):
            import pylab
            addgauges(current_data)
            t = current_data.t
            t = t / 3600.  # hours
            pylab.title('Surface at %4.2f hours' % t, fontsize=20)
            pylab.xticks(fontsize=15)
            pylab.yticks(fontsize=15)
        plotaxes.afteraxes = fixup

        # Water
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        # plotitem.plot_var = geoplot.surface
        plotitem.plot_var = geoplot.surface_or_depth
        plotitem.pcolor_cmap = geoplot.tsunami_colormap
        plotitem.pcolor_cmin = -surface_limit
        plotitem.pcolor_cmax = surface_limit
        plotitem.add_colorbar = True
        plotitem.colorbar_label = "Surface (m)"
        plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
        plotitem.patchedges_show = 0

        # Land
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        plotitem.plot_var = geoplot.land
        plotitem.pcolor_cmap = geoplot.land_colors
        plotitem.pcolor_cmin = 0.0
        plotitem.pcolor_cmax = 100.0
        plotitem.add_colorbar = False
        plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
        plotitem.patchedges_show = 0
        plotaxes.xlimits = region[4:6]
        plotaxes.ylimits = region[6:9]

        # add contour lines of bathy if desired:
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.show = False
        plotitem.plot_var = geoplot.topo
        plotitem.contour_levels = np.linspace(-3000,-3000,1)
        plotitem.amr_contour_colors = ['y']  # color on each level
        plotitem.kwargs = {'linestyles':'solid','linewidths':2}
        plotitem.amr_contour_show = [1,0,0]  
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0


    #-----------------------------------------
    # Figure for velocities
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Speeds')

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('speeds')
    plotaxes.title = 'Speeds'
    plotaxes.scaled = True

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Speeds at %4.2f hours' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)
    plotaxes.afteraxes = fixup

    # Speed
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = water_speed
    plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = speed_limit
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
    plotitem.patchedges_show = 0

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 1000.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0,0,0,0,0]
    plotitem.patchedges_show = 0
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits

    # =================================
    #  Plot Deformation with Coastline
    # =================================
    plotfigure = plotdata.new_plotfigure(name="Deformation")
    plotfigure.show = (subfault is not None)

    def add_deformation_plots(cd):
        axes = plt.gca()
        subfault.plot(axes, contours=[0.02, 0.15, 0.50])
        subfault.plot_fault_rect(axes)
        subfault.plot_rake(axes)

    if subfault is not None:
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.xlimits = [subfault.x[0], subfault.x[-1]]
        plotaxes.ylimits = [subfault.y[0], subfault.y[-1]]
        plotaxes.title = "Subfault Deformation"
        plotaxes.scaled = True
    
        plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
        plotitem.show = True
        plotitem.plot_var = geoplot.topo
        plotitem.contour_levels = [0.0]
        plotitem.amr_contour_colors = ['k']  # color on each level
        plotitem.kwargs = {'linestyles':'solid','linewidths':2}
        plotitem.amr_contour_show = [1]
        plotitem.celledges_show = 0
        plotitem.patchedges_show = 0


        plotaxes.afteraxes = add_deformation_plots

    # =======================
    #  Figure for Bathymetry
    # =======================
    plotfigure = plotdata.new_plotfigure(name="Bathymetry")
    plotfigure.show = False

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.title = "Bathymetry"
    plotaxes.scaled = True

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = b
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = -10
    plotitem.pcolor_cmax = 10
    plotitem.amr_celledges_show = [0,0,0,0,0,0,0,0,0]
    plotitem.patchedges_show = 0

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-5,5]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'

    def gauge_afteraxes(current_data):
        from pylab import plot, xticks, floor
        t = current_data.t
        #legend(('surface','topography'),loc='lower left')
        plot(t, 0*t, 'k')
        n = int(floor(t.max()/3600.) + 2)
        xticks([3600*i for i in range(n)])
        plt.xlabel("t (s)")
        plt.ylabel("Surface (m)")

        if current_data.gaugeno == 4:
            # Add gauge observations
            plot(gauge_data[:,0] * 60.0, gauge_data[:,1], 'x')

    plotaxes.afteraxes = gauge_afteraxes


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

