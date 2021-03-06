#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os

import numpy

import clawpack.geoclaw.dtopotools as dtopo

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')


    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = -118.0      # west longitude
    clawdata.upper[0] = -86.0       # east longitude

    clawdata.lower[1] = 7.0          # south latitude
    clawdata.upper[1] = 21.0         # north latitude


    # Number of grid cells
    res_factor = 1
    clawdata.num_cells[0] = 32 * res_factor
    clawdata.num_cells[1] = 14 * res_factor

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 4

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00036'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 2

    if clawdata.output_style==1:
        # Output nout frames at equally spaced times up to tfinal:
        hours = 4
        output_per_hour = 12
        clawdata.num_output_times = hours * output_per_hour
        clawdata.tfinal = float(hours) * 3600.0
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [float(time) for time in xrange(0,250,25)]
        acapulco_time_zoom = numpy.linspace(0.07, 0.2, 10) * 3600.0
        for new_time in acapulco_time_zoom:
            clawdata.output_times.append(new_time) 
        hours = 4
        output_per_hour = 6
        for time in numpy.linspace(0, hours * 3600, output_per_hour * hours + 1):
            if time > clawdata.output_times[-1]:
                clawdata.output_times.append(float(time))

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True
        

    clawdata.output_format = 'binary'      # 'ascii' or 'netcdf' 

    clawdata.output_q_components = 'all'   # need all
    clawdata.output_aux_components = 'none'  # eta=h+B is in q
    clawdata.output_aux_onlyonce = False    # output aux arrays each frame



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 3



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 2.

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000




    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 1

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 7

    # List of refinement ratios at each level (length at least mxnest-1)
    #  Level 1 - (1.0º,1.0º) - (100950.057205 m,110772.872596 m)
    #  Level 2 - (0.25º,0.25º) - (25237.5143013 m,27693.2181489 m)
    #  Level 3 - (0.125º,0.125º) - (12618.7571506 m,13846.6090744 m)
    #  Level 4 - (0.0625º,0.0625º) - (6309.37857532 m,6923.30453722 m)
    #  Level 5 - (0.0104166666667º,0.0104166666667º) - (1051.56309589 m,1153.88408954 m)
    #  Level 6 - (0.00130208333333º,0.00130208333333º) - (131.445386986 m,144.235511192 m)
    #  Level 7 - (8.13802083333e-05º,8.13802083333e-05º) - (8.21533668662 m,9.01471944951 m))
    amrdata.refinement_ratios_x = [4,2,2,6,8,8]
    amrdata.refinement_ratios_y = [4,2,2,6,8,8]
    amrdata.refinement_ratios_t = [4,2,2,6,8,8]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft','center']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0  

    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = True       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    rundata.regiondata.regions.append([7,7,0.0,1e10, -99.930021, -99.830477, 16.780640, 16.870122])
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    # Region                Long (E)            Lat (N)
    # Acapulco              -99º 52' 41.10"     16º 50' 18.19"
    #rundata.regiondata.regions.append([1, 6, 0.0, 1e10,
    #                                    -100.1, -99.66666667,
    #                                     16.7, 16.96666667])

    # Ixtapa-Zihuatanejo    -101º 33' 8.61"     17º 38' 15.15"
    # Puerto Angel          -96º 29' 35.08"     15º 39' 53.28"

    # Lázaro Cárdenas       -102º 9' 54.86"     17º 55' 30.66"
    #rundata.regiondata.regions.append([1, 6, 0.0, 1e10,
    #                                    -102.2440361, -102.0918583,
    #                                     17.89015556,  17.99216667])

    # ---------------
    # Gauges:
    # ---------------
    degminsec2dec = lambda deg, minutes, seconds: float(deg) + (float(minutes) + float(seconds) / 60.0) / 60.0
    rundata.gaugedata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    # ID      Location Name             Lat               Long
    # 1       Manzanillo, Col.          19º 3.4 N         104º 19.1 W
    # rundata.gaugedata.gauges.append([1, -104.3183333, 19.05666667, 0.0, 1.e10])
    # 2       Ixtapa, Gro.              17º 40.1 N        101º 38.7 W
    # rundata.gaugedata.gauges.append([2, -101.645, 17.66833333, 0.0, 1.e10])
    # 3       Zihuatanejo, Gro.         17º 38.2 N        101º 33.5 W
    # rundata.gaugedata.gauges.append([3, -101.5583333, 17.63666667, 0.0, 1.e10])
    # 4       Acapulco, Gro.            16º 50.3 N        99º 54.2 W
    rundata.gaugedata.gauges.append([4, -99.90333333, 16.83833333, 0.0, 1.e10])
    # 5       Isla Socorro, Col.        18º 43.5 N        110º 57.0 W
    # rundata.gaugedata.gauges.append([5, -110.95, 18.725, 0.0, 1.e10])
    # 6       Puerto Angel, Oax         15º 40 N          96º 29.5 W
    # rundata.gaugedata.gauges.append([6, -degminsec2dec(96,29.5,0), degminsec2dec(15,40,0), 0.0, 1.e10])
    # 7       Salina Cruz, Oax.         16º 19.1 N        95º 11.8 W
    rundata.gaugedata.gauges.append([7, -degminsec2dec(95,11.8,0), degminsec2dec(16,19.1,0.0), 0.0, 1.e10])
    # 8       Puerto Madero, Chis       14º 42.7 N        92º 24.1 W
    # rundata.gaugedata.gauges.append([8, -degminsec2dec(92,24.1,0), degminsec2dec(14,42.7,0), 0.0, 1.e10])
    # 9       Lazaro Cardenas, Mich     17º 56.4 N        102º 10.7 W
    # rundata.gaugedata.gauges.append([9, -degminsec2dec(102,10.7,0.0), degminsec2dec(17,56.4,0.0), 0.0, 1.e10])
    # 10      Huatulco                  15° 45'.2 N       96° 07'.8 W     
    # rundata.gaugedata.gauges.append([10, -degminsec2dec(96,7.8,0.0), degminsec2dec(15,45.2,0.0), 0.0, 1.e10])
    # 11      Acapulco                  16°51'9.00"N      99°52'50"W 
    rundata.gaugedata.gauges.append([11, -degminsec2dec(99,52,50), degminsec2dec(16,51,9), 0.0, 1.e10])
    # 12      Acapulco additional gauges
    rundata.gaugedata.gauges.append([12, -99.904294, 16.839721, 0.0, 1.e10])
    rundata.gaugedata.gauges.append([13, -99.905197, 16.840743, 0.0, 1.e10])
    rundata.gaugedata.gauges.append([14, -99.903940, 16.842113, 0.0, 1.e10])
    rundata.gaugedata.gauges.append([15, -99.902489, 16.843462, 0.0, 1.e10])
    rundata.gaugedata.gauges.append([16, -99.898397, 16.845365, 0.0, 1.e10])
    rundata.gaugedata.gauges.append([17, -99.891848, 16.851036, 0.0, 1.e10])
    rundata.gaugedata.gauges.append([18, -99.860943, 16.848830, 0.0, 1.e10])
    rundata.gaugedata.gauges.append([19, -99.856680, 16.839136, 0.0, 1.e10])
    rundata.gaugedata.gauges.append([20, -99.888627, 16.816910, 0.0, 1.e10])

    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print "*** Error, this rundata has no geo_data attribute"
        raise AttributeError("Missing geo_data attribute")
       
    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = [0.03, 0.025]
    geo_data.manning_break = [0.0]
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.25
    refinement_data.deep_depth = 200.0
    refinement_data.max_level_deep = 5

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]
    topo_data.topofiles.append([3, 1, 5, 0., 1.e10, 
                          os.path.abspath('./bathy/mexican_coast_pacific.tt3')])
    # topo_data.topofiles.append([3, 6, 7, 0., 1.e10, 
    #                       os.path.abspath('./bathy/acapulco_projected_30m.tt2')])
    topo_data.topofiles.append([2, 5, 7, 0., 1.e10, 
                          os.path.abspath('./bathy/new_bathy/new_acapulco_bathy.tt2')])

    # topo_data.topofiles.append([3, 1, 10, 0., 1.e10, 
    #                       os.path.abspath('./bathy/srtm_subsection.tt3')])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data
    # for moving topography, append lines of the form :
    #   [topotype, minlevel,maxlevel,fname]
    dtopo_data.dtopofiles.append([1, 5, 5, 'bathy/rot_gapSvr1zvT.xyzt'])
    #dtopo_data.dtopofiles.append([3, 5, 5, 'okada_1957Sm_du370.tt3'])
    # subfault = dtopo.SubFault(units={"slip":"cm", "dimensions":"km", "depth":"km"})
    # subfault.coordinates = [-99.25, 16.6]
    # subfault.coordinate_specification = 'top center'
    # subfault.slip = 200
    # subfault.rake = 90.0
    # subfault.strike = 296
    # subfault.dip = 15.0
    # subfault.depth = 4.0
    # subfault.dimensions = (320.0, 80.0)
    # subfault.my = 5e11
    # subfault.write('./dtopo.tt3')
    # dtopo_data.dtopofiles.append([3,5,5,'dtopo.tt3'])
    # Note that if the run_faults.py script is used this is overriden there



    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    # == setfixedgrids.data values ==
    fixed_grids = rundata.fixed_grid_data
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    # == fgmax.data values ==
    fgmax_files = rundata.fgmax_data.fgmax_files
    # for fixed grids append to this list names of any fgmax input files
    # fgmax_files.append(os.path.abspath(os.path.join(os.getcwd(),'fgmax_grid.txt')))

    return rundata
    # end of function setgeo
    # ----------------------



if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()

