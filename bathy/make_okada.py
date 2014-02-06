#!/usr/bin/env python

import sys

import numpy
import matplotlib.pyplot as plt

import clawpack.geoclaw.okada2 as okada
import clawpack.geoclaw.topo as topo
import clawpack.geoclaw.dtopotools as dtopotools
import clawpack.visclaw.colormaps as colormaps

plot_fault = False
if len(sys.argv) > 1:
    if sys.argv[1] == "plot":
        plot_fault = True

# Base subfault
# M 7.8
# L 90
# W 70
# depth of hypocenter  20 km (it is written in the text)
# coordinates of epicenter: 17.110 -99.100 (taken from: http://usuarios.geofisica.unam.mx/vladimir/sismos/100a%F1os.html; here they give a hypocentral depth of 33, but let's try with 20, as mentioned aove)
# strike 296
# dip 25
# rake 90
# top of the fault plane at 12 km.
# They use a value of rigidity of 5 X 10^11 dyna/cm^2 -> 5e11 * 1e2 = 5e13
subfaults = []
test_slips = [165, 200]
for (n, slip) in enumerate(test_slips):
    subfaults.append(topo.SubFault(units={"slip":"cm", "dimensions":"km", "depth":"km"}))
    subfaults[-1].coordinates = [-99.25, 16.6]
    subfaults[-1].coordinate_specification = "top center"
    subfaults[-1].slip = slip
    subfaults[-1].rake = 90.0
    subfaults[-1].strike = 296
    subfaults[-1].dip = 25.0
    subfaults[-1].depth = 12.0
    subfaults[-1].dimensions = (90.0, 70.0)
    subfaults[-1].write('./okada_1957_du%s.tt3' % slip)
    
    print "Subfault Characteristics:"
    print "  Mw = %s" % str(subfaults[-1].Mw(mu=5e11))
    print "  du = %s" % subfaults[-1].slip
    print "  Containing Rect = %s" % subfaults[-1].containing_rect()

if plot_fault:
    fig = plt.figure()
    for (n,subfault) in enumerate(subfaults):
        axes = fig.add_subplot(1,2,n+1)
        subfault.plot(axes)
        subfault.plot_fault_rect(axes, color='k')
        subfault.plot_rake(axes)
        axes.set_title(r"$M_w = %s$, $\Delta u = %s$ cm" % (subfault.Mw(mu=5e11), 
                                                            subfault.slip))
    plt.show()