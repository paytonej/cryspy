# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 16:59:33 2013

@author: epayton
"""

import ovlib as ov
from numpy import pi

# 4.3.1
uc = ov.unitcell(a=3., b=4., c=6., beta=120.)
bodydiag = ov.lattvec(1,1,1)
print bodydiag.length(uc) # result should be 6.557

# 4.3.2
uc = ov.unitcell(a=2., b=2., c=3.)
atom1 = ov.lattsite(1./2., 1./3., 1./4.)
atom2 = ov.lattsite(1./3., 1./2., 3./4.)
print atom1.distance(atom2, uc) # result should be 1.572

# 4.3.3
uc = ov.unitcell()
atom1 = ov.lattsite(1./2., 1./2., 0.)
atom2 = ov.lattsite(1./2., 0., 1./2.)
origin = ov.lattsite(0., 0., 0.)
print atom1.angle(atom2, uc, origin)*180./pi # result should be 1.572

# 4.3.4
uc = ov.unitcell(a=4., b=6., c=5., beta=120.)
vec1 = ov.lattvec(1., 0., 1.)
vec2 = ov.lattvec(-2., 0., 1.)
print vec1.angle(vec2, uc)*180./pi

# 6.5(i)
uc = ov.unitcell(a=2., b=2., c=2.)
plane = ov.miller(1., 1., 0.)
print plane.dspacing(uc)

# 6.5(ii)
uc = ov.unitcell(a=3., b=4., c=6., gamma=120.)
plane = ov.miller(1., 1., 1.)
print plane.dspacing(uc)

# 6.5(iii)
uc = ov.unitcell(a=4., b=6., c=5., beta=120.)
plane1 = ov.miller( 1., 0., 1.)
plane2 = ov.miller(-2., 0., 1.)
print plane1.angle(plane2, uc)*180./pi

# Table 7.2
uc = ov.unitcell(beta=45.)
u = [1., 0., 0., 1., 1., 0., 1.]
v = [0., 1., 0., 1., 0., 1., 1.]
w = [0., 0., 1., 0., 1., 1., 1.]
vec1 = ov.lattvec(u, v, w)
cart1 = vec1.to_cartesian(uc)
xproj1, yproj1, hemi1 = ov.stereotrans(cart1)

# Table 7.2
vec2 = ov.miller(u, v, w)
cart2 = vec2.to_cartesian(uc)
xproj2, yproj2, hemi2 = ov.stereotrans(cart2)

# Fig 7.11(a)
sg = ov.stereoproj() # initialize the stereographic projection
sg.add_lattvec(vec1,uc) # add our lattice vectors to the projection
sg.add_lattveclabels(vec1,uc)
sg.grid_greatcircles()
sg.grid_smallcircles()
sg.grid_spokes(degreestep=45)

from matplotlib.font_manager import FontProperties
font=FontProperties()
font.set_size('large')
font.set_style('italic')
sg.add_eastlabel(fontproperties=font)
sg.add_southlabel(fontproperties=font)
sg.add_centerlabel(fontproperties=font)
sg.add_coordinatereadout(uc)

# Fig 7.11(b)
sg = ov.stereoproj() # initialize the stereographic projection
sg.add_miller(vec2,uc) # add our lattice vectors to the projection
sg.add_millerlabels(vec2,uc)
sg.grid_wulffnet()

from matplotlib.font_manager import FontProperties
font=FontProperties()
font.set_size('large')
font.set_style('italic')
sg.add_eastlabel(fontproperties=font)
sg.add_southlabel(fontproperties=font)
sg.add_centerlabel(fontproperties=font)
sg.add_coordinatereadout(uc)
# NOTE THAT (111) and (110) are in the wrong places in the figure in the text!

## Fig 7.11(b) as equal area proj
sg = ov.eaproj() # initialize the stereographic projection
sg.add_miller(vec2,uc) # add our lattice vectors to the projection
sg.add_millerlabels(vec2,uc)
sg.grid_spokes()
#
from matplotlib.font_manager import FontProperties
font=FontProperties()
font.set_size('large')
font.set_style('italic')
sg.add_eastlabel(fontproperties=font)
sg.add_southlabel(fontproperties=font)
sg.add_centerlabel(fontproperties=font)
sg.add_coordinatereadout(uc)