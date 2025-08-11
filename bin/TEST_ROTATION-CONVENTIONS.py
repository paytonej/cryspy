# -*- coding: utf-8 -*-
"""
Created on Thu Jul 04 17:05:51 2013

@author: epayton
"""

import ovlib as ov

uc = ov.unitcell()

# From Lecture notes for CMU 27-750, Fall 2009, "Analysis of EBSD Data (L17)"
# by B. El-Dasher, A.D. Rollett, G.S. Rohrer, P.N. Kalu, p. 35:
# "A simple test of the frames used for Euler angles is to have the softwares
# [sic] plot pole figures for a single orientation..."

'''
Check if Euler angles are working properly in passive convention
'''
pf = ov.stereoproj()
m = ov.miller(1,0,0)
pf.add_miller(m.rotate(ov.rmat.from_bunge(ov.bunge(0.2,0,0))), uc, uppermarkerfacecolor='r', lowermarkerfacecolor='r')
pf.add_miller(m, uc)

# the red marker should appear counterclockwise from the blue
m = ov.miller(0,0,1)
pf.add_miller(m.rotate(ov.rmat.from_bunge(ov.bunge(0.0,0.2,0))), uc, uppermarkerfacecolor='r', lowermarkerfacecolor='r')
pf.add_miller(m, uc)

'''
Check if ang/ax is working properly in passive convention
'''
pf = ov.stereoproj()
# the red marker should appear in same directions as above:
m = ov.miller(1,0,0)
pf.add_miller(m.rotate(ov.rmat.from_angax(ov.angax(0.2,0,0,1))), uc, uppermarkerfacecolor='r', lowermarkerfacecolor='r')
pf.add_miller(m, uc)
m = ov.miller(0,0,1)
pf.add_miller(m.rotate(ov.rmat.from_angax(ov.angax(0.2,1,0,0))), uc, uppermarkerfacecolor='r', lowermarkerfacecolor='r')
pf.add_miller(m, uc)

'''
Check if quat is working properly in passive convention
'''
# the red marker should appear in same directions as above:
pf = ov.stereoproj()
m = ov.miller(1,0,0)
pf.add_miller(m.rotate(ov.quat.from_angax(ov.angax(0.2,0,0,1))), uc, uppermarkerfacecolor='r', lowermarkerfacecolor='r')
pf.add_miller(m, uc)
m = ov.miller(0,0,1)
pf.add_miller(m.rotate(ov.quat.from_angax(ov.angax(0.2,1,0,0))), uc, uppermarkerfacecolor='r', lowermarkerfacecolor='r')
pf.add_miller(m, uc)

pf = ov.stereoproj()
m = ov.miller(1,0,0)
pf.add_miller(m.rotate(ov.quat.from_bunge(ov.bunge(0.2,0,0))), uc, uppermarkerfacecolor='r', lowermarkerfacecolor='r')
pf.add_miller(m, uc)
m = ov.miller(0,0,1)
pf.add_miller(m.rotate(ov.quat.from_bunge(ov.bunge(0,0.2,0))), uc, uppermarkerfacecolor='r', lowermarkerfacecolor='r')
pf.add_miller(m, uc)

