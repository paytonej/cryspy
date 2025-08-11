# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 11:18:08 2012

@author: epayton
"""

from ovlib import *
from numpy import array, pi

################################################################################

print "\n check that we can initialize"
print quat()
print tvec()
print bunge()
print rodri()


################################################################################
print '-----------------'
print "\n check that we can construct individual quaternions"

# check that we can construct with a single rotation given as a tuple
q1a=quat(0.6070,-0.7043,0.0634,-0.3627)
print q1a

# check that we can construct with a tuple of lists
q2a=quat([0.6070],[-0.7043],[0.0634],[-0.3627])
print q2a

# check that we can construct with an array
q3a=quat(array([0.6070]),array([-0.7043]),array([0.0634]),array([-0.3627]))
print q3a

################################################################################
print '-----------------'
print "\n check that we can construct multiple quaternions"

# check that we can construct with multiple rotations given as tuples
a=0.6070,0.9688
b=-0.7043,0.1176
c=0.0634,0.0263
d=-0.3627,0.2166
q1b=quat(a,b,c,d)
print q1b

# check that we can construct with multiple values in a tuple of lists
a=[0.6070,0.9688]
b=[-0.7043,0.1176]
c=[0.0634,0.0263]
d=[-0.3627,0.2166]
q2b=quat(a,b,c,d)
print q2b

# check that we can construct with a tuple of arrays
a=array([0.6070,0.9688])
b=array([-0.7043,0.1176])
c=array([0.0634,0.0263])
d=array([-0.3627,0.2166])
q3b=quat(a,b,c,d)
print q3b

################################################################################
print '-----------------'
print "\n check that we can convert quaternions to tvecs"

tq1=tvec.from_quat(q1a)
print tq1

tq2=tvec.from_quat(q1b)
print tq2

################################################################################
print '-----------------'
print "\n check that we can convert tvecs to quats"

qt1=quat.from_tvec(tq1)
print qt1

qt2=quat.from_tvec(tq2)
print qt2

################################################################################
print '-----------------'
print "\n check that we can convert Bunge Euler to quats"

be1=quat.from_bunge(bunge(phi1=pi,PHI=pi/7.,phi2=0.0))

################################################################################
print '-----------------'
print "\n check that we can get symmetries"

r1 = rotationelements('m-3m') # rotation symmetry elements
r2 = pointgroupelements('6/mmm') # point group elements
cs = rotsymm('m3') # rotational symmetry class
print cs

