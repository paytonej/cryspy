import ovlib as ov
import matplotlib.pyplot as plt

fnom = 'projtutorial'
ftyp = '.png'
fnum = 0

''' PROJECTION PLOTTING TUTORIAL '''


''' 
To plot stereographic projections using ovlib, you need to import the OV 
library, define a unit cell, and define either some crystal planes or lattice
vectors to plot. 
'''

uc = ov.unitcell() # default with no arguments is cubic
u = [1, 0, 0, 1, 1, 0, 1]
v = [0, 1, 0, 1, 0, 1, 1]
w = [0, 0, 1, 0, 1, 1, 1]
v = ov.lattvec(u,v,w) # define lattice vectors
sg1 = ov.stereoproj() # initialize the stereographic projection
sg1.add_lattvec(v,uc) # add our lattice vectors to the projection

'''
Running this series of commands within Spyder produces the following figure:
'''
fnum += 1
plt.savefig(fnom+str(fnum).zfill(2)+ftyp)


''' We can add the standard grid showing the intersection of {110} and 
{100}-type planes by then executing the following:
'''

sg1.grid_standard(uc)
fnum += 1
plt.savefig(fnom+str(fnum).zfill(2)+ftyp)

'''
Let's now jump to a more complicated example, comparing the locations of the 
same poles in the cubic and monoclinic systems:
'''
plt.close('all')
from pylab import close, figure, subplot, show
close('all')
u = [1, 0, 0, 1, 1, 0, 1]
v = [0, 1, 0, 1, 0, 1, 1]
w = [0, 0, 1, 0, 1, 1, 1]
v = ov.lattvec(u,v,w)
fig1 = figure(1)
sub1 = subplot(121)
sg1 = ov.stereoproj(figure=fig1,subplot=sub1)
uc1 = ov.unitcell() # default is cubit
sg1.add_lattvec(v,uc1)

sub2 = subplot(122)
sg2 = ov.stereoproj(figure=fig1,subplot=sub2)
uc2 = ov.unitcell(beta=45) # monoclinic with a=b=c, alpha=gamma=90, beta=45
sg2.add_lattvec(v,uc2)

'''
This produces a figure with two subplots.
'''
show()
fnum += 1
plt.savefig(fnom+str(fnum).zfill(2)+ftyp)

'''
The angles between the lattice directions can now be compared by adding Wulff 
nets.
'''
sg1.grid_wulffnet()
sg2.grid_wulffnet()
show()
fnum += 1
plt.savefig(fnom+str(fnum).zfill(2)+ftyp)

'''
We can rotate the Wulff net by adding the keyword "rotation." We can also plot
non-standard projections by defining the center, south, and east directions. We
can label these directions using add_eastlabel, add_normallabel, and 
add_northlabel. The exact locations of the labels can be controlled using x and
y keyword arguments.

We can apply additional grid types, including spokes, rings, great circles, and
small circles. Great circles and spokes can cause a crowded appearance around
poles, and this can be reduced using the poleclip keyword argument. Spokes,
rings, great circles, and small circles can also be rotated using the keyword
"rotation."

We can also get a readout of the directions and the planes that correspond to
where our cursor is in the figure using the add_coordinatereadout(uc) command.

Our figures can get correspondingly more complex:
'''

close('all')
from matplotlib.font_manager import FontProperties
u = [1, 0, 0, 1, 1, 0, 1]
v = [0, 1, 0, 1, 0, 1, 1]
w = [0, 0, 1, 0, 1, 1, 1]
v = ov.lattvec(u,v,w)
fig1 = figure(1)
sub1 = subplot(121)
sg1 = ov.stereoproj(figure=fig1, subplot=sub1, center=[1,1,1],south=[-1,-1,2],\
                    east=[-1,1,0])
uc1 = ov.unitcell() # default is cubit
sg1.add_lattvec(v,uc1)
sg1.add_lattveclabels(v,uc1)
sg1.grid_wulffnet(rotation=60)
font1=FontProperties()
font1.set_size('large')
font1.set_style('italic')
sg1.add_eastlabel(fontproperties=font1)
sg1.add_southlabel(fontproperties=font1)
sg1.add_centerlabel(fontproperties=font1)
sg1.add_coordinatereadout(uc1)
sub2 = subplot(122)
sg2 = ov.stereoproj(figure=fig1,subplot=sub2, center=[1,1,1],south=[-1,-1,2],\
                    east=[-1,1,0])
uc2 = ov.unitcell(beta=45) # monoclinic with a=b=c, alpha=gamma=90, beta=45
sg2.add_lattvec(v,uc2,lowermarker='d',lowermarkerfacecolor='g')
sg2.add_lattveclabels(v,uc2)
sg2.grid_rings(degreestep=5)
sg2.grid_spokes(degreestep=5,poleclip=5)
sg2.add_lattveclabels(v,uc2)
font2=FontProperties()
font2.set_family('monospace')
font2.set_size('large')
sg2.add_eastlabel(fontproperties=font2)
sg2.add_southlabel(fontproperties=font2)
sg2.add_centerlabel(fontproperties=font2)
sg2.add_coordinatereadout(uc2)
show()
fnum += 1
plt.savefig(fnom+str(fnum).zfill(2)+ftyp)

'''One final example, showing the plotting of plane traces as well as miller
indices, so that the angles between directions and plane normals can be compared
in a non-cubic crystal:
'''
close('all')
uc = ov.unitcell(a=5,b=2,c=1,alpha=79,beta=90,gamma=108) # default with no 
                                                         # arguments is cubic
u = [3, 0, 0, 1, 1, 0,-2]
v = [2,-1,-5, 0, 1, 1, 1]
w = [1, 6, 1, 3, 1,-1, 3]
vv = ov.lattvec(u,v,w) # define lattice vectors
m = ov.miller(u,v,w) # define lattice vectors
sg = ov.stereoproj() # initialize the stereographic projection
sg.add_lattvec(vv,uc) # add our lattice vectors to the projection
sg.add_lattveclabels(vv,uc)
sg.add_miller(m,uc) # add our lattice vectors to the projection
sg.add_millerlabels(m,uc)
sg.add_trace(m,uc) # add plane traces
sg.grid_greatcircles()
sg.grid_smallcircles()
sg.grid_spokes()
font=FontProperties()
font.set_size('large')
font.set_style('italic')
sg.add_eastlabel(fontproperties=font)
sg.add_southlabel(fontproperties=font)
sg.add_centerlabel(fontproperties=font)
sg.add_coordinatereadout(uc)
fnum += 1
plt.savefig(fnom+str(fnum).zfill(2)+ftyp)
close('all')