# -*- coding: utf-8 -*-
"""
ovlib.util: Utilities
===============================================================================
"""

def vecarraycross(b, c):
    ''' cross products of arrays of vectors '''
    import numpy as np
    bx = vecarrayconvert(b[0])
    by = vecarrayconvert(b[1])
    bz = vecarrayconvert(b[2])
    cx = vecarrayconvert(c[0])
    cy = vecarrayconvert(c[1])
    cz = vecarrayconvert(c[2])

    if np.shape(bx) == np.shape(by) == np.shape(bz) == \
       np.shape(cx) == np.shape(cy) == np.shape(cz):

        ax = by * cz - bz * cy
        ay = bz * cx - bx * cz
        az = bx * cy - by * cx

        return [ax, ay, az]
    else:
        print("vecarraycross error: check that the lengths of arguments are equal.") # TODO: check into using warnings package

#-------------------------------------------------------------------------------

def vecarraynorm(a):
    ''' norm of an array of vectors '''
    import numpy as np
    ax = vecarrayconvert(a[0])
    ay = vecarrayconvert(a[1])
    az = vecarrayconvert(a[2])

    if np.shape(ax) == np.shape(ay) == np.shape(az):
        nrm=(ax*ax + ay*ay + az*az)**0.5
    else:
        print("vecarraynorm error: check that the lengths of arguments are equal.") # TODO: check into using warnings package

    return nrm

#-------------------------------------------------------------------------------

def vecarraydot(b,c):
    ''' dot products of an array of vectors  '''
    import numpy as np
    bx = vecarrayconvert(b[0])
    by = vecarrayconvert(b[1])
    bz = vecarrayconvert(b[2])
    cx = vecarrayconvert(c[0])
    cy = vecarrayconvert(c[1])
    cz = vecarrayconvert(c[2])

    if np.shape(bx) == np.shape(by) == np.shape(bz) == \
       np.shape(cx) == np.shape(cy) == np.shape(cz):

        return bx*cx + by*cy + bz*cz

    else:

        print("vecarraydot error: check that the lengths of arguments are equal.") # TODO: check into using warnings package

#-------------------------------------------------------------------------------

def vecarraygcd(a,b):
    ''' Compute greatest common denominator for values in two arrays'''
    import numpy as np

    a = vecarrayconvert(a)
    b = vecarrayconvert(b)

    if np.shape(a) == np.shape(b):

        r = np.zeros( np.size(a) ) # preallocate remainder array

        # if b is smaller than floating point error, we will make gcd=1
        a[ b < np.spacing(1) ] = 1.0
        j = b > np.spacing(1)

        if j.any():
            while j.any():
                r[j] = a[j] % b[j]
                a[j] = b[j]
                b[j] = r[j]
                j[b==0]=False

            return a
        else:
            return a

#-------------------------------------------------------------------------------

def vecarrayconvert(a):
    ''' convert any reasonable datatype into an n x 3 vector array'''
    import numpy as np

    # asanyarray makes an array unless it's already an array
    # squeeze remove unnecessary dimensions
    # atleast_1d makes it so that the arrays are indexable if they are scalars
    a = np.atleast_1d(np.squeeze(np.asanyarray(a)))

    return a

#-------------------------------------------------------------------------------

def stereotrans(xyz, center=[0,0,1], south=[1,0,0], east=[0,1,0]):
    '''
    Perform non-standard stereographic projections following
    Kosel TH, J Mater Sci 19 (1984)
    '''
    import numpy as np
    abc = vecarrayconvert(center)
    nrmabc = 1.0/vecarraynorm(abc)

    uvw = vecarrayconvert(east)
    nrmuvw = 1.0/vecarraynorm(uvw)

    fed = vecarrayconvert(south) # 'def' is reserved in python
    nrmdef = 1.0/vecarraynorm(fed)

    xyz = vecarrayconvert(xyz)
    nrmxyz = 1.0/vecarraynorm(xyz)

    xx = xyz[0]*nrmxyz
    yy = xyz[1]*nrmxyz
    zz = xyz[2]*nrmxyz

    n = np.size(xx)

    aa = np.tile(abc[0]*nrmabc,n)
    bb = np.tile(abc[1]*nrmabc,n)
    cc = np.tile(abc[2]*nrmabc,n)
    dd = np.tile(fed[0]*nrmdef,n)
    ee = np.tile(fed[1]*nrmdef,n)
    ff = np.tile(fed[2]*nrmdef,n)
    uu = np.tile(uvw[0]*nrmuvw,n)
    vv = np.tile(uvw[1]*nrmuvw,n)
    ww = np.tile(uvw[2]*nrmuvw,n)

    cosdl = vecarraydot([xx,yy,zz],[dd,ee,ff])
    cosmu = vecarraydot([xx,yy,zz],[uu,vv,ww])
    cosal = vecarraydot([xx,yy,zz],[aa,bb,cc])

    denom = 1.0/(1.0+np.absolute(cosal))

    xproj =  cosmu * denom
    yproj = -cosdl * denom

    hemis = np.tile('S',n)
    if np.array(n)>1.0:
        hemis[cosal<0.0] = 'N'
    else:
        if cosal<0:
            hemis = np.tile('N',n)

    return xproj, yproj, hemis

#-------------------------------------------------------------------------------

def stereotransstandard(xyz):
    '''
    Perform the standard stereographic transform on cartesian vectors
    following De Graef & McHenry

    Returns projected x, projected y, and the hemisphere of the projection
    '''
    import numpy as np
    x = vecarrayconvert(xyz[0])
    y = vecarrayconvert(xyz[1])
    z = vecarrayconvert(xyz[2])
    nrm = 1.0/vecarraynorm([x,y,z])
    denom = nrm/(1.0+abs(z*nrm))
    xproj =  x * denom
    yproj = -y * denom

    n = np.shape(x)
    hemis = np.tile('S',n)
    if n[0]>1.0:
        hemis[z<0.0] = 'N'
    else:
        if z<0.0:
            hemis = 'N'

    return xproj, yproj, hemis

#-------------------------------------------------------------------------------

def revstereotrans(xyh,center=[0,0,-1],south=[-1,0,0],east=[0,-1,0]):
    '''
    Reverse non-standard stereographic projection
    '''
    import numpy as np
    import numpy.linalg as npla
    abc = vecarrayconvert(center)
    nrmabc = 1.0/vecarraynorm(abc)

    uvw = vecarrayconvert(east)
    nrmuvw = 1.0/vecarraynorm(uvw)

    fed = vecarrayconvert(south) # 'def' is reserved in python
    nrmdef = 1.0/vecarraynorm(fed)

    xproj = vecarrayconvert(xyh[0])
    yproj = vecarrayconvert(xyh[1])
    hemis = np.array(np.squeeze(xyh[2]))

    n = np.shape(xproj)
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)

    aa = abc[0]*nrmabc
    bb = abc[1]*nrmabc
    cc = abc[2]*nrmabc
    dd = fed[0]*nrmdef
    ee = fed[1]*nrmdef
    ff = fed[2]*nrmdef
    uu = uvw[0]*nrmuvw
    vv = uvw[1]*nrmuvw
    ww = uvw[2]*nrmuvw

    R = xproj*xproj + yproj*yproj

    m = np.ones(np.shape(hemis))
    m[hemis=='N'] = -1.0

    cosal = (1.0 - R) / (1.0 + R)
    denom = abs(cosal) + 1.0
    cosmu = xproj * denom
    cosdl = yproj * denom

    # For each case we determine xyz by solving a system of linear equations
    for i in range(0,np.shape(cosal)[0]):
        xyzcoef = np.squeeze(np.array([[dd,ee,ff],[uu,vv,ww],[aa,bb,cc]]))
        xyzsoln = np.array([cosdl[i], cosmu[i], cosal[i]])
        xyz = npla.solve(xyzcoef, xyzsoln)
        x[i]= xyz[0]
        y[i]=-xyz[1]
        z[i]= xyz[2]

    # Correct z for the hemisphere
    z = z * m

    return x,y,z

#-------------------------------------------------------------------------------

def revstereotransstandard(xyh):
    '''
    Performs the reverse stereographic transform on the projected x and y
    coordinates measured from the stereographic projection.
    Reverse of De Graef and McHenry procedure.

    Returns the cartesian x, y, and z values of the normalized vector.
    '''
    import numpy as np

    xproj = vecarrayconvert(xyh[0])
    yproj = vecarrayconvert(xyh[1])
    hemis = np.array(np.squeeze(xyh[2]))

    R = xproj*xproj + yproj*yproj

    m = np.ones(np.shape(hemis))
    m[hemis=='N'] = -1.0
    z = (m - m * R) / (1.0 + R)
    s = m * z + 1.0
    x =  xproj * s
    y = -yproj * s

    return x,y,z

#-------------------------------------------------------------------------------

def eatrans(xyz,center=[0,0,1],south=[0,1,0],east=[1,0,0]):
    '''
    Perform non-standard equal area projections.

    Reference:
    [1] Kosel TH, J Mater Sci 19 (1984)
    '''
    import numpy as np

    abc = vecarrayconvert(center)
    nrmabc = 1.0/vecarraynorm(abc)

    uvw = vecarrayconvert(east)
    nrmuvw = 1.0/vecarraynorm(uvw)

    fed = vecarrayconvert(south) # 'def' is reserved in python
    nrmdef = 1.0/vecarraynorm(fed)

    xyz = vecarrayconvert(xyz)
    nrmxyz = 1.0/vecarraynorm(xyz)

    xx = xyz[0]*nrmxyz
    yy = xyz[1]*nrmxyz
    zz = xyz[2]*nrmxyz

    n = np.shape(xx)

    aa = np.tile(abc[0]*nrmabc,n)
    bb = np.tile(abc[1]*nrmabc,n)
    cc = np.tile(abc[2]*nrmabc,n)
    dd = np.tile(fed[0]*nrmdef,n)
    ee = np.tile(fed[1]*nrmdef,n)
    ff = np.tile(fed[2]*nrmdef,n)
    uu = np.tile(uvw[0]*nrmuvw,n)
    vv = np.tile(uvw[1]*nrmuvw,n)
    ww = np.tile(uvw[2]*nrmuvw,n)

    cosdl = vecarraydot([xx,yy,zz],[dd,ee,ff])
    cosmu = vecarraydot([xx,yy,zz],[uu,vv,ww])
    cosal = vecarraydot([xx,yy,zz],[aa,bb,cc])

    denom = 1.0/np.sqrt(1.0+np.absolute(cosal))

    xproj =  cosmu * denom
    yproj = -cosdl * denom

    hemis = np.tile('S',n)
    if np.array(n)>1.0:
        hemis[cosal<0.0] = 'N'
    else:
        if cosal<0:
            hemis = np.tile('N',n)

    return xproj, yproj, hemis

#-------------------------------------------------------------------------------

def eatransstandard(xyz):
    '''
    Perform the standard stereographic transform on cartesian vectors
    following De Graef & McHenry

    Returns projected x, projected y, and the hemisphere of the projection
    '''
    import numpy as np

    x=vecarrayconvert(xyz[0])
    y=vecarrayconvert(xyz[1])
    z=vecarrayconvert(xyz[2])
    nrm=1.0/vecarraynorm([x,y,z])
    denom=nrm/np.sqrt(1.0+abs(z*nrm))
    xproj =  x * denom
    yproj = -y * denom

    n=np.shape(x)
    hemis = np.tile('S',n)
    if n[0]>1.0:
        hemis[z<0.0] = 'N'
    else:
        if z<0.0:
            hemis = 'N'

    return xproj, yproj, hemis

#-------------------------------------------------------------------------------

def reveatrans(xyh,center=[0,0,1],south=[0,1,0],east=[1,0,0]):
    '''
    Reverse non-standard stereographic projection
    '''
    import numpy as np
    import numpy.linalg as npla

    abc = vecarrayconvert(center)
    nrmabc = 1.0/vecarraynorm(abc)

    uvw = vecarrayconvert(east)
    nrmuvw = 1.0/vecarraynorm(uvw)

    fed = vecarrayconvert(south) # 'def' is reserved in python
    nrmdef = 1.0/vecarraynorm(fed)

    xproj = vecarrayconvert(xyh[0])
    yproj = vecarrayconvert(xyh[1])
    hemis = np.array(np.squeeze(xyh[2]))

    n = np.shape(xproj)
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)

    aa = abc[0]*nrmabc
    bb = abc[1]*nrmabc
    cc = abc[2]*nrmabc
    dd = fed[0]*nrmdef
    ee = fed[1]*nrmdef
    ff = fed[2]*nrmdef
    uu = uvw[0]*nrmuvw
    vv = uvw[1]*nrmuvw
    ww = uvw[2]*nrmuvw

    R = xproj*xproj + yproj*yproj

    m = np.ones(np.shape(hemis))
    m[hemis=='N'] = -1.0

    cosal = (1.0 - R) / (1.0 + R)
    denom = (np.absolute(cosal) + 1.0)**2.0
    cosmu = xproj * denom
    cosdl = yproj * denom

    # For each case we determine xyz by solving a system of linear equations
    for i in range(0, np.shape(cosal)[0]):
        xyzcoef = np.squeeze(np.array([[dd,ee,ff],[uu,vv,ww],[aa,bb,cc]]))
        xyzsoln = np.array([cosdl[i], cosmu[i], cosal[i]])
        xyz = npla.solve(xyzcoef, xyzsoln)
        x[i]= xyz[0]
        y[i]=-xyz[1]
        z[i]= xyz[2]

    # Correct z for the hemisphere
    z = z * m

    return x,y,z

#-------------------------------------------------------------------------------

def reveatransstandard(xyh):
    '''
    Performs the reverse stereographic transform on the projected x and y
    coordinates measured from the stereographic projection.

    Reverse of De Graef and McHenry procedure.

    Returns the cartesian x, y, and z values of the normalized vector.
    '''
    import numpy as np

    xproj = vecarrayconvert(xyh[0])
    yproj = vecarrayconvert(xyh[1])
    hemis = np.array(np.squeeze(xyh[2]))

    R = xproj*xproj + yproj*yproj

    m = np.ones(np.shape(hemis))
    m[hemis=='N'] = -1.0
    z = (m - m * R) / (1.0 + R)
    s = (m * z + 1.0)**2.0
    x =  xproj * s
    y = -yproj * s

    return x,y,z

#-------------------------------------------------------------------------------

def gnomonictrans(xyz,center=[0,0,1],south=[0,1,0],east=[1,0,0]):
    '''
    Perform the gnomonic projection, allowing different projection views
    '''
    import numpy as np

    abc = vecarrayconvert(center)
    nrmabc = 1.0/vecarraynorm(abc)

    uvw = vecarrayconvert(east)
    nrmuvw = 1.0/vecarraynorm(uvw)

    fed = vecarrayconvert(south) # 'def' is reserved in python
    nrmdef = 1.0/vecarraynorm(fed)

    xyz = vecarrayconvert(xyz)
    nrmxyz = 1.0/vecarraynorm(xyz)

    xx = xyz[0]*nrmxyz
    yy = xyz[1]*nrmxyz
    zz = xyz[2]*nrmxyz

    n = np.shape(xx)

    aa = np.tile(abc[0]*nrmabc,n)
    bb = np.tile(abc[1]*nrmabc,n)
    cc = np.tile(abc[2]*nrmabc,n)
    dd = np.tile(fed[0]*nrmdef,n)
    ee = np.tile(fed[1]*nrmdef,n)
    ff = np.tile(fed[2]*nrmdef,n)
    uu = np.tile(uvw[0]*nrmuvw,n)
    vv = np.tile(uvw[1]*nrmuvw,n)
    ww = np.tile(uvw[2]*nrmuvw,n)

    cosdl = vecarraydot([xx,yy,zz],[dd,ee,ff])
    cosmu = vecarraydot([xx,yy,zz],[uu,vv,ww])
    cosal = vecarraydot([xx,yy,zz],[aa,bb,cc])

    denom = 1.0/np.absolute(cosal)

    xproj =  cosmu * denom
    yproj = -cosdl * denom

    hemis = np.tile('S',n)
    if np.array(n)>1.0:
        hemis[cosal<0.0] = 'N'
    else:
        if cosal<0:
            hemis = np.tile('N',n)

    return xproj, yproj, hemis

#-------------------------------------------------------------------------------

def gnomonictransstandard(xyz):
    '''
    Perform the standard gnomonic projection

    Returns projected x, projected y, and the hemisphere of the projection
    '''
    import numpy as np

    x=vecarrayconvert(xyz[0])
    y=vecarrayconvert(xyz[1])
    z=vecarrayconvert(xyz[2])
    denom=1.0/np.absolute(z)
    xproj =  x * denom
    yproj = -y * denom

    n=np.shape(x)
    hemis = np.tile('S',n)
    if n[0]>1.0:
        hemis[z<0.0] = 'N'
    else:
        if z<0.0:
            hemis = 'N'

    return xproj, yproj, hemis

#-------------------------------------------------------------------------------

def revgnomonictrans(xyh,center=[0,0,1],south=[0,1,0],east=[1,0,0]):
    '''
    reverse gnomonic projection allowing different projection views
    '''
    import numpy as np
    import numpy.linalg as npla

    abc = vecarrayconvert(center)
    nrmabc = 1.0/vecarraynorm(abc)

    uvw = vecarrayconvert(east)
    nrmuvw = 1.0/vecarraynorm(uvw)

    fed = vecarrayconvert(south) # 'def' is reserved in python
    nrmdef = 1.0/vecarraynorm(fed)

    xproj = vecarrayconvert(xyh[0])
    yproj = vecarrayconvert(xyh[1])
    hemis = np.array(np.squeeze(xyh[2]))

    n = np.shape(xproj)
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)

    aa = abc[0]*nrmabc
    bb = abc[1]*nrmabc
    cc = abc[2]*nrmabc
    dd = fed[0]*nrmdef
    ee = fed[1]*nrmdef
    ff = fed[2]*nrmdef
    uu = uvw[0]*nrmuvw
    vv = uvw[1]*nrmuvw
    ww = uvw[2]*nrmuvw

    R = xproj*xproj + yproj*yproj

    m = np.ones(np.shape(hemis))
    m[hemis=='N'] = -1.0

    cosal = (1.0 - R) / (1.0 + R)
    denom = np.absolute(cosal)
    cosmu = xproj * denom
    cosdl = yproj * denom

    # For each case we determine xyz by solving a system of linear equations
    for i in range(0,np.shape(cosal)[0]):
        xyzcoef = np.squeeze(np.array([[dd,ee,ff],[uu,vv,ww],[aa,bb,cc]]))
        xyzsoln = np.array([cosdl[i], cosmu[i], cosal[i]])
        xyz = npla.solve(xyzcoef, xyzsoln)
        x[i]= xyz[0]
        y[i]=-xyz[1]
        z[i]= xyz[2]

    # Correct z for the hemisphere
    z = z * m

    return x,y,z

#-------------------------------------------------------------------------------

def revgnomonictransstandard(xyh):
    '''
    performs the reverse standard gnomonic projection
    '''
    import numpy as np

    xproj = vecarrayconvert(xyh[0])
    yproj = vecarrayconvert(xyh[1])
    hemis = np.array(np.squeeze(xyh[2]))

    R = xproj*xproj + yproj*yproj

    m = np.ones(np.shape(hemis))
    m[hemis=='N'] = -1.0
    z = (m - m * R) / (1.0 + R)
    s =  m * z
    x =  xproj * s
    y = -yproj * s

    return x,y,z

#-------------------------------------------------------------------------------

def xtaldot(p1=1,p2=0,p3=0,
            g11=1,g12=0,g13=0,g21=0,g22=1,g23=0,g31=0,g32=0,g33=1,
            q1=1,q2=0,q3=0):
    ''' implements the dot product within the crystallographic reference frame:
        p*G\q where p and q are vectors and G is a metric matrix '''

    p1  = vecarrayconvert(p1)
    p2  = vecarrayconvert(p2)
    p3  = vecarrayconvert(p3)
    q1  = vecarrayconvert(q1)
    q2  = vecarrayconvert(q2)
    q3  = vecarrayconvert(q3)
    g11 = vecarrayconvert(g11)
    g12 = vecarrayconvert(g12)
    g13 = vecarrayconvert(g13)
    g21 = vecarrayconvert(g21)
    g22 = vecarrayconvert(g22)
    g23 = vecarrayconvert(g23)
    g31 = vecarrayconvert(g31)
    g32 = vecarrayconvert(g32)
    g33 = vecarrayconvert(g33)

    return (g11*q1+g12*q2+g13*q3)*p1+ \
         (g21*q1+g22*q2+g23*q3)*p2+ \
         (g31*q1+g32*q2+g33*q3)*p3

#-------------------------------------------------------------------------------

def xtalangle(p1=1,p2=0,p3=0,q1=1,q2=0,q3=0,
            g11=1,g12=0,g13=0,g21=0,g22=1,g23=0,g31=0,g32=0,g33=1):
    ''' compute the angle between two directions in a crystal'''

    import numpy as np

    p1  = vecarrayconvert(p1)
    p2  = vecarrayconvert(p2)
    p3  = vecarrayconvert(p3)
    q1  = vecarrayconvert(q1)
    q2  = vecarrayconvert(q2)
    q3  = vecarrayconvert(q3)
    g11 = vecarrayconvert(g11)
    g12 = vecarrayconvert(g12)
    g13 = vecarrayconvert(g13)
    g21 = vecarrayconvert(g21)
    g22 = vecarrayconvert(g22)
    g23 = vecarrayconvert(g23)
    g31 = vecarrayconvert(g31)
    g32 = vecarrayconvert(g32)
    g33 = vecarrayconvert(g33)

    nrm=1.0/(vecarraynorm([p1,p2,p3])*vecarraynorm([q1,q2,q3]))

    return np.arccos(nrm *
                xtaldot(p1,p2,p3,g11,g12,g13,g21,g22,g23,g31,g32,g33,q1,q2,q3))

#-------------------------------------------------------------------------------

def rationalize(v, maxval=9):
    ''' produce rational indices from fractions '''
    import numpy as np
    import copy

    v0 = vecarrayconvert(v[0])
    v1 = vecarrayconvert(v[1])
    v2 = vecarrayconvert(v[2])
    nrm=1.0/vecarraynorm([v0,v1,v2])
    v0 = v0 * nrm
    v1 = v1 * nrm
    v2 = v2 * nrm

    if np.shape(v0)==np.shape(v1)==np.shape(v2):
        if np.size(v0)==1:
            v0=np.array([v0])
            v1=np.array([v1])
            v2=np.array([v2])

        n=np.size(v[0])
        vi=np.zeros([n,maxval])
        vj=copy.copy(vi); vk=copy.copy(vi); vq=copy.copy(vi)
        for i in range(1,maxval+1):
            vx = np.around(np.float64(i) * v0)
            vy = np.around(np.float64(i) * v1)
            vz = np.around(np.float64(i) * v2)
            tmpx = np.absolute(vx); tmpy = np.absolute(vy); tmpz = np.absolute(vz)

            # calculate the greatest common divisor between tmpx, tmpy, and tmpz
            # we have to do this manually because fractions.gcd doesn't work on
            # arrays and numpy doesn't yet have a gcd function implemented
            div = 1.0/vecarraygcd(vecarraygcd(tmpx,tmpy),tmpz)

            # multipy the irrational indices by the greatest common divisor
            vi[:,i-1] = vx * div
            vj[:,i-1] = vy * div
            vk[:,i-1] = vz * div
            nrm=1.0/vecarraynorm([vi[:,i-1],vj[:,i-1],vk[:,i-1]])
            vq[:,i-1] = np.arccos(np.amin([
                          np.tile(1.0,n),
                          vecarraydot([ vi[:,i-1]*nrm,vj[:,i-1]*nrm,vk[:,i-1]*nrm],
                                   [v[0],v[1],v[2]])
                                     ],axis=0))

        # extract the best match rational values
        loc=np.argmin(vq.T, axis=0)
        vi=vi[range(0, n),loc]
        vj=vj[range(0, n),loc]
        vk=vk[range(0, n),loc]
        vq=vq[range(0, n),loc]

        return vi, vj, vk,vq

    else:
        print("rationalize error:"+ \
                "check that the lengths of arguments are equal.") # TODO: check into using warnings package

#------------------------------------------------------------------------------

def rand_vonMisesFisherM(n, kappa=0, mu=[1.0, 0.0, 0.0, 0.0]):
    """ random number generation from von Mises-Fisher matrix distribution

    Return n samples of random unit directions centered around mu with
    dispersion parameter kappa

    Parameters
    ----------
    n : int
        number of samples

    kappa : float > 0
        concentration parameter

    mu : iterable of floats
        central vector; length determines dimensionality m

    Returns
    -------
    x : n x m numpy array
        rows correspond to random unit vectors from distribution

    Notes
    -----
    This is a python translation of Sungkyu Jung's matlab code published on
    his website, version dated 3 Feb 2010 [1]_. It uses the modified Ulrich's
    algorithm from Wood [2]_.

    References
    ----------
    .. [1] S. Jung, M. Foskey, J.S. Marron, "Principal Arc Analysis on Direct
           Product Manifolds," Annnals of Applied Statistics (2011).
    .. [2] A.T.A. Wood, "Simulation of the von Mises Fisher distribution," Commun.
           Statist. 23 (1994).
    """
    import numpy as np
    import scipy.stats as stats

    # Convert mu to a 2d numpy array
    mu = np.atleast_2d(np.squeeze(np.asarray(mu).ravel()))

    # the dimensionality is always given by the size of mu
    m = mu.size

    b = (-2.0 * kappa + np.sqrt(4.0 * kappa**2.0 + (m - 1.0)**2.0)) / (m - 1.0)
    x0 = (1.0 - b) / (1.0 + b)
    c = kappa * x0 + (m - 1.0)*np.log(1.0 - x0**2.0)

    # steps 1 & 2 from [2]
    nnow = n
    ww = []
    while True:
        ntrial = np.amax([np.around(nnow * 1.2), nnow + 10.0])
        z = stats.beta.rvs((m - 1.0) / 2.0, (m - 1.0) / 2.0, size=ntrial)
        u = np.random.rand(ntrial)
        w = (1.0 - (1.0 + b) * z) / (1.0 - (1.0 - b) * z)

        indicator = kappa * w + (m - 1.0) * np.log(1.0 - x0 * w) - c >= np.log(u)
        if np.sum(indicator) >= nnow:
            w1 = w[indicator]
            ww = np.hstack([ww, w1[0:nnow]])
            break
        else:
            ww = np.hstack([ww, w[indicator]])
            nnow = nnow - np.sum[indicator]

    # step 3 from [2]: generate n uniformly distributed m dimensional random
    # directions, using the logic: "directions of normal distribution are
    # uniform on the sphere."
    v = np.zeros([m - 1, n])
    nr = stats.norm.rvs(1.0, size=[m - 1, n])
    for i in range(0, n):
        while True:
            ni = np.dot(nr[:, i], nr[:, i]) # length of ith vector
            # exclude too small values to avoid numerical discretization
            if ni < np.sqrt(np.spacing(np.float64(1))):
                # repeat randomization
                nr[:, i] = stats.norm.rvs(1.0, size=[m - 1, 1])
            else:
                v[:, i] = nr[:, i] / np.sqrt(ni)
                break

    x = np.vstack([np.tile(np.sqrt(1.0 - ww**2.0), [m - 1, 1]) * v,
                   np.atleast_2d(ww)])


    # Get the rotation matrix that rotates the data to be centered at mu
    d = mu.size
    a = np.zeros(d)
    a[-1] = 1.0
    a = np.atleast_2d(a)
    ab = np.dot(a, mu.T)
    alpha = np.arccos(ab)
    ii = np.eye(d)

    if   np.abs(ab - 1) < 1e-15:
        rot =  ii
    elif np.abs(ab + 1) < 1e-15:
        rot = -ii
    else:
        c = mu - a * ab
        c = c / np.linalg.norm(c)
        aa = np.dot(a.T, c) - np.dot(c.T, a)
        rot = ii + np.sin(alpha)*aa + (np.cos(alpha) - 1.0)*(np.dot(a.T, a) + np.dot(c.T, c))

    return np.dot(rot.T, x).T

#------------------------------------------------------------------------------

class helper(object):
    """ Provides system-specific information to cryspy
    """
    def __init__(self):

        import tempfile, os, platform, warnings, inspect

        # get platform-dependent temporary directory
        tmp_path = tempfile.gettempdir()

        # Determine the path to cryspy
        cryspy_path = os.sep.join(inspect.getfile(inspect.currentframe()).split(os.sep)[0:-2])

        # Get the bits used by the platform
        bits, linkage = platform.architecture()
        bits = bits[0:2]

        # determine the operating system platform
        plat = platform.system().lower()

        # parse the platform bits and OS
        result = []
        extension = ''
        if plat == 'windows':
            result = 'win' + bits
            extension = '.exe'
        elif plat == 'darwin' and bits == '64':
            result = 'maci' + bits
        elif plat == 'darwin' and bits == '32':
            result = 'maci'
        elif (plat == 'linux' or plat == 'linux2') and bits == '32':
            result = 'glnx86'
        elif (plat == 'linux' or plat == 'linux2') and bits == '64':
            result = 'glnxa64'
        else:
            warnings.warn('Required c binary not compiled for your system.')

        # Create the path to our platform-dependent binary executables
        prgpth = fullfile([cryspy_path, 'c', 'bin', result])

        # Store this info in our helper object
        self.tmpdir = tmp_path
        self.arch   = result
        self.ext    = extension
        self.prgpth = prgpth

#------------------------------------------------------------------------------

def uniformBungeGrid(delta):
    '''
    Produce equal volume sampling over orientation space in Bunge Euler angles

    Parameters
    ----------
    delta :: int

    Returns
    -------
    b :: Bunge Euler angle object

    References
    ----------
    .. [1] Kalidendi et al, Acta Mater 57 (2009) p3916. doi:10.1016/j.actamat.2009.04.055

    '''
    import numpy as np
    import cryspy.rot as rot

    phi1 = np.linspace(0.0, 2.0 * np.pi, delta)
    PHI  = np.cos(phi1)
    phi2 = np.copy(phi1)
    b = rot.bunge(phi1, PHI, phi2)

    return b

#------------------------------------------------------------------------------

# TODO:    def genCSL()