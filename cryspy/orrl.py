# -*- coding: utf-8 -*-
'''
ovlib.orrl: Orientation relationship characterization
===============================================================================
'''

def namedOR(name):
    """ Returns ksi values for named orientation relationships

    Parameters
    ----------
    name : {'ks', 'nw', 'bain'}
        Orientation relationship name

    Returns
    -------
    ksi_values : 1x3 numpy array of floats
        ksi values in degrees

    Notes
    -----
    TODO: Add plane parallel Greninger-Troiano, Kelly, etc.
    TODO: Allow 'Kurdjumov-Sachs' as well as 'ks', etc.
    """
    import numpy as np

    if isinstance(name, str):

        if name.lower() == 'ks':

            s6 = np.sqrt(6.0)
            s3 = np.sqrt(3.0)
            ksi1 = np.arccos((s6 + 1.0) / (2.0 * s3))
            ksi2 = np.arccos((s6 + 18.0) / (12.0 * s3))
            ksi3 = np.arccos((s6 + 12.0) / (6.0 * s6))
            ksi = np.array([ksi1, ksi2, ksi3])
            del s6, s3, ksi1, ksi2, ksi3

        elif name.lower() == 'nw':

            s6 = np.sqrt(6)
            s2 = np.sqrt(2)
            ksi0 = np.arccos((s2 + 1.0) / s6)
            ksi = np.array([0.0, ksi0, ksi0])

        elif name.lower() == 'bain':

            ksi = np.array([0.0, 0.0, 0.0])

        else:

            print('namedOR: Unrecognized named OR') # TODO: use warnings package

    else:

        print('namedOR requires a string input. Returning Bain.') # TODO: use warnings package
        ksi = np.array([0.0, 0.0, 0.0])

    return ksi * 180.0/np.pi

def getRepresentativeVariant(ksi_values):
    """ Returns quaternion of representative variant

    Parameters
    ----------
    ksi_values :

    Returns
    -------
    psiq : quat class
    """
    import numpy as np
    from cryspy.rot import quat

    if isinstance(ksi_values, str):
        # Import named OR and convert to radians
        orrl = namedOR(ksi_values) * np.pi / 180.0
    else:
        # Convert OR specification into radians
        orrl = ksi_values * np.pi / 180.0

    csum = np.cos(orrl[0])+np.cos(orrl[1])+np.cos(orrl[2])
    x0 = 0.5 * np.sqrt(csum + 1.0)
    x = 0.5 * np.sqrt(2. * np.cos(orrl) - csum + 1.0)

    ksiq = quat(x0, x[0], x[1], x[2])

    # Calculate quaternion for normalized Bain correspondence.

    qw = 0.5 * np.sqrt(2.0 + np.sqrt(2.0))
    qz = (2.0 / np.sqrt(2.0)) / (4.0 * qw)
    gammaq = quat(qw, 0.0, 0.0, qz);

    psiq = ksiq * gammaq

    return psiq

# --------------------------------------------------------------------------

def generateVariants(ksi_values):
    """ Generate variants from Kurdjumov-Sachs angles

    Returns matrices of an orientation relationship specified in Kurjumov-Sachs
    angles.

    Parameters
    ----------
    ksi_values : length 3 iterable OR {'KS', 'NW', 'Bain'}

    Returns
    -------
    vv : rmat object
        rotation matrices corresponding to variants

    """
    import numpy as np
    import cryspy.rot as rot
    import cryspy.util as util
    import righthand as rh
    import warnings

    if isinstance(ksi_values, str):

        ksi = namedOR(ksi_values)

    elif np.shape(ksi_values) == (3,):

        ksi = np.asanyarray(ksi_values)

    else:

        warnings.warn('ksi values must be a numeric iterable of length 3.')

    # convert ksi radians to rotation matrices

    ksi = ksi * np.pi / 180.0

    mb = np.zeros([2, 9])

    mb[0, 0] = np.cos(ksi[0])
    mb[0, 4] = np.cos(ksi[1])
    mb[0, 8] = np.cos(ksi[2])

    costh = 0.5 * (np.sum(np.cos(ksi)) - 1.0) # sum(cos(ksi)) is the matrix trace
    mosth = 1.0 - costh
    sinth = np.sqrt(1.0 - costh**2.0)


    r1 = np.sqrt((mb[0, 0] - costh) / mosth)
    r2 = np.sqrt((mb[0, 4] - costh) / mosth)
    r3 = np.sqrt((mb[0, 8] - costh) / mosth)
    del costh

    r1r2 = r1 * r2 * mosth
    r1r3 = r1 * r3 * mosth
    r2r3 = r2 * r3 * mosth
    r3st = r3 * sinth
    r2st = r2 * sinth
    r1st = r1 * sinth
    del r1, r2, r3, mosth, sinth

    mb[0, 5] = r2r3 - r1st
    mb[0, 7] = r2r3 + r1st
    mb[1, :] = mb[0, :]

    mb[0, 1] = -r1r2 + r3st
    mb[0, 2] = -r1r3 - r2st
    mb[0, 3] = -r1r2 - r3st
    mb[0, 6] = -r1r3 + r2st
    del r1r2, r1r3, r2r3, r3st, r2st, r1st

    mb[1, 1] = -mb[0, 1]
    mb[1, 2] = -mb[0, 2]
    mb[1, 3] = -mb[0, 3]
    mb[1, 6] = -mb[0, 6]
    # mb[0] is the 'positive' solution; mb[1] is the 'negative' solution

    # create Bain correspondence matrices
    bb = np.zeros([12, 9])
    bb[ 0, :] = [ 1.0, -1.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  1.0]
    bb[ 1, :] = [ 0.0,  1.0, -1.0,  0.0,  1.0,  1.0,  1.0,  0.0,  0.0]
    bb[ 2, :] = [-1.0,  0.0,  1.0,  1.0,  0.0,  1.0,  0.0,  1.0,  0.0]
    bb[ 3, :] = [ 0.0,  1.0,  1.0,  0.0, -1.0,  1.0,  1.0,  0.0,  0.0]
    bb[ 4, :] = [-1.0, -1.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  1.0]
    bb[ 5, :] = [ 1.0,  0.0, -1.0,  1.0,  0.0,  1.0,  0.0, -1.0,  0.0]
    bb[ 6, :] = [ 1.0,  1.0,  0.0, -1.0,  1.0,  0.0,  0.0,  0.0,  1.0]
    bb[ 7, :] = [-1.0,  0.0, -1.0, -1.0,  0.0,  1.0,  0.0,  1.0,  0.0]
    bb[ 8, :] = [ 0.0, -1.0,  1.0,  0.0,  1.0,  1.0, -1.0,  0.0,  0.0]
    bb[ 9, :] = [ 1.0,  0.0,  1.0,  1.0,  0.0, -1.0,  0.0,  1.0,  0.0]
    bb[10, :] = [ 0.0, -1.0, -1.0,  0.0,  1.0, -1.0,  1.0,  0.0,  0.0]
    bb[11, :] = [-1.0,  1.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0, -1.0]

    # normalize correspondence matrices
    bb = rot.rmat.from_array(bb / util.vecarraynorm(bb))
    mb = rot.rmat.from_array(mb)

    # produce variants
    vv = np.zeros([24, 9])
    tmp = mb[0] * bb
    vv[np.arange(0, 24, 2), :] = tmp.to_array()
    tmp = mb[1] * bb
    vv[np.arange(1, 24, 2), :] = tmp.to_array()

    # reduce redundancies, if they exist (as they do, for example, in NW)
    vv, ia, ic = rh.uniquerows(rh.sigdec(vv, 7))
    del ia, ic

    return rot.rmat.from_array(vv)

#------------------------------------------------------------------------------

def _bg_ksi1(bincenters):
    """ returns ksi1 background

    Parameters
    ----------
    bincenters : numpy array
        centers of bins for ksi1 background

    Returns
    -------
    bg : numpy array
        background number fractions for ksi1

    Notes
    -----
    Background is determined by cubic spline interpolation from a large
    numerical simulation containing 17156424 ksi values measured from
    random orientations.

    """
    import numpy as np
    import scipy.interpolate as interp

    d = np.array([ [  0.00000000e+00,   0.00000000e+00],
                   [  1.50000000e-01,   8.04400000e+03],
                   [  4.50000000e-01,   2.45720000e+04],
                   [  7.50000000e-01,   4.06200000e+04],
                   [  1.05000000e+00,   5.61680000e+04],
                   [  1.35000000e+00,   7.13330000e+04],
                   [  1.65000000e+00,   8.80860000e+04],
                   [  1.95000000e+00,   1.02050000e+05],
                   [  2.25000000e+00,   1.17010000e+05],
                   [  2.55000000e+00,   1.31760000e+05],
                   [  2.85000000e+00,   1.47720000e+05],
                   [  3.15000000e+00,   1.60180000e+05],
                   [  3.45000000e+00,   1.75440000e+05],
                   [  3.75000000e+00,   1.88830000e+05],
                   [  4.05000000e+00,   2.03510000e+05],
                   [  4.35000000e+00,   2.16520000e+05],
                   [  4.65000000e+00,   2.29790000e+05],
                   [  4.95000000e+00,   2.43770000e+05],
                   [  5.25000000e+00,   2.56080000e+05],
                   [  5.55000000e+00,   2.68910000e+05],
                   [  5.85000000e+00,   2.81760000e+05],
                   [  6.15000000e+00,   2.93690000e+05],
                   [  6.45000000e+00,   3.06460000e+05],
                   [  6.75000000e+00,   3.16310000e+05],
                   [  7.05000000e+00,   3.27680000e+05],
                   [  7.35000000e+00,   3.40790000e+05],
                   [  7.65000000e+00,   3.50150000e+05],
                   [  7.95000000e+00,   3.63730000e+05],
                   [  8.25000000e+00,   3.74700000e+05],
                   [  8.55000000e+00,   3.83250000e+05],
                   [  8.85000000e+00,   3.93360000e+05],
                   [  9.15000000e+00,   4.01330000e+05],
                   [  9.45000000e+00,   4.08800000e+05],
                   [  9.75000000e+00,   4.10170000e+05],
                   [  1.00500000e+01,   4.16810000e+05],
                   [  1.03500000e+01,   4.22100000e+05],
                   [  1.06500000e+01,   4.24130000e+05],
                   [  1.09500000e+01,   4.25640000e+05],
                   [  1.12500000e+01,   4.26720000e+05],
                   [  1.15500000e+01,   4.26860000e+05],
                   [  1.18500000e+01,   4.25850000e+05],
                   [  1.21500000e+01,   4.25010000e+05],
                   [  1.24500000e+01,   4.19670000e+05],
                   [  1.27500000e+01,   4.14030000e+05],
                   [  1.30500000e+01,   4.06510000e+05],
                   [  1.33500000e+01,   3.97900000e+05],
                   [  1.36500000e+01,   3.85480000e+05],
                   [  1.39500000e+01,   3.72070000e+05],
                   [  1.42500000e+01,   3.55970000e+05],
                   [  1.45500000e+01,   3.36740000e+05],
                   [  1.48500000e+01,   3.11690000e+05],
                   [  1.51500000e+01,   2.79820000e+05],
                   [  1.54500000e+01,   2.57130000e+05],
                   [  1.57500000e+01,   2.35370000e+05],
                   [  1.60500000e+01,   2.18230000e+05],
                   [  1.63500000e+01,   2.01140000e+05],
                   [  1.66500000e+01,   1.84260000e+05],
                   [  1.69500000e+01,   1.68320000e+05],
                   [  1.72500000e+01,   1.54590000e+05],
                   [  1.75500000e+01,   1.39760000e+05],
                   [  1.78500000e+01,   1.27400000e+05],
                   [  1.81500000e+01,   1.14020000e+05],
                   [  1.84500000e+01,   1.03180000e+05],
                   [  1.87500000e+01,   9.14750000e+04],
                   [  1.90500000e+01,   8.00980000e+04],
                   [  1.93500000e+01,   7.07470000e+04],
                   [  1.96500000e+01,   6.10580000e+04],
                   [  1.99500000e+01,   5.08750000e+04],
                   [  2.02500000e+01,   4.25510000e+04],
                   [  2.05500000e+01,   3.39610000e+04],
                   [  2.08500000e+01,   2.52270000e+04],
                   [  2.11500000e+01,   1.73040000e+04],
                   [  2.14500000e+01,   1.12930000e+04],
                   [  2.17500000e+01,   7.17200000e+03],
                   [  2.20500000e+01,   3.76000000e+03],
                   [  2.23500000e+01,   1.55200000e+03],
                   [  2.26500000e+01,   3.94000000e+02],
                   [  2.29500000e+01,   8.00000000e+00],
                   [  2.32500000e+01,   0.00000000e+00],
                   [  2.35500000e+01,   0.00000000e+00],
                   [  2.38500000e+01,   0.00000000e+00],
                   [  2.41500000e+01,   0.00000000e+00],
                   [  2.44500000e+01,   0.00000000e+00],
                   [  2.47500000e+01,   0.00000000e+00],
                   [  2.50500000e+01,   0.00000000e+00],
                   [  2.53500000e+01,   0.00000000e+00],
                   [  2.56500000e+01,   0.00000000e+00],
                   [  2.59500000e+01,   0.00000000e+00],
                   [  2.62500000e+01,   0.00000000e+00],
                   [  2.65500000e+01,   0.00000000e+00],
                   [  2.68500000e+01,   0.00000000e+00],
                   [  2.71500000e+01,   0.00000000e+00],
                   [  2.74500000e+01,   0.00000000e+00],
                   [  2.77500000e+01,   0.00000000e+00],
                   [  2.80500000e+01,   0.00000000e+00],
                   [  2.83500000e+01,   0.00000000e+00],
                   [  2.86500000e+01,   0.00000000e+00],
                   [  2.89500000e+01,   0.00000000e+00],
                   [  2.92500000e+01,   0.00000000e+00],
                   [  2.95500000e+01,   0.00000000e+00],
                   [  2.98500000e+01,   0.00000000e+00],
                   [  3.01500000e+01,   0.00000000e+00],
                   [  3.04500000e+01,   0.00000000e+00],
                   [  3.07500000e+01,   0.00000000e+00],
                   [  3.10500000e+01,   0.00000000e+00],
                   [  3.13500000e+01,   0.00000000e+00],
                   [  3.16500000e+01,   0.00000000e+00],
                   [  3.19500000e+01,   0.00000000e+00],
                   [  3.22500000e+01,   0.00000000e+00],
                   [  3.25500000e+01,   0.00000000e+00],
                   [  3.28500000e+01,   0.00000000e+00],
                   [  3.31500000e+01,   0.00000000e+00],
                   [  3.34500000e+01,   0.00000000e+00],
                   [  3.37500000e+01,   0.00000000e+00],
                   [  3.40500000e+01,   0.00000000e+00],
                   [  3.43500000e+01,   0.00000000e+00],
                   [  3.46500000e+01,   0.00000000e+00]])

    minval = np.sqrt(np.spacing(np.float32(1)))
    cj = interp.UnivariateSpline(d[:,0],
                                 d[:,1]+minval,
                                 k=1, s=0.5)

    # Interpolate the spline for our bin array into a background for our histogram
    bg = cj(bincenters)
    bg[bg < minval] = minval # do not allow division by zero

    return bg

def _bg_ksi2(bincenters):
    """ returns ksi1 background

    Parameters
    ----------
    bincenters : numpy array
        centers of bins for ksi1 background

    Returns
    -------
    bg : numpy array
        background number fractions for ksi1

    Notes
    -----
    Background is determined by cubic spline interpolation from a large
    numerical simulation containing 17156424 ksi values measured from
    random orientations. The bin spacing in this histogram is 0.3 degrees.

    """
    import numpy as np
    import scipy.interpolate as interp

    d = np.array([ [  0.00000000e+00,   0.00000000e+00],
                   [  1.50000000e-01,   3.80000000e+01],
                   [  4.50000000e-01,   2.53000000e+02],
                   [  7.50000000e-01,   6.37000000e+02],
                   [  1.05000000e+00,   1.30400000e+03],
                   [  1.35000000e+00,   2.14800000e+03],
                   [  1.65000000e+00,   3.15500000e+03],
                   [  1.95000000e+00,   4.51300000e+03],
                   [  2.25000000e+00,   6.03800000e+03],
                   [  2.55000000e+00,   7.55500000e+03],
                   [  2.85000000e+00,   9.49000000e+03],
                   [  3.15000000e+00,   1.15950000e+04],
                   [  3.45000000e+00,   1.37350000e+04],
                   [  3.75000000e+00,   1.64950000e+04],
                   [  4.05000000e+00,   1.88700000e+04],
                   [  4.35000000e+00,   2.18850000e+04],
                   [  4.65000000e+00,   2.54230000e+04],
                   [  4.95000000e+00,   2.88190000e+04],
                   [  5.25000000e+00,   3.21320000e+04],
                   [  5.55000000e+00,   3.60380000e+04],
                   [  5.85000000e+00,   3.99840000e+04],
                   [  6.15000000e+00,   4.36950000e+04],
                   [  6.45000000e+00,   4.81920000e+04],
                   [  6.75000000e+00,   5.32350000e+04],
                   [  7.05000000e+00,   5.78030000e+04],
                   [  7.35000000e+00,   6.26690000e+04],
                   [  7.65000000e+00,   6.77040000e+04],
                   [  7.95000000e+00,   7.24380000e+04],
                   [  8.25000000e+00,   7.90320000e+04],
                   [  8.55000000e+00,   8.46900000e+04],
                   [  8.85000000e+00,   9.05220000e+04],
                   [  9.15000000e+00,   9.71920000e+04],
                   [  9.45000000e+00,   1.03820000e+05],
                   [  9.75000000e+00,   1.10760000e+05],
                   [  1.00500000e+01,   1.17170000e+05],
                   [  1.03500000e+01,   1.23890000e+05],
                   [  1.06500000e+01,   1.31260000e+05],
                   [  1.09500000e+01,   1.38800000e+05],
                   [  1.12500000e+01,   1.47670000e+05],
                   [  1.15500000e+01,   1.54370000e+05],
                   [  1.18500000e+01,   1.62540000e+05],
                   [  1.21500000e+01,   1.71450000e+05],
                   [  1.24500000e+01,   1.80180000e+05],
                   [  1.27500000e+01,   1.88610000e+05],
                   [  1.30500000e+01,   1.96360000e+05],
                   [  1.33500000e+01,   2.06070000e+05],
                   [  1.36500000e+01,   2.15760000e+05],
                   [  1.39500000e+01,   2.25380000e+05],
                   [  1.42500000e+01,   2.34630000e+05],
                   [  1.45500000e+01,   2.45320000e+05],
                   [  1.48500000e+01,   2.56620000e+05],
                   [  1.51500000e+01,   2.62000000e+05],
                   [  1.54500000e+01,   2.64870000e+05],
                   [  1.57500000e+01,   2.69010000e+05],
                   [  1.60500000e+01,   2.71560000e+05],
                   [  1.63500000e+01,   2.74630000e+05],
                   [  1.66500000e+01,   2.78040000e+05],
                   [  1.69500000e+01,   2.80570000e+05],
                   [  1.72500000e+01,   2.82790000e+05],
                   [  1.75500000e+01,   2.85960000e+05],
                   [  1.78500000e+01,   2.88360000e+05],
                   [  1.81500000e+01,   2.89220000e+05],
                   [  1.84500000e+01,   2.92980000e+05],
                   [  1.87500000e+01,   2.93430000e+05],
                   [  1.90500000e+01,   2.96900000e+05],
                   [  1.93500000e+01,   2.97160000e+05],
                   [  1.96500000e+01,   2.98050000e+05],
                   [  1.99500000e+01,   3.00280000e+05],
                   [  2.02500000e+01,   3.02770000e+05],
                   [  2.05500000e+01,   3.01800000e+05],
                   [  2.08500000e+01,   3.04190000e+05],
                   [  2.11500000e+01,   3.05480000e+05],
                   [  2.14500000e+01,   3.02770000e+05],
                   [  2.17500000e+01,   3.01660000e+05],
                   [  2.20500000e+01,   2.96840000e+05],
                   [  2.23500000e+01,   2.91940000e+05],
                   [  2.26500000e+01,   2.87030000e+05],
                   [  2.29500000e+01,   2.79490000e+05],
                   [  2.32500000e+01,   2.74410000e+05],
                   [  2.35500000e+01,   2.68910000e+05],
                   [  2.38500000e+01,   2.60700000e+05],
                   [  2.41500000e+01,   2.54640000e+05],
                   [  2.44500000e+01,   2.47580000e+05],
                   [  2.47500000e+01,   2.42160000e+05],
                   [  2.50500000e+01,   2.36140000e+05],
                   [  2.53500000e+01,   2.29720000e+05],
                   [  2.56500000e+01,   2.26630000e+05],
                   [  2.59500000e+01,   2.17650000e+05],
                   [  2.62500000e+01,   2.12590000e+05],
                   [  2.65500000e+01,   2.06500000e+05],
                   [  2.68500000e+01,   2.01920000e+05],
                   [  2.71500000e+01,   1.95560000e+05],
                   [  2.74500000e+01,   1.88290000e+05],
                   [  2.77500000e+01,   1.84200000e+05],
                   [  2.80500000e+01,   1.78390000e+05],
                   [  2.83500000e+01,   1.72670000e+05],
                   [  2.86500000e+01,   1.66430000e+05],
                   [  2.89500000e+01,   1.61550000e+05],
                   [  2.92500000e+01,   1.55230000e+05],
                   [  2.95500000e+01,   1.50620000e+05],
                   [  2.98500000e+01,   1.44690000e+05],
                   [  3.01500000e+01,   1.00740000e+05],
                   [  3.04500000e+01,   6.55290000e+04],
                   [  3.07500000e+01,   4.09580000e+04],
                   [  3.10500000e+01,   2.05860000e+04],
                   [  3.13500000e+01,   3.76200000e+03],
                   [  3.16500000e+01,   0.00000000e+00],
                   [  3.19500000e+01,   0.00000000e+00],
                   [  3.22500000e+01,   0.00000000e+00],
                   [  3.25500000e+01,   0.00000000e+00],
                   [  3.28500000e+01,   0.00000000e+00],
                   [  3.31500000e+01,   0.00000000e+00],
                   [  3.34500000e+01,   0.00000000e+00],
                   [  3.37500000e+01,   0.00000000e+00],
                   [  3.40500000e+01,   0.00000000e+00],
                   [  3.43500000e+01,   0.00000000e+00],
                   [  3.46500000e+01,   0.00000000e+00]])

    minval = np.sqrt(np.spacing(np.float32(1)))
    cj = interp.UnivariateSpline(d[:,0],
                                 d[:,1]+minval,
                                 k=1, s=0.5)

    # Interpolate the spline for our bin array into a background for our histogram
    bg = cj(bincenters)
    bg[bg < minval] = minval # do not allow division by zero

    return bg

def _bg_ksi3(bincenters):
    """ returns ksi1 background

    Parameters
    ----------
    bincenters : numpy array
        centers of bins for ksi1 background

    Returns
    -------
    bg : numpy array
        background number fractions for ksi1

    Notes
    -----
    Background is determined by cubic spline interpolation from a large
    numerical simulation containing 17156424 ksi values measured from
    random orientations. The bin spacing in this histogram is 0.3 degrees.

    """
    import numpy as np
    import scipy.interpolate as interp

    d = np.array([ [  0.00000000e+00,   0.00000000e+00],
                   [  1.50000000e-01,   2.02400000e+03],
                   [  4.50000000e-01,   6.40000000e+03],
                   [  7.50000000e-01,   1.04980000e+04],
                   [  1.05000000e+00,   1.44530000e+04],
                   [  1.35000000e+00,   1.87920000e+04],
                   [  1.65000000e+00,   2.27030000e+04],
                   [  1.95000000e+00,   2.67350000e+04],
                   [  2.25000000e+00,   3.07190000e+04],
                   [  2.55000000e+00,   3.50120000e+04],
                   [  2.85000000e+00,   3.91140000e+04],
                   [  3.15000000e+00,   4.28890000e+04],
                   [  3.45000000e+00,   4.71630000e+04],
                   [  3.75000000e+00,   5.13210000e+04],
                   [  4.05000000e+00,   5.43970000e+04],
                   [  4.35000000e+00,   5.88480000e+04],
                   [  4.65000000e+00,   6.29620000e+04],
                   [  4.95000000e+00,   6.75770000e+04],
                   [  5.25000000e+00,   7.14270000e+04],
                   [  5.55000000e+00,   7.46270000e+04],
                   [  5.85000000e+00,   7.91680000e+04],
                   [  6.15000000e+00,   8.24090000e+04],
                   [  6.45000000e+00,   8.65320000e+04],
                   [  6.75000000e+00,   9.06090000e+04],
                   [  7.05000000e+00,   9.47120000e+04],
                   [  7.35000000e+00,   9.79920000e+04],
                   [  7.65000000e+00,   1.01590000e+05],
                   [  7.95000000e+00,   1.05900000e+05],
                   [  8.25000000e+00,   1.09020000e+05],
                   [  8.55000000e+00,   1.12940000e+05],
                   [  8.85000000e+00,   1.17490000e+05],
                   [  9.15000000e+00,   1.20630000e+05],
                   [  9.45000000e+00,   1.23850000e+05],
                   [  9.75000000e+00,   1.27110000e+05],
                   [  1.00500000e+01,   1.33170000e+05],
                   [  1.03500000e+01,   1.35750000e+05],
                   [  1.06500000e+01,   1.38000000e+05],
                   [  1.09500000e+01,   1.42120000e+05],
                   [  1.12500000e+01,   1.45870000e+05],
                   [  1.15500000e+01,   1.49250000e+05],
                   [  1.18500000e+01,   1.52450000e+05],
                   [  1.21500000e+01,   1.56420000e+05],
                   [  1.24500000e+01,   1.59460000e+05],
                   [  1.27500000e+01,   1.63300000e+05],
                   [  1.30500000e+01,   1.65540000e+05],
                   [  1.33500000e+01,   1.68700000e+05],
                   [  1.36500000e+01,   1.71120000e+05],
                   [  1.39500000e+01,   1.75060000e+05],
                   [  1.42500000e+01,   1.78070000e+05],
                   [  1.45500000e+01,   1.80540000e+05],
                   [  1.48500000e+01,   1.83390000e+05],
                   [  1.51500000e+01,   1.85630000e+05],
                   [  1.54500000e+01,   1.89110000e+05],
                   [  1.57500000e+01,   1.91660000e+05],
                   [  1.60500000e+01,   1.95150000e+05],
                   [  1.63500000e+01,   1.97860000e+05],
                   [  1.66500000e+01,   1.99910000e+05],
                   [  1.69500000e+01,   2.03460000e+05],
                   [  1.72500000e+01,   2.05680000e+05],
                   [  1.75500000e+01,   2.08830000e+05],
                   [  1.78500000e+01,   2.10760000e+05],
                   [  1.81500000e+01,   2.12770000e+05],
                   [  1.84500000e+01,   2.14450000e+05],
                   [  1.87500000e+01,   2.17540000e+05],
                   [  1.90500000e+01,   2.19990000e+05],
                   [  1.93500000e+01,   2.22920000e+05],
                   [  1.96500000e+01,   2.24490000e+05],
                   [  1.99500000e+01,   2.25390000e+05],
                   [  2.02500000e+01,   2.27330000e+05],
                   [  2.05500000e+01,   2.29780000e+05],
                   [  2.08500000e+01,   2.31760000e+05],
                   [  2.11500000e+01,   2.34130000e+05],
                   [  2.14500000e+01,   2.34930000e+05],
                   [  2.17500000e+01,   2.36170000e+05],
                   [  2.20500000e+01,   2.39060000e+05],
                   [  2.23500000e+01,   2.39200000e+05],
                   [  2.26500000e+01,   2.40280000e+05],
                   [  2.29500000e+01,   2.41970000e+05],
                   [  2.32500000e+01,   2.42410000e+05],
                   [  2.35500000e+01,   2.44470000e+05],
                   [  2.38500000e+01,   2.44630000e+05],
                   [  2.41500000e+01,   2.45470000e+05],
                   [  2.44500000e+01,   2.46940000e+05],
                   [  2.47500000e+01,   2.46660000e+05],
                   [  2.50500000e+01,   2.48020000e+05],
                   [  2.53500000e+01,   2.49590000e+05],
                   [  2.56500000e+01,   2.50060000e+05],
                   [  2.59500000e+01,   2.49730000e+05],
                   [  2.62500000e+01,   2.49300000e+05],
                   [  2.65500000e+01,   2.49770000e+05],
                   [  2.68500000e+01,   2.51230000e+05],
                   [  2.71500000e+01,   2.50670000e+05],
                   [  2.74500000e+01,   2.50110000e+05],
                   [  2.77500000e+01,   2.53140000e+05],
                   [  2.80500000e+01,   2.50100000e+05],
                   [  2.83500000e+01,   2.49790000e+05],
                   [  2.86500000e+01,   2.49210000e+05],
                   [  2.89500000e+01,   2.48700000e+05],
                   [  2.92500000e+01,   2.47540000e+05],
                   [  2.95500000e+01,   2.47380000e+05],
                   [  2.98500000e+01,   2.47910000e+05],
                   [  3.01500000e+01,   1.88290000e+05],
                   [  3.04500000e+01,   1.38640000e+05],
                   [  3.07500000e+01,   1.07030000e+05],
                   [  3.10500000e+01,   8.20110000e+04],
                   [  3.13500000e+01,   5.92320000e+04],
                   [  3.16500000e+01,   4.13020000e+04],
                   [  3.19500000e+01,   2.75660000e+04],
                   [  3.22500000e+01,   1.75590000e+04],
                   [  3.25500000e+01,   9.84800000e+03],
                   [  3.28500000e+01,   4.43200000e+03],
                   [  3.31500000e+01,   1.59800000e+03],
                   [  3.34500000e+01,   1.18000000e+02],
                   [  3.37500000e+01,   0.00000000e+00],
                   [  3.40500000e+01,   0.00000000e+00],
                   [  3.43500000e+01,   0.00000000e+00],
                   [  3.46500000e+01,   0.00000000e+00]])

    minval = np.sqrt(np.spacing(np.float32(1)))
    cj = interp.UnivariateSpline(d[:,0],
                                 d[:,1]+minval,
                                 k=1, s=0.5)

    # Interpolate the spline for our bin array into a background for our histogram
    bg = cj(bincenters)
    bg[bg < minval] = minval # do not allow division by zero

    return bg

def ksi1hist(ksi1values, binwidth=0.5, return_everything=False):
    """ histogram ksi1 data

    Parameters
    ----------
    ksi1values : numpy array of floats
        ksi1values, in degrees

    binwidth : numpy float (optional)
        width of bins, in degrees

    return_everything : bool (optional)
        flag for whether

    Returns
    -------
    fcorr : numpy array
        background-corrected ksi1 number fractions

    bm : numpy array
        mean values for the histogram bins

    ba : (optional, with return_everything flag) numpy array
        histogram bin edges

    fraw : (optional, with return_everything flag) numpy array
        ksi1 number fractions before background correction
    """
    import numpy as np

    # Histogram ksi data
    ba = np.r_[0:35:binwidth]
    nb = ba.size
    bm = 0.5 * (ba[0:nb-1] + ba[1:nb])
    fraw, tmp = np.histogram(ksi1values, ba, density=True)
    fraw = fraw / np.sum(fraw)

    bg = _bg_ksi1(bm)

    # Correct the histogram for measurement
    fcorr = fraw / bg
    fcorr = fcorr / np.sum(fcorr)

    if return_everything == False:
        return fcorr, bm
    else:
        return fcorr, bm, ba, fraw

def ksi2hist(ksi2values, binwidth=0.5, return_everything=False):
    """ histogram ksi2 data

    Parameters
    ----------
    ksi1values : numpy array of floats
        ksi1values, in degrees

    binwidth : numpy float (optional)
        width of bins, in degrees

    return_everything : bool (optional)
        flag for whether

    Returns
    -------
    fcorr : numpy array
        background-corrected ksi1 number fractions

    bm : numpy array
        mean values for the histogram bins

    ba : (optional, with return_everything flag) numpy array
        histogram bin edges

    fraw : (optional, with return_everything flag) numpy array
        ksi1 number fractions before background correction
    """
    import numpy as np

    # Histogram ksi data
    ba = np.r_[0:35:binwidth]
    nb = ba.size
    bm = 0.5 * (ba[0:nb-1] + ba[1:nb])
    fraw, tmp = np.histogram(ksi2values, ba, density=True)
    fraw = fraw / np.sum(fraw)

    bg = _bg_ksi2(bm)

    # Correct the histogram for measurement
    fcorr = fraw / bg
    fcorr = fcorr / np.sum(fcorr)

    if return_everything == False:
        return fcorr, bm
    else:
        return fcorr, bm, ba, fraw

def ksi3hist(ksi3values, binwidth=0.5, return_everything=False):
    """ histogram ksi3 data

    Parameters
    ----------
    ksi3values : numpy array of floats
        ksi1values, in degrees

    binwidth : numpy float (optional)
        width of bins, in degrees

    return_everything : bool (optional)
        flag for whether

    Returns
    -------
    fcorr : numpy array
        background-corrected ksi1 number fractions

    bm : numpy array
        mean values for the histogram bins

    ba : (optional, with return_everything flag) numpy array
        histogram bin edges

    fraw : (optional, with return_everything flag) numpy array
        ksi1 number fractions before background correction
    """
    import numpy as np

    # Histogram ksi data
    ba = np.r_[0:35:binwidth]
    nb = ba.size
    bm = 0.5 * (ba[0:nb-1] + ba[1:nb])
    fraw, tmp = np.histogram(ksi3values, ba, density=True)
    fraw = fraw / np.sum(fraw)

    bg = _bg_ksi3(bm)

    # Correct the histogram for measurement
    fcorr = fraw / bg
    fcorr = fcorr / np.sum(fcorr)

    if return_everything == False:
        return fcorr, bm
    else:
        return fcorr, bm, ba, fraw

def fit_foldnorm(f, bm, return_info=False, accuracy=None, iterations=5000,
                 method = 'slsqp', mub=None, sgb=None):
    """ fit ksi histogram to folded normal distribution

    Parameters
    ----------
    fcorr : numpy float array
        histogram bin counts or frequencies

    bm : numpy float array
        middle values of bins

    return_info : bool {false}
        return info from fitting procedure

    accuracy : float or None {None}
        numerical accuracy of fitting procedure, default is floating point

    iterations : int
        maximum number of iterations in fitting minimization

    method : {'slsqp', 'lsq'}
        method for constrained function minimization, either sequential
        least squares programming (default) or standard least squares.

    mub : float iterable of size 2
        lower bound search value for folded normal mu at location 0;
        upper bound at location 1

    sgb : float iterable of size 2
        lower bound search value for folded normal mu at location 0;
        upper bound at location 1

    Returns
    -------
    mu   : float
        folded normal fit location parameter
    sig  : float
        folded normal fit shape parameter
    rss  : float
        final residual sum of squares difference between the fit and
        the histogram
    info : tuple
        (1) the number of iterations required for fitting
        (2) the exit mode from the optimizer (see notes)
        (3) message describing exit mode

    Notes
    -----
    - For very narrow distributions, it can be necessary to adjust the upper
      bound for mu. When the distribution is narrow and well-behaved, np.inf
      works fine (i.e., use mub=[np.spacing(np.float(1.0)), np.inf]).
    - Minimization is performed by default with scipy.optimize.fmin_slsqp.
      Refer to the documentation for this function for further details on exit
      mode returns.
    """
    import numpy as np
    from scipy.special import erf
    import scipy.optimize as spop

    if accuracy == None: # default to floating point accuracy
        accuracy = np.spacing(np.float32(0))

    # Estimate the distribution fitting parameters and then minimize the RMS
    # Note that RMS error is linearly related to log likelihood, so we should get
    # comparable results in both minimizations
    m1 = np.dot(bm, f)          # first moment
    m2 = np.dot(bm**2.0, f)     # second moment
    m4 = np.dot(bm**4.0, f)     # fourth moment
    s = np.dot((bm - m1)**2.0, f) # standard deviation
    x = np.arange(0.0,
                  np.amax(bm) + np.amax(bm) / 1000.0,
                  np.amax(bm) / 1000.0)

    # Elandt finds that the 1st and 2nd moments give a better approximation
    # when m1/s is greater than about 1.35, and the 2nd and 4th moment
    # method is better otherwise.
    if m1/s > 1.35: # estimate parameters by 1st and 2nd moments
        a  = m1**2.0 / m2 # ratio of moments
        def i0(x): # Eq 2 [2]
            return 1.0 - 0.5 * (1.0 + erf(x / np.sqrt(2.0)))
        def g(th, i0, a):
            return (np.sqrt(2.0/np.pi) * \
                    np.exp(-0.5 * th**2.0) - \
                    th * \
                    (1.0 - 2.0 * i0(-th)))**2.0 / (a * (1.0 + th**2.0)) - 1.0

        # Find the value of G closest to zero, limiting iterations
        # (fminbnd is not sufficiently robust in our case because G
        # sometimes is close to zero but does not necessarily intersect the
        # x axis in the limit of the half normal distribution). Assuming
        # the function will only come close to zero once, then we will just
        # search for the the solution to the equation iteratively
        lb  = 0.0
        ub  = np.amax(x)
        res = np.inf
        nit = 0
        j = 0
        while (res > accuracy and
                nit < iterations and
                (ub - lb) / 1000.0 > accuracy):

            # set a search range
            qtm = np.arange(lb, ub+(ub-lb)/1000.0, (ub-lb)/1000.0)

            tmp  = np.absolute(g(qtm, i0, a))
            resn = np.amin(tmp) # get value closest to zero in range
            j    = np.argmin(tmp)

            ub = qtm[np.amin([j + 1, qtm.size - 1])] # get new upper bound
            lb = qtm[np.amax([j - 1, 0])] # get new lower bound
            nit += 1 # increment number of iterations counter
        th=qtm[j]
        del qtm, j, ub, lb, res, g, a, i0, nit, tmp

    else: # Estimate Parameters by second and fourth moments

        b = m4 / m2**2.0 # ratio of moments

        def h(q, b): # Eq. 22 [2]
            return b * (1.0 + q**2.0)**2.0 - (3.0 + 6.0 * q**2.0 + q**4.0)

        # Find the value of h closest to zero, limiting iterations.
        # Same method as above.
        lb  = 0.0
        ub  = np.amax(x)
        res = np.inf
        nit = 0
        j   = 0
        while (res > accuracy and
                nit < iterations and
                np.absolute(ub-lb) / 1000.0 > accuracy):

            # set a search range
            qtm = np.arange(lb, ub+(ub-lb)/1000.0, (ub-lb)/1000.0)

            tmp = np.absolute(h(qtm, b))
            resn = np.amin(tmp) # get value closest to zero in range
            j    = np.argmin(tmp)
            if np.isreal(resn):
                res = resn # in case min cannot be found
            ub  = qtm[np.amin([j + 1, qtm.size - 1])] # get new upper bound
            lb  = qtm[np.amax([j - 1, 0])] # get new lower bound
            nit += 1 # increment number of iterations counter

        th = qtm[j]
        del qtm, j, ub, lb, res, h, b, nit, tmp
    del m2, m4

    # Compute initial guess values
    sgi = np.sqrt((s**2.0 + m1**2.0) / (1.0 + th**2.0))
    mui = th * sgi

    # Define minimization function as a least squares problem
    def ff(params):
        import numpy as np
        import scipy.stats as stats
        fx = stats.foldnorm.pdf(bm, params[0], loc=0, scale=params[1])
        fx = fx/(np.sum(fx)+np.spacing(1.0))
        return np.sum((f-fx)**2.0)

    # Set upper and lower bounds
    if mub==None:
        # np.inf would be a nicer upper bound here but realistically
        # we need to hold it under 15 for distributions that
        # are not well-behaved.
        mub = (np.sqrt(accuracy), 15.0)

    if sgb==None:
        # np.inf would be a nicer upper bound here but realistically
        # we need to hold it under 5 for distributions that
        # are not well-behaved
        # Zero would be a nicer lower bound but realistically it
        # we shouldn't allow it to be less than the bin width.
        sgb = (np.mean(np.diff(bm)), 5.0)

    if method=='slsqp':
    # ----- Begin perform constrained minimization, option #1 -----------------
        out, rss, its, imode, smode = spop.fmin_slsqp(ff, x0=[mui, sgi],
                                                      bounds=[mub, sgb],
                                                      full_output=True,
                                                      iter=iterations,
                                                      acc=accuracy,
                                                      iprint=1)
        mu  = out[0]
        sig = out[1]
    # ----- End perform constrained minimization, option #1 -----------------
    else:
    # ----- Begin perform constrained minimization, option #2 -----------------
        out = spop.least_squares(ff, x0=[mui, sgi],
                                 bounds=([mub[0], sgb[0]],
                                         [mub[1], sgb[1]]),
                                 verbose=1,
                                 max_nfev=iterations,
                                 gtol=accuracy,
                                 xtol=accuracy,
                                 ftol=accuracy)
        mu = out.x[0]
        sig = out.x[1]
        rss = out.cost
        imode = out.status
        smode = out.message
        its = out.nfev
    # ----- End perform constrained minimization, option #2 -----------------

    if return_info==False:
        return mu, sig
    else:
        return mu, sig, rss, (its, imode, smode)

def foldnorm_mode(mu, sig, accuracy=None):
    """ mode/peak of the folded normal distribution

    Parameters
    ----------
    mu : float
        location parameter
    sig : float
        shape parameter

    Returns
    -------
    expected_value : float
        E(Y)

    Notes
    -----
    Equation comes from the documentation for the folded normal distribution in
    the R statistics package.
    http://hosho.ees.hokudai.ac.jp/~kubo/Rdoc/library/VGAM/html/fnormal1.html
    """
    import scipy.stats as stats
    import numpy as np
    import scipy.optimize as spop

    if accuracy == None: # default to floating point accuracy
        accuracy = np.spacing(np.float32(1.0))

    fn = stats.foldnorm(mu, loc=0, scale=sig)
    def chk(x):
        val = -fn.pdf(x)
        return val

    res = spop.minimize(chk, 0.0, method='nelder-mead',
                        options={'xtol': 1e-8, 'disp': False})

    return res.x[0]


def fit_ksivals(ksi1vals, ksi2vals, ksi3vals, binwidth=0.5, makeplot=True, \
                accuracy=None, plotdpi=96, fmax=15.0, bmax=15.0,
                method = 'slsqp', mub=None, sgb=None, plotlegcontains='all'):
    """
    Parameters
    ----------

    ksi1vals : numpy float array

    ksi2vals : numpy float array

    ksi3vals : numpy float array

    makeplot : bool (optional, default=True)
        return a plot of the ksi value distributions

    binwidth : float (optional, default=0.5)
        width of bins on plot (if makeplot=True)

    fmax : float (optional, default=15.0)
        maximum value of y axis of plots (if makeplot=True)

    bmax : float (optional, default=15.0)
        maximum value of y axis of plots (if makeplot=True)

    plotdpi : int (optional, default=96)
        dpi for plot output. Use for screen resolution.

    accuracy : float (optional)
        defined accuracy of least squares optimization on folded normal fit.
        Default value is floating point accuracy.

    method : {'slsqp', 'lsq'}
        method for constrained function minimization, either sequential
        least squares programming (default) or standard least squares.

    mub : float iterable of size 2
        lower bound search value for folded normal mu at location 0;
        upper bound at location 1

    sgb : float iterable of size 2
        lower bound search value for folded normal mu at location 0;
        upper bound at location 1

    plotlegcontains : list of strings
                      {'all', 'mo', 'sig', 'mu', 'ksibar', 'std', 'rss'}
        defines which parameters to include in the plot legend

    Returns
    -------
    result : numpy float array
        np.array([[mo1, mean1, mu1, sig1, rss1],
                  [mo2, mean2, mu2, sig2, rss2],
                  [mo3, mean3, mu3, sig3, rss3]]); i.e.,
        - result[0, 0] : modal value of folded normal fit of ksi1
        - result[0, 1] : mean value of ksi1 histogram
        - result[0, 2] : mu value for folded normal fit of ksi1
        - result[0, 3] : sigma value for folded normal fit of ksi1
        - result[0, 4] : rss for folded normal fit of ksi1
        - result[1, 0] : modal value of folded normal fit of ksi2
        - result[1, 1] : mean value of ksi2 histogram
        - result[1, 2] : mu value for folded normal fit of ksi2
        - result[1, 3] : sigma value for folded normal fit of ksi2
        - result[1, 4] : rss for folded normal fit of ksi2
        - result[2, 0] : modal value of folded normal fit of ksi3
        - result[2, 1] : mean value of ksi3 histogram
        - result[2, 2] : mu value for folded normal fit of ksi3
        - result[2, 3] : sigma value for folded normal fit of ksi3
        - result[2, 4] : rss for folded normal fit of ksi3

    fig : (optional, returned if makeplot=True)
        figure matplotlib handle

    ax1 : (optional, returned if makeplot=True)
        ksi1 matplotlib axis handle

    ax2 : (optional, returned if makeplot=True)
        ksi2 matplotlib axis handle

    ax3 : (optional, returned if makeplot=True)
        ksi2 matplotlib axis handle

    Notes
    -----
    - If plot is output, the 'observation probability corrected' histogram is
      shown in white in the foreground. The darker gray background histogram is
      the uncorrected ksi values, and the lighter gray is the observation
      probability distribution.

    """

    import numpy as np
    import scipy.stats as stats
    import matplotlib.pyplot as plt

    if plotlegcontains=='all':
        plotlegcontains = ['mo', 'sig', 'mu', 'ksibar', 'std', 'rss']

    f1c, bm, ba, f1r = ksi1hist(ksi1vals,
                                 binwidth, return_everything=True)

    f2c, bm, ba, f2r = ksi2hist(ksi2vals,
                                 binwidth, return_everything=True)

    f3c, bm, ba, f3r = ksi3hist(ksi3vals,
                                 binwidth, return_everything=True)

    mu1, sig1, rss1, dat1 = fit_foldnorm(f1c, bm, return_info=True,
                                         accuracy=accuracy,
                                         method=method,
                                         mub=mub,
                                         sgb=sgb)

    mu2, sig2, rss2, dat2 = fit_foldnorm(f2c, bm, return_info=True,
                                         accuracy=accuracy,
                                         method=method,
                                         mub=mub,
                                         sgb=sgb)

    mu3, sig3, rss3, dat3 = fit_foldnorm(f3c, bm, return_info=True,
                                         accuracy=accuracy,
                                         method=method,
                                         mub=mub,
                                         sgb=sgb)

    mo1 = foldnorm_mode(mu1, sig1)
    mo2 = foldnorm_mode(mu2, sig2)
    mo3 = foldnorm_mode(mu3, sig3)

    mean1 = np.dot(bm, f1c)
    mean2 = np.dot(bm, f2c)
    mean3 = np.dot(bm, f3c)

    returnarray = np.array([[mo1, mean1, mu1, sig1, rss1],
                            [mo2, mean2, mu2, sig2, rss2],
                            [mo3, mean3, mu3, sig3, rss3]])

    if makeplot==False:

        return returnarray

    else:

        # parameters we might want to adjust
        xlimv = [0.0, bmax]
        ylimv = [0.0, fmax]
        xtickstep = 3
        ytickstep = 3
        bgdistclr = [0.85, 0.85, 0.85]
        uncorrclr = [0.5, 0.5, 0.5]
        ksi1faceclr = [1.0, 1.0, 1.0]
        ksi2faceclr = [1.0, 1.0, 1.0]
        ksi3faceclr = [1.0, 1.0, 1.0]
        fontsz=12

        xv = np.arange(0.0, np.amax(bm) + np.amax(bm) / 1000.0,
                      np.amax(bm) / 1000.0)

        # Create an awesome plot
        fig = plt.figure(figsize=(15, 4), dpi=plotdpi)
        fig.patch.set_facecolor('w')

        # ksi1 subplot
        ax1 = fig.add_subplot(131)
        f = stats.foldnorm.pdf(xv, mu1, scale=sig1)
        f = xv.size*f/np.sum(f)
        ax1.plot(xv, f, '-k', linewidth=1)
        ax1.bar(bm, f1r*bm.size, width=binwidth,
                facecolor=uncorrclr, edgecolor=uncorrclr)
        ax1.bar(bm, f1c*bm.size, width=binwidth,
                facecolor=ksi1faceclr, edgecolor=[0.0, 0.0, 0.0])
        bg = _bg_ksi1(xv)
        bg = bg/np.sum(bg)
        ax1.fill_between(x=xv, y1=bg*bg.size, y2=0, facecolors=bgdistclr, edgecolors=bgdistclr)
        legstring = []
        if any('mo' in s for s in plotlegcontains):
            legstring.append(
                 r'$Mo_{\xi1}$' +
                 r' = {0:5.2f}'.format(mo1) +
                 r'$^\circ$' +'\n')
        if any('mu' in s for s in plotlegcontains):
            legstring.append(
                 r'$\mu_{\xi1}$' + r' = {0:5.2f}'.format(mu1) +
                 r'$^\circ$' +'\n')
        if any('sig' in s for s in plotlegcontains):
            legstring.append(
                 r'$\sigma_{\xi1}$' + r' = {0:5.2f}'.format(sig1) +
                 r'$^\circ$' + '\n')
        if any('ksibar' in s for s in plotlegcontains):
            legstring.append(
                 r'$<\xi_1>$' + r' = {0:5.2f}'.format(np.mean(ksi1vals)) +
                 r'$^\circ$' +'\n')
        if any('std' in s for s in plotlegcontains):
            legstring.append(
                 r'$St.Dev.(\xi_1)$' + r' = {0:5.2f}'.format(np.std(ksi1vals)) +
                 r'$^\circ$' + '\n')
        if any('rss' in s for s in plotlegcontains):
            legstring.append(
                 r'$RSS$' + r' = {0:5.2e}'.format(rss1))
        ax1.text(0.97, 0.97,
                 ''.join(legstring),
                 fontsize=fontsz,
                 bbox=dict(edgecolor='k', facecolor='w'),
                 transform = ax1.transAxes,
                 horizontalalignment='right',
                 verticalalignment='top')
        ax1.set_xlim(xlimv)
        ax1.set_ylim(ylimv)
        ax1.set_xticks(np.arange(xlimv[0], xlimv[1]+xtickstep, xtickstep))
        ax1.set_yticks(np.arange(ylimv[0], ylimv[1]+ytickstep, ytickstep))

        # ksi2 subplot
        ax2 = fig.add_subplot(132)
        f = stats.foldnorm.pdf(xv, mu2, scale=sig2)
        f = xv.size*f/np.sum(f)
        ax2.plot(xv, f, '-k', linewidth=1)
        ax2.bar(bm, f2r*bm.size, width=binwidth,
                facecolor=uncorrclr, edgecolor=uncorrclr)
        ax2.bar(bm, f2c*bm.size, width=binwidth,
                facecolor=ksi2faceclr, edgecolor=[0.0, 0.0, 0.0])
        bg = _bg_ksi2(xv)
        bg = bg/np.sum(bg)
        ax2.fill_between(x=xv, y1=bg*bg.size, y2=0, facecolors=bgdistclr, edgecolors=bgdistclr)
        legstring = []
        if any('mo' in s for s in plotlegcontains):
            legstring.append(
                 r'$Mo_{\xi2}$' +
                 r' = {0:5.2f}'.format(mo2) +
                 r'$^\circ$' +'\n')
        if any('mu' in s for s in plotlegcontains):
            legstring.append(
                 r'$\mu_{\xi2}$' + r' = {0:5.2f}'.format(mu2) +
                 r'$^\circ$' +'\n')
        if any('sig' in s for s in plotlegcontains):
            legstring.append(
                 r'$\sigma_{\xi_2}$' + r' = {0:5.2f}'.format(sig2) +
                 r'$^\circ$' + '\n')
        if any('ksibar' in s for s in plotlegcontains):
            legstring.append(
                 r'$<\xi_2>$' + r' = {0:5.2f}'.format(np.mean(ksi2vals)) +
                 r'$^\circ$' +'\n')
        if any('std' in s for s in plotlegcontains):
            legstring.append(
                 r'$St.Dev.(\xi2)$' + r' = {0:5.2f}'.format(np.std(ksi2vals)) +
                 r'$^\circ$' + '\n')
        if any('rss' in s for s in plotlegcontains):
            legstring.append(
                 r'$RSS$' + r' = {0:5.2e}'.format(rss2))
        ax2.text(0.03, 0.97,
                 ''.join(legstring),
                 fontsize=fontsz,
                 bbox=dict(edgecolor='k', facecolor='w'),
                 transform = ax2.transAxes,
                 horizontalalignment='left',
                 verticalalignment='top')
        ax2.set_xlim(xlimv)
        ax2.set_ylim(ylimv)
        ax2.set_xticks(np.arange(xlimv[0], xlimv[1]+xtickstep, xtickstep))
        ax2.set_yticks([])


        # ksi2 subplot
        ax3 = fig.add_subplot(133)
        f = stats.foldnorm.pdf(xv, mu3, scale=sig3)
        f = xv.size*f/np.sum(f)
        ax3.plot(xv, f, '-k', linewidth=1)
        ax3.bar(bm, f3r*bm.size, width=binwidth,
                facecolor=uncorrclr, edgecolor=uncorrclr)
        ax3.bar(bm, f3c*bm.size, width=binwidth,
                facecolor=ksi3faceclr, edgecolor=[0.0, 0.0, 0.0])
        bg = _bg_ksi3(xv)
        bg = bg/np.sum(bg)
        ax3.fill_between(x=xv, y1=bg*bg.size, y2=0, facecolors=bgdistclr, edgecolors=bgdistclr)
        legstring = []
        if any('mo' in s for s in plotlegcontains):
            legstring.append(
                 r'$Mo_{\xi3}$' +
                 r' = {0:5.2f}'.format(mo3) +
                 r'$^\circ$' +'\n')
        if any('mu' in s for s in plotlegcontains):
            legstring.append(
                 r'$\mu_{\xi3}$' + r' = {0:5.2f}'.format(mu3) +
                 r'$^\circ$' +'\n')
        if any('sig' in s for s in plotlegcontains):
            legstring.append(
                 r'$\sigma_{\xi3}$' + r' = {0:5.2f}'.format(sig3) +
                 r'$^\circ$' + '\n')
        if any('ksibar' in s for s in plotlegcontains):
            legstring.append(
                 r'$<\xi_3>$' + r' = {0:5.2f}'.format(np.mean(ksi3vals)) +
                 r'$^\circ$' +'\n')
        if any('std' in s for s in plotlegcontains):
            legstring.append(
                 r'$St.Dev.(\xi_3)$' + r' = {0:5.2f}'.format(np.std(ksi3vals)) +
                 r'$^\circ$' + '\n')
        if any('rss' in s for s in plotlegcontains):
            legstring.append(
                 r'$RSS$' + r' = {0:5.2e}'.format(rss3))
        ax3.text(0.03, 0.97,
                 ''.join(legstring),
                 fontsize=fontsz,
                 bbox=dict(edgecolor='k', facecolor='w'),
                 transform = ax3.transAxes,
                 horizontalalignment='left',
                 verticalalignment='top')
        ax3.set_xlim(xlimv)
        ax3.set_ylim(ylimv)
        ax3.set_xticks(np.arange(xlimv[0], xlimv[1]+xtickstep, xtickstep))
        ax3.set_yticks([])

        ax1.set_ylabel(r'Relative Number Fraction', fontsize=fontsz*1.2)
        ax1.set_xlabel(r'$\xi_1$ ($^\circ$)', fontsize=fontsz*1.2)
        ax2.set_xlabel(r'$\xi_2$ ($^\circ$)', fontsize=fontsz*1.2)
        ax3.set_xlabel(r'$\xi_3$ ($^\circ$)', fontsize=fontsz*1.2)

        for label in (ax1.get_yticklabels() + ax1.get_xticklabels() +
                      ax2.get_xticklabels() + ax3.get_xticklabels()):
            label.set_fontsize(fontsz)

        fig.subplots_adjust(bottom=0.15)
        plt.tight_layout()
        plt.draw()

        return returnarray, \
                fig, ax1, ax2, ax3

# ----------------------------------------------------------------------------

def CSLcrit(name='Brandon', sigmas=None):
    ''' calculate CSL criteria

    **CSLcrit** calculates the angular deviation from the ideal misorientation
    relationship allowed by different *named* coincident site lattice (CSL)
    criteria of the form:

    ..math ::
                   \theta_c=\theta_0 / sigma^eta

    Parameters
    ----------
    name : {'Brandon', 'Pumphrey', 'Ishida-McLean', 'Palumbo-Aust'}

        CSL criteria name. Defaults to Brandon. Also accepts 'B' for Brandon,
        'P' for Pumphrey, 'IM' for Ishida-McLean, and 'PA' for Palumbo-Aust.
        String argument is case insensitive. Dash can be replaced by '&' or
        'and' with or without spaces around it.

    sigmas : None

        Sigma values of interest. If None, returns all values from 1 to 51.

    Returns
    -------
    thc : numpy float array

    Notes
    -----
    This is a Python translation of a Matlab function of the same name written
    by EJP 2011-05-02 while at Ruhr-Universitaet Bochum. This function was used
    for the analysis in references [1]_ and [2]_.

    For general reference on Coincident Site Lattice theory, see
    reference [3]_. The named CSL criteria used here from the following:
    Brandon [4]_, Pumphrey [5]_, Palumbo & Aust [6]_, and Ishida & McLean [7]_.

    References
    ----------
    .. [1] Otto et al, J Mater Sci 47 (2012) p2915-2927.
    .. [2] Otto et al, Acta Mater 60 (2012) p2982-2998.
    .. [3] Warrington & Boon, Acta Metall 23 (1975) p599
    .. [4] Brandon, Acta Metall 14 (1966) p1479
    .. [5] Pumphrey, in Grain Boundary Structure & Properties, London: Academic Press (1976) p139
    .. [6] Palumbo & Aust, Acta Metall 38 (1990) p2343
    .. [7] Ishida & McLean, Philos Mag 27 (1973) p1125
    '''
    import numpy as np
    import warnings

    # Remove
    name = name.replace('&', '').replace('and', '').replace('-','')
    name = name.replace(' ', '').lower()

    if sigmas == None:

        sig = np.arange(1, 53, 2, dtype=np.float)

    else:

        sig = np.asanyarray(sigmas, dtype=np.float)

        if np.any(sig % 2):

            warnings.warn('Only odd values of sigma are valid. ' + \
                          'Odd values have been removed.')

    if (name == 'brandon') or (name == 'b'):

        th0 = 15.0
        eta = 1.0 / 2.0

    elif (name == 'pumphrey') or (name =='p'):

        th0 = 15.0
        eta = 2.0 / 3.0

    elif (name == 'palumboaust') or (name == 'pa'):

        th0 = 15.0
        eta = 5.0 / 6.0

    elif (name == 'ishidamclean') or (name == 'im'):

        th0 = 8.0
        eta = 1.0

    else:

        warnings.warn('CSLcrit: unrecognized named criterion. ' + \
                      'Proceeding with Brandon criterion.')

        th0 = 15.0
        eta = 1.0 / 2.0

    # Convert theta to radians
    th0 = th0 * np.pi / 180.0

    # Calculate criteria
    thc = th0 / sig**eta

    return thc

# ----------------------------------------------------------------------------

#def genCSLs(sigmas=None):
#function [CSL,SQ,SV,SN]=GenCSL(varargin)
#% GenCSL generates the exact misorientations of the coincident site lattice
#% theory as rotation matrices and as quaternions.
#%
#% INPUTS:  varargin can be a list of CSL values (as in, for example,
#%          [3, 9, 27]). If left blank, it generates the CSLs up to Sigma49.
#%
#% OUTPUTS: the CSL in a structural format (CSL), or as lists of
#%          quaternions (SQ), their corresponding sigma values (SV), and
#%          their standard names (SN).
#%
#% REFERENCES:
#% [1] H. Grimmer, W. Bollman, D. H. Warrington. "Coincident-Site Lattices
#% and Complete Pattern-Shift Lattices in Cubic Crystals." Acta Cryst. A30
#% (1974) p197.
#% [2] D. H. Warrington, P. Bufalini. "The Coincident Site Lattice and Grain
#% Boundaries." Scripta Metall. 5 (1971) p771.
#%
#% NOTES: Function exactly reproduces Table 1 from Grimmer et al, except for
#% Sigma31b, which appears to be in error in the paper.
#%
#% -------------------------------------------------------------------------
#% 2011-05-01 | Eric Payton, Ruhr-Universitaet Bochum (payton.28[at]osu.edu)
#% -------------------------------------------------------------------------
#% This program is provided without any guarantee of correctness.
#% If you modify it and/or improve it, please kindly share with me your new
#% and improved version to the email address above. Thanks!
#%----------------------------------------------------------------------------



#    % Remove symmetric redundancies by putting into fundamental region
#    O=RMat2Quat(SM);RR=zeros(length(SM),4);
#    for ind=1:length(SM)
#        q=QuatProd(repmat(O(ind,:),size(SYM,1),1),SYM);
#        [an,ax]=Quat2AngAx(q);
#        ax=sort(abs(ax(an==min(an),:)),2,'descend');an=min(an);ax=ax(1,:);
#        RR(ind,:)=[an,ax];
#    end
#    [~,b]=unique(sigdec(RR,5),'rows');RR=RR(b,:);
#
#    % Sort the results by increasing angle in axis/angle description, such
#    % that, in the event of multiple results, th.a<th.b<th.c
#    [~,b]=sort(RR(:,1),'ascend');RR=RR(b,:);
#    RR=AngAx2Quat(RR(:,1),RR(:,2:4));
#
#    % Only allow minimum sigma solutions
#    % (otherwise sigma3 is also a sigma9 -- but not vice versa)
#    [~,b]=ismember(sigdec(RR,5),sigdec(TT,5),'rows');
#    if ~isempty(b),RR=RR(~b,:);end;TT=vertcat(TT,RR);
#
#    % Store results in structural array of TVecs; store names for output
#    if size(RR,1)>1
#        for ind=1:size(RR,1)
#            tmp=Quat2RMat(RR(ind,:));
#            eval(['CSL.S' num2str(S) char(96+ind) ...
#                '=tmp{1};']);
#            SN{ee}=['\Sigma' num2str(S) char(96+ind)];ee=ee+1;
#        end
#    else
#        tmp=Quat2RMat(RR);
#        eval(['CSL.S' num2str(S) '=tmp{1};']);
#        SN{ee}=['\Sigma' num2str(S)];ee=ee+1;
#    end
#
#    % Store results as a list of quaternions corresponding Sigma values
#    SQ=vertcat(SQ,RR);
#    SV=vertcat(SV,repmat(S,size(RR,1),1));
#
#    clear SM O an ax RR b ind  % clean up
#end % loop over all desired N
#clear Ndex S SYM TT X h1 k1 l1 h2 k2 l2 hkl1 hkl2 hkl3 q N ee % clean up
#
#% This is it.

# ----------------------------------------------------------------------------

def genKitahara_KS():
    ''' KS OR matrices

    **genKitahara_KS** returns the rotation matrices for the 24 variants of the
    Kurdjumov-Sachs orientation relationship in the order presented in Table 2
    of Kitahara [1]_.

    Parameters
    ----------
    None

    Returns
    -------
    V : 24 x 1 list of 3 x 3 numpy matrices

    Notes
    -----
    Rounding the results of Kitahara_KS and removing redundancies gives the
    list of Bain correspondence matrices.

    See also
    --------
    genKitahara_NW, genBainCorrMatrices

    References
    ----------
    .. [1] Kitahara et al, Acta Mater 54 (2006) p1279.
    '''
    import numpy as np

    Pg = [0] * 24
    Pa = [0] * 24
    Dg = [0] * 24
    Da = [0] * 24

    i =  0
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  0.0,  1.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i =  1
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  0.0,  1.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i =  2
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 0.0,  1.0, -1.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i =  3
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 0.0,  1.0, -1.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i =  4
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0, -1.0,  0.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i =  5
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0, -1.0,  0.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i =  6
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0,  0.0, -1.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i =  7
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0,  0.0, -1.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i =  8
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0, -1.0,  0.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i =  9
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0, -1.0,  0.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i = 10
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 0.0,  1.0,  1.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i = 11
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 0.0,  1.0,  1.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i = 12
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 0.0, -1.0,  1.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i = 13
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 0.0, -1.0,  1.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i = 14
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  0.0, -1.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i = 15
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  0.0, -1.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i = 16
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0,  1.0,  0.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i = 17
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0,  1.0,  0.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i = 18
    Pg[i] = [ 1.0,  1.0, -1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  1.0,  0.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i = 19
    Pg[i] = [ 1.0,  1.0, -1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  1.0,  0.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i = 20
    Pg[i] = [ 1.0,  1.0, -1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 0.0, -1.0, -1.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i = 21
    Pg[i] = [ 1.0,  1.0, -1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 0.0, -1.0, -1.0]
    Da[i] = [-1.0,  1.0, -1.0]

    i = 22
    Pg[i] = [ 1.0,  1.0, -1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0,  0.0,  1.0]
    Da[i] = [-1.0, -1.0,  1.0]

    i = 23
    Pg[i] = [ 1.0,  1.0, -1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0,  0.0,  1.0]
    Da[i] = [-1.0,  1.0, -1.0]

    V = []

    for i in range(0, 24):

        a = Dg[i] / np.linalg.norm(Dg[i])
        b = np.cross(Dg[i], Pg[i]) / np.linalg.norm(np.cross(Dg[i], Pg[i]))
        c = Pg[i] / np.linalg.norm(Pg[i])
        M = np.vstack([a, b, c])

        a = Da[i] / np.linalg.norm(Da[i])
        b = np.cross(Da[i], Pa[i]) / np.linalg.norm(np.cross(Da[i], Pa[i]))
        c = Pa[i] / np.linalg.norm(Pa[i])
        A = np.vstack([a, b, c])

        V.append(np.linalg.solve(A, M))

    return V

# ---------------------------------------------------------------------------

def genKitahara_NW():
    ''' NW OR matrices

    **genKitahara_NW** returns the rotation matrices for the 12 variants of the
    Nishiyama-Wasserman orientation relationship in the order presented in
    Table 2 of Kitahara [1]_.

    Parameters
    ----------
    None

    Returns
    -------
    V : 21 x 1 list of 3 x 3 numpy matrices

    Notes
    -----
    Table in Reference [1]_ has an error in the gamma direction for
    variant #10! For table 10 to agree with Table 6, the gamma direction must
    be [ 2 -1  1] and not [-1  2  1].

    See also
    --------
    genKitahara_KS, genBainCorrMatrices

    References
    ----------
    .. [1] Kitahara et al, Mater Charact (2006) p378.
    '''
    import numpy as np

    Pg = [0] * 12
    Pa = [0] * 12
    Dg = [0] * 12
    Da = [0] * 12

    i =  0
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 2.0, -1.0, -1.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i =  1
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  2.0, -1.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i =  2
    Pg[i] = [ 1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0, -1.0,  2.0]
    Da[i] = [ 0.0, -1.0,  1.0]


    i =  3
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-2.0, -1.0, -1.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i =  4
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0,  2.0, -1.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i =  5
    Pg[i] = [-1.0,  1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 1.0, -1.0,  2.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i =  6
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 2.0,  1.0, -1.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i =  7
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0, -2.0, -1.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i =  8
    Pg[i] = [ 1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  1.0,  2.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i =  9
    Pg[i] = [-1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [ 2.0, -1.0,  1.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i = 10
    Pg[i] = [-1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0,  2.0,  1.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    i = 11
    Pg[i] = [-1.0, -1.0,  1.0]
    Pa[i] = [ 0.0,  1.0,  1.0]
    Dg[i] = [-1.0, -1.0, -2.0]
    Da[i] = [ 0.0, -1.0,  1.0]

    V = []

    for i in range(0, 12):

        a = Dg[i] / np.linalg.norm(Dg[i])
        b = np.cross(Dg[i], Pg[i]) / np.linalg.norm(np.cross(Dg[i], Pg[i]))
        c = Pg[i] / np.linalg.norm(Pg[i])
        M = np.vstack([a, b, c])

        a = Da[i] / np.linalg.norm(Da[i])
        b = np.cross(Da[i], Pa[i]) / np.linalg.norm(np.cross(Da[i], Pa[i]))
        c = Pa[i] / np.linalg.norm(Pa[i])
        A = np.vstack([a, b, c])

        V.append(np.linalg.solve(A, M))

    return V

# ---------------------------------------------------------------------------

def genZirconia():
    
    from cryspy import rot
    import numpy as np
    #%% Lattice correspondence A

    Abca1 = rot.rmat(g11=0.0, g12=0.0, g13=1.0,
                    g21=1.0, g22=0.0, g23=0.0,
                    g31=0.0, g32=1.0, g33=0.0)
    
    Abca2 = rot.rmat(g11=0.0, g12=0.0, g13=1.0,
                    g21=-1.0, g22=0.0, g23=0.0,
                    g31=0.0, g32=-1.0, g33=0.0)
    
    Acba1 = rot.rmat(g11=0.0, g12=0.0, g13=1.0,
                    g21=0.0, g22=1.0, g23=0.0,
                    g31=-1.0, g32=0.0, g33=0.0)
    
    Acba2 = rot.rmat(g11=0.0, g12=0.0, g13=1.0,
                    g21=0.0, g22=-1.0, g23=0.0,
                    g31=1.0, g32=0.0, g33=0.0)
    
    Abca3 = rot.rmat(g11=0.0, g12=0.0, g13=-1.0,
                    g21=-1.0, g22=0.0, g23=0.0,
                    g31=0.0, g32=1.0, g33=0.0)
    
    Abca4 = rot.rmat(g11=0.0, g12=0.0, g13=-1.0,
                    g21=1.0, g22=0.0, g23=0.0,
                    g31=0.0, g32=-1.0, g33=0.0)
    
    Acba3 = rot.rmat(g11=0.0, g12=0.0, g13=-1.0,
                    g21=0.0, g22=-1.0, g23=0.0,
                    g31=-1.0, g32=0.0, g33=0.0)
    
    Acba4 = rot.rmat(g11=0.0, g12=0.0, g13=-1.0,
                    g21=0.0, g22=1.0, g23=0.0,
                    g31=1.0, g32=0.0, g33=0.0)
    
    corrA = [Abca1, Abca2, Acba1, Acba2, Abca3, Abca4, Acba3, Acba4]
    corrtA = ['Abca1', 'Abca2', 'Acba1', 'Acba2', 'Abca3', 'Abca4', 'Acba3', 'Acba4']


    Bcab1 = rot.rmat(g11=0.0, g12=1.0, g13=0.0,
                    g21=0.0, g22=0.0, g23=1.0,
                    g31=1.0, g32=0.0, g33=0.0)
    
    Bcab2 = rot.rmat(g11=0.0, g12=-1.0, g13=0.0,
                    g21=0.0, g22=0.0, g23=1.0,
                    g31=-1.0, g32=0.0, g33=0.0)
    
    Bacb1 = rot.rmat(g11=-1.0, g12=0.0, g13=0.0,
                    g21=0.0, g22=0.0, g23=1.0,
                    g31=0.0, g32=1.0, g33=0.0)
    
    Bacb2 = rot.rmat(g11=1.0, g12=0.0, g13=0.0,
                    g21=0.0, g22=0.0, g23=1.0,
                    g31=0.0, g32=-1.0, g33=0.0)
    
    Bcab3 = rot.rmat(g11=0.0, g12=1.0, g13=0.0,
                    g21=0.0, g22=0.0, g23=-1.0,
                    g31=-1.0, g32=0.0, g33=0.0)
    
    Bcab4 = rot.rmat(g11=0.0, g12=-1.0, g13=0.0,
                    g21=0.0, g22=0.0, g23=1.0,
                    g31=1.0, g32=0.0, g33=0.0)
    
    Bacb3 = rot.rmat(g11=-1.0, g12=0.0, g13=0.0,
                    g21=0.0, g22=0.0, g23=-1.0,
                    g31=0.0, g32=-1.0, g33=0.0)
    
    Bacb4 = rot.rmat(g11=1.0, g12=0.0, g13=0.0,
                    g21=0.0, g22=0.0, g23=-1.0,
                    g31=0.0, g32=1.0, g33=0.0)
    
    corrB = [Bcab1, Bcab2, Bacb1, Bacb2, Bcab3, Bcab4, Bacb3, Bacb4]
    corrtB = ['Bcab1', 'Bcab2', 'Bacb1', 'Bacb2', 'Bcab3', 'Bcab4', 'Bacb3', 'Bacb4']

    Cabc1 = rot.rmat(g11=1.0, g12=0.0, g13=0.0,
                    g21=0.0, g22=1.0, g23=0.0,
                    g31=0.0, g32=0.0, g33=1.0)
    
    Cabc2 = rot.rmat(g11=-1.0, g12=0.0, g13=0.0,
                    g21=0.0, g22=-1.0, g23=0.0,
                    g31=0.0, g32=0.0, g33=1.0)
    
    Cbac1 = rot.rmat(g11=0.0, g12=1.0, g13=0.0,
                    g21=-1.0, g22=0.0, g23=0.0,
                    g31=0.0, g32=0.0, g33=1.0)
    
    Cbac2 = rot.rmat(g11=0.0, g12=-1.0, g13=0.0,
                    g21=1.0, g22=0.0, g23=0.0,
                    g31=0.0, g32=0.0, g33=1.0)
    
    Cabc3 = rot.rmat(g11=-1.0, g12=0.0, g13=0.0,
                    g21=0.0, g22=1.0, g23=0.0,
                    g31=0.0, g32=0.0, g33=-1.0)
    
    Cabc4 = rot.rmat(g11=1.0, g12=0.0, g13=0.0,
                    g21=0.0, g22=-1.0, g23=0.0,
                    g31=0.0, g32=0.0, g33=-1.0)
    
    Cbac3 = rot.rmat(g11=0.0, g12=-1.0, g13=0.0,
                    g21=-1.0, g22=0.0, g23=0.0,
                    g31=0.0, g32=0.0, g33=-1.0)
    
    Cbac4 = rot.rmat(g11=0.0, g12=1.0, g13=0.0,
                    g21=1.0, g22=0.0, g23=0.0,
                    g31=0.0, g32=0.0, g33=-1.0)
    
    corrC = [Cabc1, Cabc2, Cbac1, Cbac2, Cabc3, Cabc4, Cbac3, Cbac4]
    corrtC = ['Cabc1', 'Cabc2', 'Cbac1', 'Cbac2', 'Cabc3', 'Cabc4', 'Cbac3', 'Cbac4']

    return np.hstack([corrA, corrB, corrC]), np.hstack([corrtA, corrtB, corrtC])

# ---------------------------------------------------------------------------  

def lattParamAus_Onink(at_pct_carbon, temperature_degC):
    ''' estimate austenite lattice parameter at temperature

    **latticeParamAus_Onink** estimates the lattice parameter of the
    face-centered cubic austenite phase at high temperature as a function of
    carbon content based on the best fit equations of Onink et al [1]_.
    Lattice parameter is returned in angstroms.

    Parameters
    ----------
    at_pct_carbon : float
        atomic percent of carbon

    temperature_degC : float
        temperature in degrees celsius

    Returns
    -------
    a0 : float
        lattice parameter of austenite (gamma) phase in angstroms

    Notes
    -----
    Verified validity range varies with temperature and carbon content
    (Ref [1]_, Table 2):

    =======  ======================
    at % C   Temperature Range in K
    =======  ======================
     0.05          1080-1250
     1.30          1080-1250
     1.75          1060-1250
     2.60          1030-1250
     3.79          1000-1250
    =======  ======================

    References
    ----------
    .. [1] M Onink, C. M. Brakman, F. D. Tichelaar, E. J. Mittemeijer,
           S. van der Zwaag. "The lattice parameters of austenite and ferrite
           in Fe-C alloys as functions of carbon concentration and
           temperature." Scripta Metall Mater 29 (1993) p1011-1016.

    Examples
    --------
    >>> # Reproduce Figure 2 from Onink et al [1]_. The x-axis will be similar
    >>> # but not *exactly* the same because we are plotting atomic percent
    >>> # instead of number C atoms / 100 Fe atoms.
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy import constants
    >>> temperature_degC = constants.K2C(np.arange(1000.0, 1275.0, 25.0))
    >>> for at_pct_carbon in [0.0, 1.3, 1.75, 2.6, 3.65]:
    >>>     print(at_pct_carbon)
    >>>     a0 = latticeParamAus_Onink(at_pct_carbon, temperature_degC)
    >>>     # For consistency in cryspy, our inputs are in degC
    >>>     # and output is in angstroms.
    >>>     # Converting for comparison to the published plot...
    >>>     plt.plot(constants.C2K(temperature_degC), a0/10.0)
    >>> plt.legend(['Fe', 'Fe-0.3C', 'Fe-0.4C', 'Fe-0.6C', 'Fe-0.8C'])
    >>> plt.ylim([0.3630, 0.3690])
    >>> plt.xlim([950, 1350])
    >>> plt.xlabel('Temperature (K)')
    >>> plt.ylabel('Lattice parameter of austenite (nm)')
    >>> plt.show()
    '''
    from scipy import constants
    tempK = constants.C2K(temperature_degC)
    a0 = (0.363067 + 0.000783 * at_pct_carbon) * \
        (1.0 + (24.92 - 0.51 * at_pct_carbon) * 1.0E-6 * (tempK - 1000.0))
    return a0 * 10.0

# ---------------------------------------------------------------------------

def lattParamMart_Kurdjumov(c_wt_pct, a0_angstroms):
    ''' RTP Fe alpha prime parameters as function of wt%C

    **lattParams_Kurdjumov** calculates room temperature lattice parameters
    for austenite and martensite from carbon content.

    Parameters
    ----------
    c_wt_pct : float
        weight percent carbon
    a0 : float
        austenite lattice parameter in angstroms

    Returns
    -------
    a : float
        martensite a axis lattice parameter in angstroms
    c : float
        martensite c axis lattice parameter in angstroms

    Notes
    -----
    Equations come from Ref. [1]_ with reference to Ref. [2]_ and [3]_ for
    discovery, [4]_ and [5]_ for validation, and [6]_ and [7]_ for showing that
    the same equations are acceptable for alloy steels.

    References
    ----------
    .. [1] Z. Fan et al. Phys Rev B 52(14) (1995) p9979.
    .. [2] G. V. Kurdjumow and E. Kaminsky, Nature 122 (1928) p425.
    .. [3] G. V. Kurdjumov and E. Kaminsky, Z Phys 53 (1929) p696.
    .. [4] G. V. Kurdjumov, J ISIL 195 (1960) p26.
    .. [5] C. S. Roberts, Trans AIME 197 (1953) p203
    .. [6] Z. Nishiyama and M. Doi, J. JIM 8 (1944) p305
    .. [7] P. G. Winchell and M. Cohen, Trans AIME 55 (1962) p347.
    '''
    # c in atomic percent

    a = a0_angstroms - 0.013 * c_wt_pct
    c = a0_angstroms + 0.116 * c_wt_pct

    return a, c

# ----------------------------------------------------------------------------

def lattParams_RobertsCohen(c_at_pct):
    ''' RTP Fe alpha prime & gamma parameters as function of at%C

    **lattParams_RobertsCohen** calculates room temperature lattice parameters
    for austenite and martensite from carbon content.

    Parameters
    ----------
    c_at_pct : float
        atomic percent carbon

    Returns
    -------
    a0 : float
        austenite lattice parameter in angstroms
    a : float
        martensite a axis lattice parameter in angstroms
    c : float
        martensite c axis lattice parameter in angstroms

    Notes
    -----
    Atomic pct equations come from Ref. [1]_ using the data obtained in
    Ref. [2]_.

    References
    ----------
    .. [1] M. Cohen, "The Strengthening of Steel." Transactions of the
           Metallurgical Society of AIME 224 (1962) p637-657.
    .. [2] C.S. Roberts, "Technical Note: Effect of Carbon on the Volume
           Fractions and Lattice Parameters Of Retained Austenite and
           Martensite." Journal of Metals, AIME Transactions (Institute of
           Metals Division) 197 (1953). p203-204.
    '''
    # c in atomic percent

    a0 = 3.548 + 0.0094 * c_at_pct # units of angstroms
    a = 2.861 - 0.0028 * c_at_pct
    c = 2.861 + 0.0247 * c_at_pct

    return a0, a, c

# ----------------------------------------------------------------------------

def lattParams_Lee(c_at_pct, temp_degC):
    ''' Temperature dependent alpha & gamma prime parameters as function of %C

    Parameters
    ----------
    c_at_pct : float
        atomic percent carbon in austenite
    temp_degC : float
        temperature in degrees C

    Returns
    -------
    a0 : float
        austenite lattice parameter in angstroms
    a : float
        martensite a axis lattice parameter in angstroms
    c : float
        martensite c axis lattice parameter in angstroms

    Notes
    -----
    We assume in this function that the carbon content in austenite and
    martensite has to be the same, since the martensite forms from the
    austenite.

    References
    ----------
    .. [1] S-J Lee et al. Acta Mater 55 (2007) p875
    '''

    tt = temp_degC + 273.15

    betaG = (24.9 - 0.5 * c_at_pct) * 1.0E-6
    a0 = (0.36306 + 7.83E-4 * c_at_pct) * (1.0 + betaG * (tt - 1000.0))

    betaM = (14.9 - 1.9 * c_at_pct) * 1.0E-6
    c = (0.28610 + 0.0025855 * c_at_pct) * (1.0 + betaM * (tt - 273.0))
    a = (0.28610 - 0.0002898 * c_at_pct) * (1.0 + betaM * (tt - 273.0))

    return a0, a, c

# ---------------------------------------------------------------------------

def martstart_YangBhadeshiaAdjustment(G, Ms0=363.5, f=0.01, m=0.05, a=1, b=0.2689):
    '''
    Predicts reduction in martensite start temperature with reduction in
    austenite grain size according to geometric model.

    Parameters
    ----------
    Ms0 : float
        Martensite start temperature at limit of infinitely large grains, which
        is a function of composition.
    f : float
        Detectable martensite fraction
    m : float
        Martensite plate aspect ratio
    a : float
        fitting constant, units 1/mm^3
    b : float
        fitting constant, unitless
    G : float
        ASTM grain size number

    Returns
    -------
    Ms : float
        Martensite start temperature, in degC

    Notes
    -----
    All input values default to those from [1]_. The user should note that Ms0
    will vary with composition.


    References
    ----------
    ..[1] Yang & Bhadeshia, Scripta Mater 60 (2009) p493
    ..[2] p90-91 of Underwood, Quantitative Stereology. Reading, Mass: Addison-
          Wesley Publishing Co., 1970.
    '''
    import numpy as np

    # Convert the ASTM grain size to the mean grain section area in sq microns
    A_gamma = 100.0 * 25.4**2 / 2.0**(G - 1.0)
    
    # Convert to sq millimeters
    A_gamma = A_gamma / 1000.0**2

    # For lack of a better option at this time, we will assume spherical grains
    # and use the equivalent spherical volume for the mean cross sectional
    # area. The mean sectional area of a sphere is 2/3*pi*r^2. If we set this
    # equal to unity, then r=sqrt(3/2/pi). The volume of a sphere with this
    # radius is 4*pi/3*(3/2/pi)^3/2
    V_gamma = A_gamma * 4.0 * np.pi / 3.0 * (3.0 / 2.0 / np.pi)**(3.0 / 2.0)

    Ms = Ms0 - (1.0 / b) * np.log((1.0 / (a * V_gamma)) * \
                (np.exp(-np.log(1.0 - f) / m) - 1.0) + 1.0)

    return Ms

# ---------------------------------------------------------------------------

def martstart_vanBohemen2017(G, C=0.0, Mn=0.0, Cr=0.0, Ni=0.0,
                               Mo=0.0, Si=0.0, K_C=370.0, K_HP=350.0,
                               K_1=1015.0, D_C=15.0):
    ''' predicts the martensite start temperature from composition and austenite grain size

    Parameters
    ----------
    G : float
        Austenite ASTM grain size number
    C : float
        wt% C
    Mn : float
        wt% Mn
    Cr : float
        wt% Cr
    Ni : float
        wt% Ni
    Mo : float
        wt% Mo
    Si : float
        wt% Si
    K_C : float
        J/mol
    K_HP : float
        J*(microns^0.5)/mol
    D_C : float
        Critical austenite equivalent diameter grain size for single packet
        microstructures, with units of microns
    K_1 : float {1015.0}
        A constant depending on defect size, strain, and interfacial energyies,
        with units of J/mol


    Returns
    -------
    Ms : float
        martensite start temperature in degC

    Notes
    -----
    The model does not account for the compositional effects of Nb, Ti, V, Al,
    Cu,W and Co, since limited data is available for a thorough validation to
    evaluate the coefficients Km for these elements [1]_.


    References
    ----------
    .. [1] van Bohemen & Morsdorf, Acta Mater 125 (2017) p401.
    '''
    import numpy as np

    # The following are the fitting parameters and constants
    # used by van Bohemen & Morsdorf

    # Convert ASTM GS to "diameter" in microns.
    # Note that this is not really a diameter, but the square root of area.
    D_gamma = np.sqrt(100.0 * 25.4**2 / 2.0**(G - 1.0))

    # Non-chemical free energy term from the reduction in martensite lath
    # aspect ratio that occurs as grain size is reduced
    W_C = K_C * np.exp(-6.0 * D_gamma / D_C)

    # Non-chemical free energy term associated with Hall-Petch hardening
    # of the austenite
    W_HP = K_HP / np.sqrt(D_gamma)

    # Compositional dependence of athermal frictional work
    # (as fit to experimental data)
    W_mu = 670.0 * np.sqrt(C) + \
          np.sqrt((195.0 * np.sqrt(Mn))**2 + \
                  (140.0 * np.sqrt(Si))**2 + \
                  (170.0 * np.sqrt(Cr))**2 + \
                  (  5.0 * np.sqrt(Ni))**2 + \
                  (205.0 * np.sqrt(Mo))**2)

    # Estimated change in Gibbs energy due to defect content, athermal
    # frictional work, Hall-Petch hardening, and martensite lath aspect ratio
    delGc = K_1 + W_mu + W_HP + W_C

    # van Bohemen's best fit to FactSage prediction of martensite start
    # temperatures
    T1 = 718.3 - 291.0*C - 24.0*Mn - 1.8*Si - 5.6*Cr - 18.4*Ni + 3.5*Mo

    #
    ms = T1 - delGc / 7.22

    return ms

# ---------------------------------------------------------------------------

def martstart_vanBohemen2012(C=0.0, Mn=0.0, Cr=0.0, Ni=0.0,
                               Mo=0.0, Si=0.0):
    ''' predicts the martensite start temperature from composition

    Parameters
    ----------
    C : float
        wt% C
    Mn : float
        wt% Mn
    Cr : float
        wt% Cr
    Ni : float
        wt% Ni
    Mo : float
        wt% Mo
    Si : float
        wt% Si

    Returns
    -------
    Ms : float
        martensite start temperature in degC

    Notes
    -----
    Non-linear fit may be better for high C than Steven-Haynes.


    References
    ----------
    .. [1] van Bohemen, Mater Sci Tech 28 (2012) p487.
    '''
    import numpy as np
    sumKixi = 31.*Mn + 13.*Si + 10.*Cr + 18.*Ni + 12.*Mo
    ms = 565.0 - sumKixi - 600.0*(1.0-np.exp(-0.96*C))
    return ms

# ----------------------------------------------------------------------------

def martstart_StevenHaynesMod(C=0.0, Mn=0.0, Cr=0.0, Ni=0.0,
                               Mo=0.0, Co=0.0, Si=0.0):
    ''' predicts the martensite start temperature from composition

    Parameters
    ----------
    C : float
        wt% C
    Mn : float
        wt% Mn
    Cr : float
        wt% Cr
    Ni : float
        wt% Ni
    Mo : float
        wt% Mo
    Co : float
        wt% Co
    Si : float
        wt% Si

    Returns
    -------
    Ms : float
        martensite start temperature in degC

    Notes
    -----
    - Implements modification from Ref. _[2] to equation originally from
    Ref. _[1].
    - The equation should work reasonably well under the limits of
    0.6 wt%C, 4.9 wt% Mn, 12.2 wt% Cr, 12.5 wt% Ni, and 5.4 wt% Mo, 3.2 wt% Si,
    and 9 wt% Co.
    - Underestimates start temperature at high Mo contents

    References
    ----------
    .. [1] W. Steven & A.G. Haynes. JISI 183 (1956) p349.
    .. [2] C.Y. Kung & J.J. Rayment, MMTA 13A (1982) p328.
    .. [3] G. Krauss, Steels: Processing, Structure, and Performance.
           Materials Park, OH: ASM International (2008). (Ch 5, p64)
    '''
    Ms = 561.0 - 474.0 * C - 33.0 * Mn - 17.0 * Cr - 17.0 * Ni - 21.0 * Mo +\
         (10.0 * Co - 7.5 * Si)
    return Ms


def martstart_AndrewsLinearMod(C=0.0, Mn=0.0, Cr=0.0, Ni=0.0,
                               Mo=0.0, Co=0.0, Si=0.0):
    ''' predicts the martensite start temperature from composition

    Parameters
    ----------
    C : float
        wt% C
    Mn : float
        wt% Mn
    Cr : float
        wt% Cr
    Ni : float
        wt% Ni
    Mo : float
        wt% Mo
    Co : float
        wt% Co
    Si : float
        wt% Si

    Returns
    -------
    Ms : float
        martensite start temperature in degC

    Notes
    -----
    - Implements modification from Ref. _[2] to equation originally from
    Ref. _[1].
    - The equation should work reasonably well under the limits of
    0.6 wt%C, 4.9 wt% Mn, 12.2 wt% Cr, 12.5 wt% Ni, and 5.4 wt% Mo, 3.2 wt% Si,
    and 9 wt% Co.
    - Overestimates martensite start temperature for high chromium alloys

    References
    ----------
    .. [1] K.W. Andrews, JISI 203 (1965) p721.
    .. [2] C.Y. Kung & J.J. Rayment, MMTA 13A (1982) p328.
    .. [3] G. Krauss, Steels: Processing, Structure, and Performance.
           Materials Park, OH: ASM International (2008). (Ch 5, p64)
    '''
    Ms = 539.0 - 423.0 * C - 30.4 * Mn - 12.1 * Cr - 17.7 * Ni - 7.5 * Mo +\
         (10.0 * Co - 7.5 * Si)
    return Ms

def martfract_HarrisCohen(deltaT):
    ''' returns transformed fraction martensite as a function of undercooling

    Parameters
    ----------
    deltaT : float
        Amount of undercooling in degC

    Notes
    -----
    - Developed for steels containing 1.1% C

    References
    ----------
    .. [1] G. Krauss, Steels: Processing, Structure, and Performance.
           Materials Park, OH: ASM International (2008). (Ch 5, p65)
    .. [2] W.H. Harris & M. Cohen. Trans AIME 180 (1949) p447.

    See Also
    --------
    cryspy.orrl.martstart_AndrewsLinearMod
    cryspy.orrl.martstart_StephenHaynesMod
    '''
    f = 1.0 - 6.96 * 1.0E-15 * (455.0 - deltaT)**5.32

    return f

# ----------------------------------------------------------------------------

def martfract_KoistenenMarburger(deltaT):
    ''' returns transformed fraction martensite as a function of undercooling

    Parameters
    ----------
    deltaT : float
        Amount of undercooling in degC

    Notes
    -----
    - Developed for steels containing 0.37-1.1 wt% C

    References
    ----------
    .. [1] G. Krauss, Steels: Processing Structure and Performance.
           Materials Park, OH: ASM International (2008). (Ch 5, p65)
    .. [2] D.P. Koistenen & R.E. Marburger. Acta Metall 7 (1959) p59.

    See Also
    --------
    cryspy.orrl.martstart_AndrewsLinearMod
    cryspy.orrl.martstart_StephenHaynesMod
    '''
    import numpy as np

    f = 1.0 - np.exp(-1.1E-2 * deltaT)

    return f

# ----------------------------------------------------------------------------

def genBainCorrMatrices():
    ''' Bain corr matrices

    **genBainCorrMatrices** returns the rotation matrices for the 12 Bain
    correspondence matrices, numbered according to parallelism of the c-axis
    of the B phase.

    Parameters
    ----------
    None

    Returns
    -------
    V : 12 x 1 list of 3 x 3 numpy matrices

    Notes
    -----
    Rounding the results of Kitahara_KS and removing redundancies gives the
    list of Bain correspondence matrices.

    See also
    --------
    genKitahara_NW, genBainCorrMatrices

    References
    ----------
    .. [1] Kitahara et al, Acta Mater 54 (2006) p1279.
    '''
    import numpy as np

    # Numbered according to parallelism of c-axis of B phase
    c = 12 * [0]

    c[ 0] =  np.mat([[ 1.0, -1.0,  0.0],
                     [ 1.0,  1.0,  0.0],
                     [ 0.0,  0.0,  1.0]])

    c[ 1] =  np.mat([[-1.0, -1.0,  0.0],
                     [ 1.0, -1.0,  0.0],
                     [ 0.0,  0.0,  1.0]])

    c[ 2] =  np.mat([[ 1.0,  1.0,  0.0],
                     [-1.0,  1.0,  0.0],
                     [ 0.0,  0.0,  1.0]])

    c[ 3] =  np.mat([[-1.0,  1.0,  0.0],
                     [ 1.0,  1.0,  0.0],
                     [ 0.0,  0.0, -1.0]])

    c[ 4] =  np.mat([[-1.0,  0.0,  1.0],
                     [ 1.0,  0.0,  1.0],
                     [ 0.0,  1.0,  0.0]])

    c[ 5] =  np.mat([[ 1.0,  0.0, -1.0],
                     [ 1.0,  0.0,  1.0],
                     [ 0.0, -1.0,  0.0]])

    c[ 6] =  np.mat([[-1.0,  0.0, -1.0],
                     [-1.0,  0.0,  1.0],
                     [ 0.0,  1.0,  0.0]])

    c[ 7] =  np.mat([[ 1.0,  0.0,  1.0],
                     [ 1.0,  0.0, -1.0],
                     [ 0.0,  1.0,  0.0]])

    c[ 8] =  np.mat([[ 0.0,  1.0, -1.0],
                     [ 0.0,  1.0,  1.0],
                     [ 1.0,  0.0,  0.0]])

    c[ 9] =  np.mat([[ 0.0,  1.0,  1.0],
                     [ 0.0, -1.0,  1.0],
                     [ 1.0,  0.0,  0.0]])

    c[10] =  np.mat([[ 0.0, -1.0,  1.0],
                     [ 0.0,  1.0,  1.0],
                     [-1.0,  0.0,  0.0]])

    c[11] =  np.mat([[ 0.0, -1.0, -1.0],
                     [ 0.0,  1.0, -1.0],
                     [ 1.0,  0.0,  0.0]])

    return c

# ----------------------------------------------------------------------------

def calcPTMT(a0, a, c, delta=1.0):
    ''' perform PTMT calculations

    **calcPTMT** performs calculations of the phenomenological theory of
    martensitic transformations using inputs of the lattice parameters of the
    face-centered cubic (gamma) and body-centered-tetragonal (alpha prime)
    phases.

    Parameters
    ----------
    a0 : float
        lattice parameter of the FCC austenite phase in angstroms.
    a : float
        lattice parameter for the a and b axes of the BC(C/T) martensite phase
        in angstroms.
    c : float
        lattice parameter for the c axis of the BC(C/T) martensite phase in
        angstroms.

    Returns
    -------
    ksivals : n x 3 numpy float array
        ksi angles for the orientation relationship

    detailtuple : tuple of floats
        detailed results of the ptmt calculation. In order, the returns are:
        (0) invariantPlaneNormal, (1) displacementVector,
        (2) shapeDeformationMatrix, (3) shapeDeformationMagnitude,
        (4) complimentaryShearDirection, (5) complimentaryShearMagnitude,
        (6) complimentaryShearAngleRadians

    References
    ----------
    .. [1] Z. Nishiyama, Martensitic Transformation. New York: Academic Press (1978).
    .. [2] J. S. Bowles and J. K. Mackenzie. "The crystallography of martensite
           transformations I." Acta Metall 2 (1954).
    .. [3] J. K. Mackenzie and J. S. Bowles. "The crystallography of martensite
           transformations II." Acta Metall 2 (1954).
    .. [4] H. K. D. H. Bhadeshia, Worked Examples in the Geometry of Crystals,
           2nd Ed. Brookfield, VT: The Institute of Metals, 1987.
    .. [5] C. M. Wayman, Introduction to the Crystallography of Martensitic
           Transformations. New York: Macmillan Company (1964).

    Notes
    -----
    Provided for ease of copy-and-paste, due to the unweildy outputs

    #from cryspy import orrl
    #ksivals, resultsTuple = orrl.calcPTMT(a0, a, c)
    #invariantPlaneNormal = resultsTuple[0]
    #displacementVector = resultsTuple[1]
    #shapeDeformationMatrix = resultsTuple[2]
    #shapeDeformationMagnitude = resultsTuple[3]
    #complimentaryShearDirection = resultsTuple[4]
    #complimentaryShearMagnitude = resultsTuple[5]
    #complimentaryShearAngleRadians = resultsTuple[6]

    #from cryspy import orrl
    #import numpy as np
    #a0, a, c = orrl.lattParams_RobertsCohen(np.arange(0, 5, 0.5))
    #ksivals = []
    #invariantPlaneNormal = []
    #displacementVector = []
    #shapeDeformationMatrix = []
    #shapeDeformationMagnitude = []
    #complimentaryShearDirection = []
    #complimentaryShearMagnitude = []
    #complimentaryShearAngleRadians = []
    #for i in np.arange(0, np.shape(a0)[0]):
    #    tmp, resultsTuple = orrl.calcPTMT(a0[i], a[i], c[i])
    #    ksivals.append(tmp)
    #    invariantPlaneNormal.append(resultsTuple[0])
    #    displacementVector.append(resultsTuple[1])
    #    shapeDeformationMatrix.append(resultsTuple[2])
    #    shapeDeformationMagnitude.append(resultsTuple[3])
    #    complimentaryShearDirection.append(resultsTuple[4])
    #    complimentaryShearMagnitude.append(resultsTuple[5])
    #    complimentaryShearAngleRadians.append(resultsTuple[6])

    Examples
    -------
    >>> # Reproduce example from Nishiyama (with lattice parameters from Wayman)
    >>> a0 = 3.591 # angstroms, gamma
    >>> a  = 2.875 # angstroms, alpha prime
    >>> c  = 2.875 # angstroms, alpha prime
    >>> ksivals, resultsTuple = calcPTMT(a0, a, c)
    >>> invariantPlaneNormal = resultsTuple[0]
    >>> displacementVector = resultsTuple[1]
    >>> shapeDeformationMatrix = resultsTuple[2]
    >>> shapeDeformationMagnitude = resultsTuple[3]
    >>> complimentaryShearDirection = resultsTuple[4]
    >>> complimentaryShearMagnitude = resultsTuple[5]
    >>> complimentaryShearAngleRadians = resultsTuple[6]

    >>> from cryspy import orrl
    >>> import numpy as np
    >>> a0, a, c = orrl.lattParams_RobertsCohen(np.arange(0, 5, 0.5))
    >>> ksivals = []
    >>> invariantPlaneNormal = []
    >>> displacementVector = []
    >>> shapeDeformationMatrix = []
    >>> shapeDeformationMagnitude = []
    >>> complimentaryShearDirection = []
    >>> complimentaryShearMagnitude = []
    >>> complimentaryShearAngleRadians = []
    >>> for i in np.arange(0, np.shape(a0)[0]):
    >>>     tmp, resultsTuple = orrl.calcPTMT(a0, a, c)
    >>>     ksivals.append(tmp)
    >>>     invariantPlaneNormal.append(resultsTuple[0])
    >>>     displacementVector.append(resultsTuple[1])
    >>>     shapeDeformationMatrix.append(resultsTuple[2])
    >>>     shapeDeformationMagnitude.append(resultsTuple[3])
    >>>     complimentaryShearDirection.append(resultsTuple[4])
    >>>     complimentaryShearMagnitude.append(resultsTuple[5])
    >>>     complimentaryShearAngleRadians.append(resultsTuple[6])
    '''

    import numpy as np
    from scipy.linalg import lstsq as sp_lstsq
    from warnings import warn

    # set up dilatation parameters from inputs
    # I added the anisotropic dilatation later, so it is implemented differently.
    # The anisotropic case (del1 and del2) follows Otte, Acta Cryst 16 (1963) p8.
    # The isotropic case follows the implementation description in Wayman's text.
    if np.size(delta) == 1:
        del1 = delta
        del2 = delta
    elif np.size(delta) == 2:
        del1 = delta[0]
        del2 = delta[1]

#    # In Python 3.6, I'm getting kernel death failures for bad parameter inputs.
#    # Putting forth these warnings.
#    if a0 <= 0:
#        warn('Input parameters outside of acceptable ranges.\n'+
#             '(a0 <= 0)')
#
#    if a0 / a <= c / a:
#        warn('Input parameters outside of acceptable ranges.\n'+
#             '(a0/a <= c/a)')
#
#    tmpval = 4.0 * (a**2 / a0**2)**2 + \
#             (4.0 * (c**2 / a0**2) - 8.0) * (a**2 / a0**2) + \
#             2.0 - (c**2 / a0**2)**2
#    if tmpval > 0:
#        warn('Input parameters outside of acceptable ranges.\n'+
#             '4x^2 + (4y-8)x + 2-y^2) > 0, where x=a^2/a0^2 and y=c^2/a0^2')
#
#    if a / a0 < -0.27 * c / a + 1.086:
#        warn('Input parameters outside of acceptable ranges.\n'+
#             '(a/a0 < -0.27 * ca + 1.086)')

    '''
    During PTMT calculations, several steps impose limits on the ranges of the
    ratios a/a0, c/a, and c/a0. For example, to generate the Bain distortion
    matrix, a0 must be greater than zero; in order to be able to calculate an
    invariant line, a0/a must be less than or equal to c/a and
    4x^2 + (4y-8)x + 2-y^2) must be less than or equal to zero, where x is a^2/a0^2
    and y is c^2/a0^2; for calculation of the invariant plane normal, both c/a0 and
    a/a0 must be greater than zero, and so on. The border between the valid and
    invalid ranges is found to follow a/a0 < -0.27 * c/a + 1.086 with an r-squared
    value in excess of 0.9999.
    '''

    # Using a0=3.591 and a=c=2.875 does not agree exactly with top of page 358 in
    # Nishiyama (values from Wayman). Better with a=c=2.8747428847860585?
    # Should we change both a and a0?
    # a0 = 3.591 # angstroms, gamma
    # a  = 2.875 # angstroms, alpha prime
    # c  = 2.875 # angstroms, alpha prime

    slipPlaneFamilyBCT = [1, 1, 2]
    slipDirectionFamilyBCT = [1, 1, -1]

    #%% Prepare inputs
    cc = genBainCorrMatrices()

    slipPlaneFamilyBCT = np.mat(slipPlaneFamilyBCT / np.linalg.norm(slipPlaneFamilyBCT))
    slipDirectionFamilyBCT = np.mat(slipDirectionFamilyBCT / np.linalg.norm(slipDirectionFamilyBCT))

    ## Below not used, but could probably be done...
    ## Loop through Bain matrices
    # O = [] # empty list to contain variant rotations
    # for baindex in range(0, np.shape(cc)[0]):
    cm = cc[0]#cc[baindex]

    #%% Create the Bain distortion matrix, bb
    # an anisotropic dilatation would be added here, a la Otte, Acta Cryst 1963
    # Wayman discusses this possibility on p122-123.
    eta1 = np.sqrt(2.0) * a / a0
    eta3 = c / a0

    # Override for checking against Bhadeshia's example calculation
    #eta1 = 1.136071
    #eta3 = 0.803324

    # In agreement with p123 Wayman; p358 Nishiyama

    bb = np.mat([[eta1, 0.0, 0.0],
                   [0.0, eta1, 0.0],
                   [0.0, 0.0, eta3]])
    # this is bold B in Nishiyama; (fBf) in Wayman

    #%% Calculation of invariant lines.
    # The equations below are hand-solved from the system of equations shown
    # as Eq. 6 on page 358 of Nishiyama.
    # Resulting values agree with Equations 6' in Nishiyama
    x1 = -np.sqrt((1.0 - eta1**2) / (eta3**2 - eta1**2))
    x2a = -np.sqrt(1.0 - 2.0 * x1**2)
    x2b = -x2a
    x3 = -x1

    xin = [np.mat([x1, x2a, x3]), np.mat([x1, x2b, x3])] # two possibilities
    # results are in agreement with x1 and x2 in Wayman; xi1 and xi2 in Nishiyama

    #%% Calculation of invariant normal
    # Resulting values agree with Equations 7' in Nishiyama
    eta1i = 1.0 / eta1**2
    eta3i = 1.0 / eta3**2
    n1 = np.sqrt((1.0 - eta1i) / (eta3i - eta1i))
    n2a = np.sqrt(1.0 - 2.0 * n1**2)
    n2b = -n2a
    n3 = n1

    nin = [np.mat([n1, n2a, n3]), np.mat([n1, n2b, n3])] # two possibilities
    # results in agreement with ni1' and ni2' in Nishiyama; n1' and n2' in Wayman
    # and the solutions on page 62 of Bhadeshia (with different values of the etas)

    #%% Lines after Bain distortion
    # There are four possible combinations of x and n (top of page 359).
    # Nishiyama only looks at one as an example in his textbook.
    # Resulting value for xui agrees with equation 8, p359
    # Resulting value for pp2 agrees with equation 9, p359

    # One could choose one of the four combinations this way
    # for xindex in np.arange(0,2):
    #    for nindex in np.arange(0,2):
    # I'm choosing not to do this here because we can get the ksi angles from
    # just the first index, then populate the set of variants from the ksi angles.
    xindex = 0
    nindex = 0

    xi = xin[xindex]
    ni = nin[nindex]

    # Convert BCT slip plane and family into FCC coordinates
    slipPlaneFCC = slipPlaneFamilyBCT * cm
    slipDirectionFCC = np.linalg.solve(cm, slipDirectionFamilyBCT.T)

    p2 = slipPlaneFCC / np.linalg.norm(slipPlaneFCC)
    d2 = slipDirectionFCC / np.linalg.norm(slipPlaneFCC)
    xui = xi * bb # x underscore sub i
    pp2 = p2 * bb.I
    pp2 = pp2 / np.linalg.norm(pp2) # p prime underscore sub 2

    #%% Invariant line strain
    # The numerical absolute values of u and v agree with bottom of page 359.
    # The middle value of u is negative in the text but positive here.
    # It appears that this is probably an error in the text.
    u = np.cross(xi, p2)
    v = np.cross(xui, pp2)

    # R1 and R2 seem to be transposed in Nishiyama vs. Wayman
    rr1 = np.mat(np.vstack([xi,   p2, u]).T)
    rr2 = np.mat(np.vstack([xui, pp2, v]).T)

    # iS_0i, Eq 13 in Nishiyama. Agrees with values in Wayman, p126.
    yy = rr2.T * bb * rr1

    # Nishiyama Eq. 15. His result is obtained with ni[0] and xi[0].
    # Agrees with result in Wayman, p126.
    zz = np.dot(ni, rr1).A1

    # Now, to produce the solution, we can use the associative properties of
    # the matrices to multiply (n'_i;i) by the matrix containing the rotation
    # of beta about x_i (see Nishiyama Eqs 15 and 14, respectively.)
    # When we multiply through and collect the necessary terms, we find:
    tmp1 = np.array([[zz[1] * yy[1,0] + zz[2] * yy[2,0],
                      zz[2] * yy[1,0] - zz[1] * yy[2,0]],
                     [zz[1] * yy[1,1] + zz[2] * yy[2,1],
                      zz[2] * yy[1,1] - zz[1] * yy[2,1]],
                     [zz[1] * yy[1,2] + zz[2] * yy[2,2],
                      zz[2] * yy[1,2] - zz[1] * yy[2,2]]])

    tmp2 = np.array([[zz[0] - zz[0] * yy[0, 0]],
                     [zz[1] - zz[0] * yy[0, 1]],  # Note that the problem is
                     [zz[2] - zz[0] * yy[0, 2]]]) # now overspecified!


    # Now we can solve the (overspecified) system for cosB and sinB
    # If our system was not overspecified, we could use
    #       np.linalg.solve(tmp1, tmp2)
    # (as we would in Matlab.)

    # It appears that matlab must convert to a least squares method
    # for an overspecified system.

    # The numpy least squares function at the time of this writing has
    # a bug that results in a segfault, so we are using the scipy version.

    # Numpy segfault issue documented here:
    #    https://github.com/numpy/numpy/issues/9891
    if (~np.any(np.isnan(tmp1)) or ~np.any(np.isnan(tmp2))):
        tmp3 = sp_lstsq(tmp1, tmp2)[0]
        cosB = tmp3[0]
        sinB = tmp3[1]
    else:
        cosB = 1.0
        sinB = 0.0
    # This result is in agreement with Wayman as well.

    betamx = np.mat([[1.0, 0.0, 0.0], [0.0, cosB, -sinB], [0.0, sinB, cosB]])

    # Create a dilatation parameter matrix
    # See page 122-123 of Wayman text. Here we allow a tetragonal anisotropic
    # dilation.
    deltamx = np.mat(np.eye(3))
    deltamx[0,0] = del1
    deltamx[1,1] = del1
    deltamx[2,2] = del2
    #deltamx = np.mat([[1.2*dil, 0., 0.], [0., 1.2*dil, 0.], [0., 0., 0.9*dil]])

    # The results for ss here agree with Eq 17.
    # Also sgrees with (fSf) in Wayman p126.
    # And agrees with (F S F) in Bhadeshia p63 (with different etas)
    ss = deltamx * rr1 * betamx * yy * rr1.T
    ## Check results. Should be <= round off error, approx 1E-15
    #print(np.mat(zz) * np.mat(betamx) * np.mat(yy) - np.mat(zz))

    # Note that ss is also given by (fJf) * bb, where (fJf) is given by
    # J =  rr1 * betamx * rr2.T

    #%% Calculate the useful stuff

    # Calculate the invariant plane normal
    p1 = p2 * np.linalg.inv(ss) - p2
    p1p = p1 / np.linalg.norm(p1) # agrees with Eq. 18 in Nishiyama and p127 Wayman

    # Calculate the displacement vector
    d1 = (ss * d2 - d2) / (p1p * d2) # agrees with Eq. 19 Nishiyama and p127 Wayman

    # Magnitude of the shape deformation
    m1 = np.linalg.norm(d1) # agrees with Eq. 20

    # Direction of complimentary shear
    y = np.mat(np.cross(p1p, np.mat([1.0, 0.0, 0.0]))).T # [1 0 0] in text p361
    d2 = (y - np.linalg.solve(ss, y)) / (p2 * y)

    # Magnitude of complimentary shear
    m2 = np.linalg.norm(d2)

    # Angle of complimentary shear
    alpha = np.arctan(0.5 * m2)

    # Shape deformation matrix (Bhadeshia)
    delta = m1 * p1p * d1
    s = np.sqrt(m1**2 - delta**2) # FIXME: Refer back to Bhadeshia, where is this variable to be used?
    pp = np.eye(3, 3) + m1 * d1 * p1p

    #%% Get orientation relationship between martensite and austenite
    # Note that this gives the CORRECT orientation relationship for the c/a and
    # a/a0 ratios for K-S.

    # Victoria said: "KSI is similarity-transformed version of J (i.e.,
    #       expressed in martensite basis)"
    # This would then be solved as follows:
    #      gamma = cm.A / np.tile(np.linalg.norm(cm.A, axis=1), [3, 1]).T
    #      vx = gamma * (a / a0) * aJg * gamma.T # some variant of the OR
    # I don't find that this returns the expected result.


    # Bhadeshia writes the following, but the determinant isnt unity
    aJg = cm * ss.I
    # However, if we multiply by a/a0, the expected det=1 matrix is returned
    # Matrices compare well to the T matrices in Kitahara

    vx = (a / a0) * aJg # some variant of the orientation relationship
    #print(np.linalg.det(vx)) # check that determinant is unity


    # This works for only certain variants. For example, it gets the right
    # ksi angles for Kitahara T4 and T24 in the KS paper; however, it does not
    # get the expected result for Kitahara T18 in the KS paper.
    # Returns expected results for Kitahara T1 and T11 and in NW paper.
    #vx = np.matrix([[0.742, -0.667, -0.075],
    #                [0.650, 0.742, -0.167], # check using Kitahara matrix
    #                [0.167, 0.075, 0.983]]) # for KS (#T24)
    #vx = np.matrix([[0.667, -0.742, 0.075],
    #                [0.742, 0.650, -0.167], # check using Kitahara matrix
    #                [0.075, 0.167, 0.983]]) # for KS (#T4)
    #vs = np.matrix([[-0.075, -0.742, 0.667],
    #                [0.167, 0.650, 0.742],  # check using Kitahara matrix
    #                [-0.983, 0.167, 0.075]])# for KS (#T18)
    #vx = np.matrix([[0.000, 0.707, -0.707],
    #                [-0.169, 0.697, 0.697], # check using a Kitahara matrix
    #                [0.986, 0.120, 0.120]]) # for NW
    #vx = np.matrix([[-0.707, 0.000, -0.707],
    #                [-0.120, -0.986, 0.120], # check using a Kitahara matrix
    #                [-0.697, 0.169, 0.697]]) # for NW

    # Now, to get the ksi angles, we first need the Bain correspondence for
    # the representative variant
    bt = np.round(vx)

    # The orientation of the variant matrix matters. NOT transposing gives the
    # expected result for the Kitahara KS and NW matrices, so it appears that
    # I am using the correct orientation (aus -> mart).
    # bt = bt.T
    # vx = vx.T
    # print(vx)
    # FIXME: Eventually this would be much more robust if we were to orient
    #        all variants against a standard variant and then get the ksi
    #        angles.

    g1 = vx.ravel()[0, 0:3]
    b1 = bt.ravel()[0, 0:3]
    g2 = vx.ravel()[0, 3:6]
    b2 = bt.ravel()[0, 3:6]
    g3 = vx.ravel()[0, 6:9]
    b3 = bt.ravel()[0, 6:9]

    # individual ksi angles
    tmpksi1 = np.real(np.arccos(np.dot(g1 / np.linalg.norm(g1),
                                       b1.T / np.linalg.norm(b1))))
    tmpksi2 = np.real(np.arccos(np.dot(g2 / np.linalg.norm(g2),
                                       b2.T / np.linalg.norm(b2))))
    ksi1 = np.amin([tmpksi1, tmpksi2]) * 180. / np.pi
    ksi2 = np.amax([tmpksi1, tmpksi2]) * 180. / np.pi
    ksi3 = np.real(np.arccos(np.dot(g3 / np.linalg.norm(g3),
                                    b3.T / np.linalg.norm(b3)))) * 180. / np.pi

    #print(ksi1, ksi2, ksi3)
    # TODO: Generate full set of variant matrices as an output

    #%% Compile results together
    ksivals = np.array([ksi1, ksi2, ksi3])
    invariantPlaneNormal = p1p
    displacementVector = d1
    shapeDeformationMatrix = pp
    shapeDeformationMagnitude = m1
    complimentaryShearDirection = d2
    complimentaryShearMagnitude = m2
    complimentaryShearAngleRadians = alpha

    return ksivals, \
            (invariantPlaneNormal, displacementVector, \
             shapeDeformationMatrix, shapeDeformationMagnitude, \
             complimentaryShearDirection, complimentaryShearMagnitude, \
             complimentaryShearAngleRadians)

# ----------------------------------------------------------------------------

def solvePTMTparams(ksi1, ksi2, ksi3, accuracy='default',
                    initial_guess=(0.8035, 1.047),
                    max_iter=1000):
    """
    **solvePTMTparams** solves the PTMT using a least squares approach to
    return the apparent a/a0 and c/a ratios for a given set of ksi parameters.

    Parameters
    ----------
    ksi1 : float
        ksi1 value

    ksi2 : float
        ksi2 value

    ksi3 : float
        ksi3 value

    accuracy : float (optional)
        desired accuracy of the result.

    initial_guess : 2-tuple of floats (optional)
        Initial guess values. initial_guess[0] should be the a/a0 ratio and
        initial_guess[1] should be the c/a ratio. Defaults to (0.8035, 1.047).

    max_iter : int (optional) {1000}
        Maximum number of search iterations.

    Returns
    -------
    res : object
        results of the optimization routine. See scipy.optimize.least_squares
        for more outputs.
        - res.x = optimized values of [a/a0, c/a]
        - res.fun = residual sum of squares of solution
        - res.nfev = number of function evalutions performed
        - res.success = True if convergence criteria satisfied
        - res.status = reason code for algorithm termination
        - res.message = string reason for algorithm termination
    """
    import numpy as np
    import scipy.optimize as spop

    if accuracy=='default':
        accuracy = np.spacing(np.float64(1.0))

    # Initial guesses for a/a0 and c/a ratios
    aa0i = initial_guess[0]
    cai  = initial_guess[1]

    # Set up objective function for minimization
    def rmsPTMT(x):

        a0 = 3.6 # this value doesn't actually matter when there is no dilation
        a = a0 * x[0]
        c = x[1] * a

        try:
            tmp, resultsTuple = calcPTMT(a0, a, c)
        except:
            tmp = 999999999999999999999999. * np.ones(3)

        return np.sum(np.sqrt((ksi1 - tmp[0])**2 +
                              (ksi2 - tmp[1])**2 +
                              (ksi3 - tmp[2])**2))

    res = spop.minimize(rmsPTMT, x0=[aa0i, cai],
                        method='Nelder-Mead',
                        options={'maxiter': max_iter,
                                 'xatol': accuracy,
                                 'fatol': accuracy})

    # Also tried the following minimizations, and neither worked as reliably...

    #res = spop.least_squares(rmsPTMT, x0=[aa0i, cai], bounds=[aa0b, cab],
    #                         verbose=1, max_nfev=max_iter, gtol=accuracy,
    #                         xtol=accuracy, ftol=accuracy)

    #out, rss, its, imode, smode = spop.fmin_slsqp(rmsPTMT, x0=[aa0i, cai],
    #                                              bounds=[aa0b, cab],
    #                                              full_output=True,
    #                                              iter=iterations,
    #                                              acc=accuracy,
    #                                              #epsilon=1.E-12,
    #                                              iprint=1)

    return res

# ----------------------------------------------------------------------------

def solvePTMTparams_dil(ksi1, ksi2, ksi3, a0, accuracy='default',
                    initial_guess=(0.8035, 1.047, 0.996, 1.052),
                    max_iter=1000):
    """
    **solvePTMTparams** solves the PTMT using a least squares approach to
    return the apparent a/a0 and c/a ratios for a given set of ksi parameters.

    Parameters
    ----------
    ksi1 : float
        ksi1 value

    ksi2 : float
        ksi2 value

    ksi3 : float
        ksi3 value

    accuracy : float (optional)
        desired accuracy of the result.

    initial_guess : 2-tuple of floats (optional)
        Initial guess values. initial_guess[0] should be the a/a0 ratio and
        initial_guess[1] should be the c/a ratio. Defaults to (0.8035, 1.047).

    max_iter : int (optional) {1000}
        Maximum number of search iterations.

    Returns
    -------
    res : object
        results of the optimization routine. See scipy.optimize.least_squares
        for more outputs.
        - res.x = optimized values of [a/a0, c/a]
        - res.fun = residual sum of squares of solution
        - res.nfev = number of function evalutions performed
        - res.success = True if convergence criteria satisfied
        - res.status = reason code for algorithm termination
        - res.message = string reason for algorithm termination
    """
    import numpy as np
    import scipy.optimize as spop

    if accuracy=='default':
        accuracy = np.spacing(np.float64(1.0))

    # Initial guesses for a/a0 and c/a ratios
    aa0i = initial_guess[0]
    cai  = initial_guess[1]
    del1i = initial_guess[2]
    del2i = initial_guess[3]

    # Set up objective function for minimization
    def rmsPTMT(x):

        a = a0 * x[0]
        c = x[1] * a
        del1 = x[2]
        del2 = x[3]

        try:
            tmp, resultsTuple = calcPTMT(a0, a, c, delta=[del1, del2])
        except:
            tmp = 999999999999999999999999. * np.ones(3)

        # I am applying a simple weighting here because
        # experience is suggesting that ksi1 is harder to fit.
        return np.sum(np.sqrt(3.0*(ksi1 - tmp[0])**2 +
                              2.0*(ksi2 - tmp[1])**2 +
                              1.0*(ksi3 - tmp[2])**2))

    res = spop.minimize(rmsPTMT, x0=[aa0i, cai, del1i, del2i],
                        method='Nelder-Mead',
                        options={'maxiter': max_iter,
                                 'xatol': accuracy,
                                 'fatol': accuracy})

    return res
