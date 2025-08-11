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

def latticeParamAus_Onink(at_pct_carbon, temperature_degC):
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
    tempK = constants.convert_temperature(temperature_degC, 'Celsius', 'Kelvin')
    a0 = (0.363067 + 0.000783 * at_pct_carbon) * \
        (1.0 + (24.92 - 0.51 * at_pct_carbon) * 1.0E-6 * (tempK - 1000.0))
    return a0 * 10.0

# ---------------------------------------------------------------------------

def lattParams_RobertsCohen(c_at_pct):
    ''' RTP Fe alpha & gamma prime parameters as function of at%C

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

def calcPTMT(a0, a, c):
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
    eta1 = np.sqrt(2.0) * a / a0
    eta3 = c / a0

    bb = np.mat([[eta1, 0.0, 0.0],
                   [0.0, eta1, 0.0],
                   [0.0, 0.0, eta3]])

    #%% Calculation of invariant lines.
    # Resulting values agree with Equations 6' in Nishiyama
    x1 = -np.sqrt((1.0 - eta1**2) / (eta1**2 + eta3**2 - 2*eta1**2))
    x2a = -np.sqrt(1.0 - 2.0 * x1**2)
    x2b = -x2a
    x3 = -x1

    xin = [np.mat([x1, x2a, x3]), np.mat([x1, x2b, x3])] # two possibilities

    #%% Calculation of invariant normal
    # Resulting values agree with Equations 7' in Nishiyama
    eta1i = 1.0 / eta1**2
    eta3i = 1.0 / eta3**2
    n1 = np.sqrt((1.0 - eta1i) / (eta1i + eta3i - 2.0 * eta1i))
    n2a = np.sqrt(1.0 - 2.0 * n1**2)
    n2b = -n2a
    n3 = n1

    nin = [np.mat([n1, n2a, n3]), np.mat([n1, n2b, n3])] # two possibilities

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
    rr1 = np.mat(np.vstack([xi,   p2, u]).T)
    rr2 = np.mat(np.vstack([xui, pp2, v]).T)

    # iS_0i, Eq 13
    yy = rr2.T * bb * rr1

    # Nishiyama Eq. 15. His result is obtained with ni[0] and xi[0]
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
    # np.linalg.solve(tmp1, tmp2) as we would in Matlab.
    # It appears that matlab must convert to a least squares method
    # for an overspecified system
    tmp3 = np.linalg.lstsq(tmp1, tmp2)[0]
    cosB = tmp3[0][0]
    sinB = tmp3[1][0]
    
    print(cosB, sinB)

    betamx = np.mat([[1.0, 0.0, 0.0], [0.0, cosB, -sinB], [0.0, sinB, cosB]])

    # The results for ss here agree with Eq 17
    ss = rr1 * betamx * yy * rr1.T

    ## Check results. Should be <= round off error, approx 1E-15
    #print(np.mat(zz) * np.mat(betMX) * np.mat(yy) - np.mat(zz))

    #%% Calculate the useful stuff

    # Calculate the invariant plane normal
    p1 = p2 * np.linalg.inv(ss) - p2
    p1p = p1 / np.linalg.norm(p1) # agrees with Eq. 18 in Nishiyama

    # Calculate the displacement vector
    d1 = (ss * d2 - d2) / (p1p * d2) # agrees with Eq. 19

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

    # KSI is similarity-transformed version of J (i.e., expressed in martensite basis)
    J =  rr1 * betamx * rr2.T
    gamma = cm.A / np.tile(np.linalg.norm(cm.A, axis=1), [3, 1]).T
    KSI = gamma * J * gamma.T

    # individual ksi angles
    ksi1 = np.arccos(KSI[0, 0]) * 180.0 / np.pi
    ksi2 = np.arccos(KSI[1, 1]) * 180.0 / np.pi
    ksi3 = np.arccos(KSI[2, 2]) * 180.0 / np.pi

    # TODO: Generate full set of variant matrices
    OR = (a / a0) * cm * np.linalg.inv(ss)
    # The above should be equivalent to the T matrix in Kitahara

    # Check: With the a/a0 in front, it gives a normalized matrix
    # print(np.linalg.det(OR))

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
