# -*- coding: utf-8 -*-

def returnPhaseNames(filepath):
    '''
    Returns the phase names within a TSL Grain File Type 1
    
    Parameters
    ----------
    filepath : str
        path to Grain File type 1
    
    Returns
    -------
    phasenames : list
        phase names    
    '''

    # open the file and separate out the header from the data
    dat=[];
    with open(filepath,'r') as f:
        for line in f:
            if line[0] != '#':
                dat.append(line.strip())
    
    phasenames = []
    for line in dat:
        s = line.split()
        name = ' '.join(s[10:])
        if name not in phasenames:
            phasenames.append(name)
    
    return phasenames
    

def loadGFtype1(filepath, phasedict, notation='international'):
    '''
    Parameters
    ----------
    filepath : str
        path to Grain File type 1
    
    phasedict : dict
        dictionary for interpreting phase results
    
    Returns
    -------
    ebsd : obj of ebsd class
    
    Notes
    -----
    TSL Grain File Type 1 includes most of the data needed to interpret an
    EBSD scan, including the euler angles and xy locations for each point.
    However, the phase is not listed as an integer. Rather, it is a string for
    the phase name (which may contain spaces). This is impossible to interpret
    without 
    '''
    
    import numpy as np
    from cryspy import ebsd
    from cryspy import rot
    from cryspy import xtal
    import righthand as rh
    
    # open the file and separate out the header from the data
    hdr=[];dat=[];
    with open(filepath,'r') as f:
        for line in f:
            if line[0]=='#':
                hdr.append(line.strip()[1:].strip()) # remove #'s and whitespace
            else:
                dat.append(line.strip())
        
       
    # We should obtain the grid, the x step, the y step, the numbers of 
    # columns, and the number of rows from the x and y data themselves, 
    # because WE CAN... just in case the header was copied and pasted 
    # from a different file. Then compare to what the header says and 
    # print something if there is a discrepancy
    
    pb = rh.progbar(finalcount=np.size(dat), message='LOADING EBSD DATA')
    # Parse the data into lists

    eul1=np.zeros(np.shape(dat))
    eul2=np.zeros(np.shape(dat))
    eul3=np.zeros(np.shape(dat))
    x=np.zeros(np.shape(dat))
    y=np.zeros(np.shape(dat))
    iq=np.zeros(np.shape(dat))
    ci=np.zeros(np.shape(dat))
    fit=np.zeros(np.shape(dat))
    gid=np.zeros(np.shape(dat))
    edge=np.zeros(np.shape(dat))
    phasetmp=[]
        
    ndex = 0
    for line in dat:
        ls = line.split()
        phasetmp.append(' '.join(ls[10:]))
        
        s = [float(item) for item in ls[0:10]]
        
        eul1[ndex] = s[0]
        eul2[ndex] = s[1]
        eul3[ndex] = s[2]
        x[ndex] = s[3]
        y[ndex] = s[4]
        iq[ndex] = s[5]
        ci[ndex] = s[6]
        fit[ndex] = s[7]
        gid[ndex] = s[8]
        edge[ndex] = s[9]

        pb.update(ndex)  
        ndex += 1    

    numangcols = len(s)
           
    phasedict2 = {}
    index = 0
    for item in [key for key in phasedict.keys()]:
        phasedict2[item] = index
        index += 1

    s = np.zeros(np.shape(dat))
    phase=np.zeros(np.shape(dat))    
    index = 0
    for item in phasetmp:
        phase[index] = phasedict2[item]
        s[index] = xtal.interpret_point_group_name(phasedict[item], notation=notation)
        index += 1 
    
    # Start working with what we've interpreted now
    o = xtal.orientation(quaternions = rot.quat.from_bunge(rot.bunge(eul1,eul2,eul3)),\
                    pointgroupnumbers = np.squeeze(s.astype(np.int)))
    pb.update(-1)
    
    ebsd_data = ebsd(orientations=o, x=x, y=y, phaseid=phase)
    ebsd_data.calc_scanstepsize()
    ebsd_data.prepare_for_plotting()
    
    # Add other fields. Note that some fields are TSL-specific, and that two of
    # the fields (fit and ...) have not been added yet.
    ebsd_data.iq = iq
    ebsd_data.ci = ci
    
    ebsd_data._anghdr = u''.join(['# {0}\n'.format(line) for line in (hdr)])
    ebsd_data._originalheaderdata = hdr
    ebsd_data._oimversion_numangcols = numangcols
    ebsd_data._original_phase = phase - 1
    
    return ebsd_data
