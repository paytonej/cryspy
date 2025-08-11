# -*- coding: utf-8 -*-

def loadebsd(filepath):
    
    from cryspy import cio
    '''
    Load EBSD data based on filename extension.
    '''
    test = filepath.split('.')[-1]
    if test == 'ctf':
        ovdat = cio.loadctf(filepath)
    elif test == 'ang':
        ovdat = cio.loadang(filepath)
    elif test == 'osc':
        ovdat = cio.loadosc(filepath, create_ebsd_object = True, \
                        ang_output = False, ang_output_filename = None)
    else:
        'EBSD data can be read in ctf, ang, and osc formats.'
    
    return ovdat
        
    