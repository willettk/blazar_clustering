import os
import numpy as np
from numpy import recfromcsv
from numpy.lib import recfunctions
import astropy.io.fits as pf
from astroquery.ned import Ned

#def readBZB(bzb_filename):
#    '''reads Kyle's catalog, expecting the csv file to be in this directory'''
#    print '** Reading Kyle\'s catalog'
#    dat = recfromcsv(bzb_filename, unpack=True)
#    return dat

def readCat(filename):
    '''reads fits catalog and returns the data table'''
    f = pf.open(filename)
    dat = f[1].data
    return dat

def queryNED(names,cat):
    '''
    takes a list of names, queries them in NED and returns the 'Object Name'
    from the NED query results.
    '''

    ned_names = []

    if cat == 'bzb':
        failstr = '---bzcat---'
    elif cat == 'fermi':
        failstr = '---fermi---'
    else:
        failstr = '---'

    for name in names:
        try:
            ned_query = Ned.query_object(name)
            ned_names.append(ned_query["Object Name"][0])
        except:
            ned_names.append(failstr)

    return ned_names


def supplementBZBWithNED(bzb_filename, overwrite = False):
    '''
    Adds a column to Kyle's catalog with NED names for the bzcat objects
    '''
    if overwrite:
        out_name = bzb_filename
    else:
        orig_file = os.path.splitext(bzb_filename)
        out_name = orig_file[0]+'_sup'+orig_file[1]

    bzb_dat = readCat(bzb_filename)

    name_col = queryNED(bzb_dat['bname'],'bzb')

    bzb_dat = recfunctions.append_fields(bzb_dat, 'nedname', name_col, usemask=False)
    
    np.savetxt(out_name, bzb_dat, delimiter=",",fmt="%s",comments='', header='bname,btype,ra,dec,catalogue,nt,nb,z,fieldsize,countingmag,bgb,bgb_err,nedname')

    return None


def supplementFermiWithNED(fermi_filename, overwrite = False):
    '''
    Adds a column to 3FGL catalog with NED names for the 3FGL associations
    '''

    if overwrite:
        out_file = fermi_filename
    else:
        orig_file = os.path.splitext(fermi_filename)
        out_file = orig_file[0]+'_sup'+orig_file[1]


    fermi_dat = readCat(fermi_filename)
    fermi_cols = fermi_dat.columns

    name_col = queryNED(fermi_dat['ASSOC1'],'fermi')

    ned_names = pf.ColDefs([pf.Column(name="nedname",format='A20',array=name_col)])
    hdu = pf.BinTableHDU.from_columns(fermi_cols + ned_names)
    hdu.writeto(out_file)

    return None
