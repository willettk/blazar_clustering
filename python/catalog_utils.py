import numpy as np
from numpy import recfromcsv
import astropy.io.fits as pf
from astroquery.ned import Ned

def readBZB(bzb_filename):
    '''reads Kyle's catalog, expecting the csv file to be in this directory'''
    dat = recfromcsv(bzb_filename, unpack=True)
    return dat

def readFermi(fermi_filename):
    '''reads the 3FGL catalog, expecting file to be in this directory'''
    dat = pf.getdata(fermi_filename)
    return dat

def supplementBZBWithNed():
    bzb_dat = readBZB('../data/blazars_bgb.csv')
    return "hmm"


def supplementFermiWithNed():
    fermi_dat = readFermi('../data/gll_psc_v16.fit')
    print fermi_dat 
    return "hmm"
