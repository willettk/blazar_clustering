import astropy.io.fits as pf
import numpy as np
from numpy.lib import recfunctions
from numpy import recfromcsv
from astropy.coordinates import ICRS
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.optimize import curve_fit

from astroquery.simbad import Simbad
from astroquery.ned import Ned

def const_func(x):
    return x*0+72.44

def addAltNamesToFermi(fermi_cat):
    read

def readBZB(bzb_filename):
    '''reads Kyle's catalog, expecting the csv file to be in this directory'''
    dat = recfromcsv(bzb_filename, unpack=True)
    return dat

def readFermi(fermi_filename):
    '''reads the 3FGL catalog, expecting file to be in this directory'''
    #Getting data from 2FGL
    dat = pf.getdata(fermi_filename)
    #hdr = pf.getheader(fermi_filename)
    return dat

def matchCats(bzb,fermi):
    '''
    matches 3FGL with Kyle's catalog;
    returns catalog index of Fermi source matched with Kyle's
    '''

    #first check if the fermi name is 

    name_fermi = fermi['ASSOC1']
    name_fermi2 = fermi['ASSOC2']
    name_bzb = bzb['bname']
    #print name_bzb[1]
    #ned_names = Ned.query_refcode("2010AJ....139..390P")#plotkin
    #ned_names = Ned.query_refcode("2009MNRAS.397.1713C")#Chen

    ned_names = {}
    for name in name_bzb:
        try:
            ned_query = Ned.query_object(name)
            ned_names[name] = ned_query["Object Name"]
            #ned_names.append(ned_query["Object Name"])
        except:
            ned_names[name] = "NO NED"
            
            
    #wt = ned_names['Object Name']
    #wt = ned_names['Object Name', 'RA(deg)', 'DEC(deg)','Redshift']
    #wt.write("plotkin_objnames.csv")
    #print ned_names



    ned_names_fermi = {}
    #simbad_names_fermi = [Simbad.query_objectids(name) for name in name_fermi]
    for name in name_fermi:
        try:
            ned_query =  Ned.query_object(name)
            ned_names_fermi[name] = ned_query["Object Name"]
        except:
            ned_names_fermi[name] = "NO NED"

    print ned_names_fermi
    #print simbad_fermi
    #theta = Simbad.query_object(name1[0])
    #simbad_assoc = simbad_names["ID"]
    ra_fermi = fermi['RAJ2000']
    dec_fermi = fermi['DEJ2000']
    ra_bzb = bzb['ra']
    dec_bzb = bzb['dec']

    #does the coordinate matching
    c_fermi = ICRS(ra_fermi, dec_fermi, unit=(u.degree,u.degree))
    c_bzb = ICRS(ra_bzb, dec_bzb, unit=(u.degree,u.degree))

    idx, d2d, d3d = c_bzb.match_to_catalog_sky(c_fermi)

    matches=c_fermi[idx]

    for name,ned in ned_names.iteritems():
        if ned in name_fermi or ned in name_fermi2:
            print ned


    dra = (matches.ra - c_bzb.ra).arcmin
    ddec = (matches.dec - c_bzb.dec).arcmin
    ddist = (dra*dra + ddec*ddec)**0.5

    smaja95 = fermi['Conf_95_SemiMajor']
    smina95 = fermi['Conf_95_SemiMinor']
    pos95 = fermi['Conf_95_PosAng']

    distmask =  [1 if d[0] < d[1]*60. else 0 for d in zip(ddist,smaja95[idx])]
    match_idx = np.where(distmask)

    #can keep this the way we're doing but put in a flag for match confidence
    return [ fermi[idx][match_idx], bzb[match_idx], ddist]

def writeMatchedCat(fermi,bzb,dist,outname='matched3fgl.csv'):

    fermi_name = fermi['Source_Name']

    var = fermi['Variability_Index']
    sindex = fermi['Spectral_Index']
    class1 = fermi['CLASS1']
    sig = fermi['Signif_Avg']
    #flux=fermi['Energy_Flux100']
    #flux_err=fermi['Unc_Energy_Flux100']
    #peak_flux=fermi['Flux_Peak']
    #peak_flux_err=fermi['Unc_Flux_Peak']

    #name1 = fermi['ASSOC1']
    #name2 = fermi['ASSOC2']


    bzb = recfunctions.append_fields(bzb, 'fermi_name', fermi_name, usemask=False)
    bzb = recfunctions.append_fields(bzb, 'distance', dist, usemask=False)
    bzb = recfunctions.append_fields(bzb, 'variability', var, usemask=False)
    bzb = recfunctions.append_fields(bzb, 'index', sindex, usemask=False)
    bzb = recfunctions.append_fields(bzb, 'class', class1, usemask=False)
    bzb = recfunctions.append_fields(bzb, 'sig', sig, usemask=False)

    np.savetxt(outname, bzb, delimiter=",",fmt="%s", header='bname,btype,ra,dec,catalogue,nt,nb,z,fieldsize,countingmag,bgb,bgb_err,ferminame,fermidist,fermivar,index,class,sig')

def plotParamRelations(bzb,fermi):

    #plotting clustering vs.
    fc, (ax1c, ax2c, ax3c) = plt.subplots(3,1,sharex=True)
    
    #variability
    bgb = bzb['bgb']
    bgb_err = bzb['bgb_err']
    var = fermi['Variability_Index']
    pflux = fermi['Flux_Peak']
    pflux_err=fermi['Unc_Flux_Peak']

    ax1c.set_ylabel('variability index / peak flux')
    ax1c.set_yscale("log")

    ax1c.errorbar(bgb,np.divide(var,pflux),xerr=bgb_err,fmt='.')


    #variability
    ax2c.set_ylabel('variability index')
    ax2c.set_yscale("log")

    ax2c.errorbar(bgb,var,xerr=bgb_err,fmt='.')
    ax2c.axhline(y=72.44,color='g')

    #peak flux
    ax3c.set_ylabel('Peak Flux')
    ax3c.set_xlabel('Clustering')
    ax3c.errorbar(bgb,pflux,xerr=bgb_err,yerr=pflux_err,fmt='.')

    #plt.show()


def main():

    
    #for the first pass at this add a column of alternative names in the Fermi-LAT catalog

    #reading in catalogs
    fermi_dat = readFermi('gll_psc_v16.fit')
    bzb_dat = readBZB('blazars_bgb.csv')

    #matching and reading of results
    match_results = matchCats(bzb_dat,fermi_dat)
    fermi_matched = match_results[0]
    bzb_matched = match_results[1]
    match_dist = match_results[2] #in arcsec

    #write out matched catalog
    writeMatchedCat(fermi_matched,bzb_matched,match_dist)

    
    #fermi
    var = fermi_matched['Variability_Index']
    sig = fermi_matched['Signif_Avg']
    sindex = fermi_matched['Spectral_Index']
    flux=fermi_matched['Energy_Flux100']
    flux_err=fermi_matched['Unc_Energy_Flux100']
    peak_flux=fermi_matched['Flux_Peak']
    peak_flux_err=fermi_matched['Unc_Flux_Peak']
    fermi_name = fermi_matched['Source_Name']
    #temp_n = np.where(fermi_name == '3FGL J1743.9+1934')
    #print var[temp_n]
    name1 = fermi_matched['ASSOC1']
    #name2 = dat['ASSOC2']
    class1 = fermi_matched['CLASS1']
    #class2 = dat['CLASS2']

    #Getting data from bzb
    bname = bzb_matched['bname']
    btype = bzb_matched['btype']
    catalog = bzb_matched['catalogue']
    nt = bzb_matched['nt']
    nb = bzb_matched['nb']
    z = bzb_matched['z']
    fieldsize = bzb_matched['fieldsize']
    countingmag = bzb_matched['countingmag']
    bgb = bzb_matched['bgb']
    bgb_err = bzb_matched['bgb_err']

    
    #applying additional cuts
    temp_idx = np.where(abs(bgb) >= 0) #cut on clustering 
    z_cut = np.where(z[temp_idx] > 0.15) #cut on redshift
    sig_cut = np.where(sig[temp_idx][z_cut] > 10.) #cut on fermi significance

    plotParamRelations(bzb_matched[temp_idx][z_cut][sig_cut], fermi_matched[temp_idx][z_cut][sig_cut])
    
    print 'num of obj less than 0 clustering:',len(np.where(bgb[temp_idx][z_cut] > 0)[0])
    print 'num of obj more than 0 clustering:',len(np.where(bgb[temp_idx][z_cut] < 0)[0])
    print 'printing mean of clustering:',np.mean(bgb[temp_idx][z_cut])
    print 'printing std of clustering:',np.std(bgb[temp_idx][z_cut])

    x = pearsonr(bgb[temp_idx][z_cut][sig_cut], var[temp_idx][z_cut][sig_cut])
    print "Pearson's correlation coefficient:", x[0], "; p-value:",x[1]

    #print np.corrcoef(bgb[match_idx], var[match_idx])[0, 1]

    #plt.xlabel('Clustering')
    #plt.xlabel('redshift')
    #plt.ylabel('Spectral index')
    #plt.xlabel('Peak Flux')
    #plt.ylabel('Variability index')
    #plt.yscale('log')

    #plt.axis([0, 1e-6, 10, 20000])#var vs flux
    #plt.axis([0.1, 0.8, 10, 20000])#var vs flux
    #plt.axis([-450, 450, 20, 20000])#var vs clustering
    #plt.errorbar(bgb[match_idx][temp_idx][z_cut],sindex[idx][match_idx][temp_idx][z_cut],fmt='.')
    #plt.errorbar(z[temp_idx][z_cut][sig_cut],var[temp_idx][z_cut][sig_cut],fmt='.')
    #plt.errorbar(bgb[match_idx][temp_idx][z_cut][sig_cut],var[idx][match_idx][temp_idx][z_cut][sig_cut],xerr=bgb_err[match_idx][temp_idx][z_cut][sig_cut],fmt='.')
    #plt.errorbar(peak_flux[idx][match_idx][temp_idx][z_cut][sig_cut],var[idx][match_idx][temp_idx][z_cut][sig_cut],fmt='.')#,xerr=peak_flux_err[idx][match_idx][temp_idx][z_cut][sig_cut],fmt='.')
    #plt.errorbar(flux[idx][match_idx][temp_idx][z_cut][sig_cut],var[idx][match_idx][temp_idx][z_cut][sig_cut],fmt='.')#,xerr=peak_flux_err[idx][match_idx][temp_idx][z_cut][sig_cut],fmt='.')
    #plt.errorbar(var_line, const_func(var_line))
    #plt.axhline(y=72.44,color='g')
    #plt.errorbar(z[match_idx],var[match_idx],fmt='.')
    #plt.axis([0,500,0,100])
    #plt.hist(smaja95[idx])
    #plt.show()

#all_data = np.concatenate((bgb_dat, ddist), 1)

#all_data = np.hstack((bgb_dat,ddist))
#print all_data

#bname, btype, ra_bgb, dec_bgb, catalog, nt, nb, z , fieldsize ,countingmag bgb bgb_err 


if __name__ == '__main__':
    main()
