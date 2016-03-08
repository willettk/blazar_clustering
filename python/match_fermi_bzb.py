import sys
import argparse
import numpy as np
from numpy.lib import recfunctions

from itertools import izip as zip, count

from astropy import units as u
from astropy.coordinates import ICRS
from astropy.table import join, Table
from scipy.stats.stats import pearsonr
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#from astroquery.simbad import Simbad
#from astroquery.ned import Ned

import catalog_utils as utils

class CatalogMatcher(object):
    """This class contains functions for matching
    a catalog of blazars with Fermi-LAT 3FGL catalog members
    and extracting some fun data...
    """
    def matchCats(self,bzb,fermi):
        '''
        matches 3FGL with Kyle's catalog;
        returns catalog index of Fermi source matched with Kyle's
        '''
        #fermi_tab = Table.read(fermi,format='fits')
    
        #first check if the fermi name is 
    
        name_fermi = fermi['ASSOC1']
        name_fermi2 = fermi['ASSOC2']
        name_bzb = bzb['bname']


        #print name_bzb[1]
        #ned_names = Ned.query_refcode("2010AJ....139..390P")#plotkin
        #ned_names = Ned.query_refcode("2009MNRAS.397.1713C")#Chen
        ''' 
        ned_names = {}
        for name in name_bzb:
            try:
                ned_query = Ned.query_object(name)
                ned_names[name] = ned_query["Object Name"][0]
                #ned_names.append(ned_query["Object Name"])
            except:
                continue
                #ned_names[name] = "NO NED"
        '''      
                
        #wt = ned_names['Object Name']
        #wt = ned_names['Object Name', 'RA(deg)', 'DEC(deg)','Redshift']
        #wt.write("plotkin_objnames.csv")
        #print ned_names
        ''' 
        ned_names_fermi = {}
        #simbad_names_fermi = [Simbad.query_objectids(name) for name in name_fermi]
        for name in name_fermi:
            try:
                ned_query = Ned.query_object(name)
                ned_names_fermi[name] = ned_query["Object Name"][0]
            except:
                continue
                #ned_names_fermi[name] = "NO NED"
        '''
        #print ned_names_fermi
        #print simbad_fermi
        #theta = Simbad.query_object(name1[0])
        #simbad_assoc = simbad_names["ID"]

        ra_fermi = fermi['RAJ2000']
        dec_fermi = fermi['DEJ2000']
        ra_bzb = bzb['ra']
        dec_bzb = bzb['dec']
        
        ned_names_bzb = bzb['nedname']
        ned_names_fermi = fermi['nedname']
        #print ned_names_bzb
        #print ned_names_fermi



        ned_name_idx = []
        for name in ned_names_bzb:
            if name == '---bzcat---':
                continue
                #ned_name_idx.append(None)
            else:
                try:
                    ind = [i for i, j in zip(count(), ned_names_fermi) if j == name]
                    ned_name_idx.append(ind[0])
                except IndexError:
                    continue
                    #ned_name_idx.append(None)

                #    ned_name_idx.append(ned_names_fermi.index("---fermi---"))
                    #ned_name_idx.append(ned_names_fermi.index(name))
                #except:
                #    ned_name_idx.append(None)

        print ned_name_idx
    
        #does the coordinate matching
        c_fermi = ICRS(ra_fermi, dec_fermi, unit=(u.degree,u.degree))
        c_bzb = ICRS(ra_bzb, dec_bzb, unit=(u.degree,u.degree))
    
        idx, d2d, d3d = c_bzb.match_to_catalog_sky(c_fermi)
        print len(idx),len(ned_names_bzb)
    
        matches=c_fermi[idx]
    
        c3fgl_name = 0
        cfermi_ned = 0
        
        #blah = join(bzb,fermi_tab,join_type='left')
        ''' 
        for name in ned_names_bzb:
            if name in name_fermi or name in name_fermi2:
                c3fgl_name+=1
                print "found bzb match to fermi_assoc_field:", name
            elif name in ned_names_fermi:
                cfermi_ned+=1
                print "found bzb match to fermi-ned query:", name
            else:
                print "no match found"
        print "Total number of bzb objects:", len(ned_names_bzb)
        print "Matched %d names using 3fgl and %d with NED queries." %(c3fgl_name, cfermi_ned)
        '''
        
        dra = (matches.ra - c_bzb.ra).arcmin
        ddec = (matches.dec - c_bzb.dec).arcmin
        ddist = (dra*dra + ddec*ddec)**0.5
    
        smaja95 = fermi['Conf_95_SemiMajor']
        smina95 = fermi['Conf_95_SemiMinor']
        pos95 = fermi['Conf_95_PosAng']
    
        distmask =  [1 if d[0] < d[1]*60. else 0 for d in zip(ddist,smaja95[idx])]
        match_idx = np.where(distmask)


        #index_dict = dict((value, ind) for ind, value in enumerate(idx))
        #matched_names_idx = [index_dict[name] for name in ned_name_idx]
        #print len(matched_names_idx)

    
        #can keep this the way we're doing but put in a flag for match confidence
        return [ fermi[idx][match_idx], bzb[match_idx], ddist]
    
    def writeMatchedCat(self,fermi,bzb,dist,outname='matched3fgl.csv'):
    
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
    
        np.savetxt(outname, bzb, delimiter=",",fmt="%s",comments='', header='bname,btype,ra,dec,catalogue,nt,nb,z,fieldsize,countingmag,bgb,bgb_err,ferminame,fermidist,fermivar,index,class,sig')
    
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
        
        plt.show()
    
    
def main():

    parser = argparse.ArgumentParser(description="Code that matches Kyle's bzb catalog with the Fermi 3FGL catalog.")
    parser.add_argument('--write','-w',action='store_true',help='Enables writing out of the matched catalog.')
    parser.add_argument('--plot','-p',action='store_true',help='Enables displaying plots.')
    parser.add_argument('--supplement','-s',action='store_true',help='Supplementing catalogs with NED names and exits.')
    args = parser.parse_args()
    #for the first pass at this add a column of alternative names in the Fermi-LAT catalog

    if args.supplement:
        print "--- supplement option enabled ---"
        utils.supplementBZBWithNED('../data/blazars_bgb.csv')
        utils.supplementFermiWithNED('../data/gll_psc_v16.fit')
        sys.exit("Added a column with NED names to 3FGL and bzcat catalog files. Exiting...")
    
    #reading in catalogs
    try:
        #fermi_dat = Table.read('../data/gll_psc_v16_sup.fit',format='fits')
        fermi_dat = utils.readFermi('../data/gll_psc_v16_sup.fit')
        bzb_dat = utils.readBZB('../data/blazars_bgb_sup.csv')
    except:
        sys.exit("ERROR: Run with the supplement option first:  'python match_fermi_bzb.py --supplement'")
    #fermi_dat = utils.readFermi('../data/gll_psc_v16.fit')
    #bzb_dat = utils.readBZB('../data/blazars_bgb.csv')

    matcher = CatalogMatcher()
    
    #tbl = utils.reformatFermiFits(fermi_dat)
    #print Table.read(tbl)
    #sys.exit()

    #matching and reading of results
    match_results = matcher.matchCats(bzb_dat,fermi_dat)
    fermi_matched = match_results[0]
    bzb_matched = match_results[1]
    match_dist = match_results[2] #in arcsec

    #write out matched catalog
    if args.write == True:
        matcher.writeMatchedCat(fermi_matched,bzb_matched,match_dist)
    
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

    if args.plot == True:
        matcher.plotParamRelations(bzb_matched[temp_idx][z_cut][sig_cut], fermi_matched[temp_idx][z_cut][sig_cut])
    
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
