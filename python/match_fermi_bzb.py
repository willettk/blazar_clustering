import pyfits as pf
import numpy as np
from numpy.lib import recfunctions
from numpy import recfromcsv
from astropy.coordinates import ICRS
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from scipy.optimize import curve_fit



#fermi_filename = 'gll_psc_v08.fit'
fermi_filename= '../data/gll_psc_v14.fit'
bgb_filename = '../data/blazars_bgb.csv'

#Getting data from 2FGL
dat = pf.getdata(fermi_filename)
hdr = pf.getheader(fermi_filename)

ra_fermi = dat['RAJ2000']
#print ra_fermi
dec_fermi = dat['DEJ2000']
#print dec_fermi
var = dat['Variability_Index']
sig = dat['Signif_Avg']
sindex = dat['Spectral_Index']
flux=dat['Energy_Flux100']
flux_err=dat['Unc_Energy_Flux100']
peak_flux=dat['Flux_Peak']
peak_flux_err=dat['Unc_Flux_Peak']

fermi_name = dat['Source_Name']
temp_n = np.where(fermi_name == '3FGL J1743.9+1934')
#print var[temp_n]

#name1 = dat['ASSOC1']
#name2 = dat['ASSOC2']

smaja95 = dat['Conf_95_SemiMajor']
smina95 = dat['Conf_95_SemiMinor']
pos95 = dat['Conf_95_PosAng']

print max(smaja95), min(smaja95)

class1 = dat['CLASS1']
#class2 = dat['CLASS2']

#Getting data from bgb
bgb_dat = recfromcsv(bgb_filename, unpack=True)


bname = bgb_dat['bname']
btype = bgb_dat['btype']
ra_bgb = bgb_dat['ra']
dec_bgb = bgb_dat['dec']
catalog = bgb_dat['catalogue']
nt = bgb_dat['nt']
nb = bgb_dat['nb']
z = bgb_dat['z']
fieldsize = bgb_dat['fieldsize']
countingmag = bgb_dat['countingmag']
bgb = bgb_dat['bgb']
bgb_err = bgb_dat['bgb_err']


c_fermi = ICRS(ra_fermi, dec_fermi, unit=(u.degree,u.degree))
c_bgb = ICRS(ra_bgb, dec_bgb, unit=(u.degree,u.degree))

idx, d2d, d3d = c_bgb.match_to_catalog_sky(c_fermi)

matches=c_fermi[idx]

dra = (matches.ra - c_bgb.ra).arcmin
ddec = (matches.dec - c_bgb.dec).arcmin
ddist = (dra*dra + ddec*ddec)**0.5

distmask =  [1 if d[0] < d[1]*60. else 0 for d in zip(ddist,smaja95[idx])]
match_idx = np.where(distmask)
#ddist.reshape(778,1)

bgb_dat = bgb_dat[match_idx]
bgb_dat = recfunctions.append_fields(bgb_dat, 'fermi_name', fermi_name[idx][match_idx], usemask=False)
bgb_dat = recfunctions.append_fields(bgb_dat, 'distance', ddist[match_idx], usemask=False)
bgb_dat = recfunctions.append_fields(bgb_dat, 'variability', var[idx][match_idx], usemask=False)
bgb_dat = recfunctions.append_fields(bgb_dat, 'index', sindex[idx][match_idx], usemask=False)
bgb_dat = recfunctions.append_fields(bgb_dat, 'class', class1[idx][match_idx], usemask=False)
bgb_dat = recfunctions.append_fields(bgb_dat, 'sig', sig[idx][match_idx], usemask=False)

np.savetxt('matched3fgl.csv', bgb_dat, delimiter=",",fmt="%s", header='bname,btype,ra,dec,catalogue,nt,nb,z,fieldsize,countingmag,bgb,bgb_err,ferminame,fermidist,fermivar,index,class,sig')


temp_idx = np.where(abs(bgb[match_idx]) >= 0)
z_cut = np.where(z[match_idx][temp_idx] > 0.15)
sig_cut = np.where(sig[idx][match_idx][temp_idx][z_cut] > 10.)
print len(np.where(bgb[match_idx][temp_idx][z_cut] > 0)[0])
print len(np.where(bgb[match_idx][temp_idx][z_cut] < 0)[0])
print np.mean(bgb[match_idx][temp_idx][z_cut])
print np.std(bgb[match_idx][temp_idx][z_cut])
#print z_cut



x = pearsonr(bgb[match_idx][temp_idx][z_cut][sig_cut], var[idx][match_idx][temp_idx][z_cut][sig_cut])
print "Pearson's correlation coefficient:", x[0], "; p-value:",x[1]

#print np.corrcoef(bgb[match_idx], var[match_idx])[0, 1]

#plt.xlabel('Clustering')
#plt.xlabel('redshift')
#plt.ylabel('Spectral index')
plt.xlabel('Flux')
plt.ylabel('Variability index')
plt.yscale('log')

#plotting line at variabilty cutoff
var_line = np.arange(-1000,1000,0.1)

plt.axis([0, 3e-10, 10, 20000])#var vs flux
#plt.axis([-450, 450, 20, 20000])#var vs clustering
#plt.errorbar(bgb[match_idx][temp_idx][z_cut],sindex[idx][match_idx][temp_idx][z_cut],fmt='.')
#plt.errorbar(z[match_idx][temp_idx][z_cut][sig_cut],var[idx][match_idx][temp_idx][z_cut][sig_cut],fmt='.')
#plt.errorbar(bgb[match_idx][temp_idx][z_cut][sig_cut],var[idx][match_idx][temp_idx][z_cut][sig_cut],xerr=bgb_err[match_idx][temp_idx][z_cut][sig_cut],fmt='.')
plt.errorbar(flux[idx][match_idx][temp_idx][z_cut][sig_cut],var[idx][match_idx][temp_idx][z_cut][sig_cut],fmt='.')#,xerr=peak_flux_err[idx][match_idx][temp_idx][z_cut][sig_cut],fmt='.')
#plt.errorbar(var_line, const_func(var_line))
plt.axhline(y=72.44,color='g')
#plt.errorbar(z[match_idx],var[match_idx],fmt='.')
plt.show()

#all_data = np.concatenate((bgb_dat, ddist), 1)

#all_data = np.hstack((bgb_dat,ddist))
#print all_data

#bname, btype, ra_bgb, dec_bgb, catalog, nt, nb, z , fieldsize ,countingmag bgb bgb_err 


