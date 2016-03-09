# Goal: port the computation of B into Python.

import numpy as np

from astropy.cosmology import FlatLambdaCDM,WMAP7,WMAP9
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

from scipy import integrate

from matplotlib import pyplot as plt

import warnings 
warnings.filterwarnings("ignore",category=RuntimeWarning,append=True)

# Default cosmology for final paper should be WMAP9
c = WMAP7

blazardir = '/Users/willettk/Astronomy/Research/blazars'
btype_dict = {'bllac':'BL Lac','fsrq':'FSRQ','uncertain':'Blazar (uncertain type)','candidate':'BL Lac candidate','bllac_gal':'BL Lac (galaxy dominated)'}

def schechter_mag(phi_star,mr_star,alpha,m):

    # Differential form of the Schechter luminosity function in absolute magnitudes 
    sch_mag = 0.4 * np.log(10)  * phi_star* (10.**(0.4 * (mr_star - m)))**(alpha + 1.) * np.exp(-1. * (10.**(0.4*(mr_star - m)))) 

    return sch_mag

def get_luminosity_function(lf,z):

    # Luminosity function parameters
    
    lfd_static = {
            'dr1': {
            		'phi_star' : 1.49e-2 * c.h**3,	# Blanton et al. (2003)
            		'Mr_star' : -20.44+5*np.log10(c.h),	# SDSS, pre-DR1
            		'alpha' : -1.05			# 
            	    },
            'dr6': {
            		'phi_star' : 0.90e-2 * c.h**3,	# Montero-Dorta & Prada (2006)
            		'Mr_star' : -20.73+5*np.log10(c.h),	# SDSS DR6
            		'alpha' : -1.23			# 
                    },
            '6dfgs': {
            		'phi_star' : 10.**(-2.081)	* c.h**3,	# Jones et al. (2006)
            		'Mr_star' : -20.98+5*np.log10(c.h),	# 6dFGS
            		'alpha' : -1.21			# 
                    }
          }

    if lf in lfd_static.keys():
        phi_star = lfd_static[lf]['phi_star']
        Mr_star = lfd_static[lf]['Mr_star']
        alpha = lfd_static[lf]['alpha']

    if lf == 'dr6ages':          # Use DR6 for low-z, AGES for increasing redshift
        if z < 0.15:
            phi_star = lfd_static['dr6']['phi_star']
            Mr_star = lfd_static['dr6']['Mr_star']
            alpha = lfd_static['dr6']['alpha']
        
        else:
            # Evolving LF from Table 4 in Cool et al.
        
            zarr = np.array([0.15,0.25,0.35,0.45,0.55,0.65])
            phi_star_arr = np.array([0.0159,0.0152,0.0124,0.0144,0.0108,0.0105]) * c.h**3
            Mr_star_arr = np.array([-20.58,-20.81,-20.81,-20.99,-21.29,-21.38]) + 5*np.log10(c.h)

            phi_star = phi_star_arr[abs(zarr-z).argmin()]
            Mr_star = Mr_star_arr[abs(zarr-z).argmin()]
            alpha = -1.05
            
    if lf == 'gama':
        phi_star = 0.90 * c.h**3/100. 	# Loveday et al. (2012)
        Mr_star = -20.73+5*np.log10(c.h)	# GAMA
        alpha = -1.26			# 
        
        p = 1.6				# step-wise max. likelihood values from Table 5
        q = 0.2
        
        z0 = 0.1			# normalization redshift of the LF
        
        phi_star = phi_star * 10.**(0.4 * p * (z - z0))	# Evolution according to parameterization of Lin et al. (1999)
        Mr_star = Mr_star - q*(z - z0)
        alpha += z

    if lf == 'ramosalmeida':
    
        alpha = -1.3
        
        # Evolving LF from Table 4
        
        zarr = (np.arange(5)+1)/5.
        phi_star_arr  = np.array([0.0038,0.0037,0.0035,0.0033,0.0031])
        if z > 0.4:
            Mr_star_arr = np.array([-21.43,-22.08,-22.77,-22.62,-22.87])
        else:
            Mr_star_arr = np.array([-21.76,-22.47,-23.27,-23.59,-23.84])

        phi_star = phi_star_arr[abs(zarr-z).argmin()]
        Mr_star = Mr_star_arr[abs(zarr-z).argmin()]
        alpha = -1.05

    return phi_star, Mr_star, alpha

def get_bgb(nt, nb, field_size, z, counting_mag, lf = 'dr6ages', verbose=False):

    # Angular covariance amplitude (Yee & Green 1987)
    
    gamma = 1.77			# Power-law index of galaxy-galaxy covariance relationship
    theta = field_size.to(u.radian)
    A_gB = (nt - nb) / nb * (3. - gamma) / 2. * theta**(gamma-1)
    
    i_gamma = 3.87
    da = c.angular_diameter_distance(z)
    
    if lf == 'ramosalmeida':
        lowerlim = counting_mag - 1 
    else:
        lowerlim = -np.inf

    phi_star, Mr_star, alpha = get_luminosity_function(lf,z)
    psi,psi_err = integrate.quad(lambda m: schechter_mag(phi_star,Mr_star,alpha,m), lowerlim, counting_mag)
    
    n_net = (nt - nb)
    a_theta = np.pi * theta**2

    # Muzzin et al. (2007)
    B_gB = n_net * (3. - gamma) * da**(gamma - 3.) * theta**(gamma - 1.) / (2. * a_theta * i_gamma * np.array(psi))
    
    # Wold et al, Best et al, Zauderer et al. method
    #ng = 1e8		# Taken from Best paper - surface density of all galaxies per sr
    #B_gB = A_gB * ng / i_gamma * da**(gamma - 3.) / psi	# YG87 have the exponent wrong.
    
    # Compute uncertainty on bgb (Yee & Lopez-Cruz 1999)
    
    if B_gB != 0:
        err_B_gB = B_gB * np.sqrt(abs(nt - nb) + 1.3**2 * nb) / abs(nt - nb)
    else:
        err_B_gB = 0 * B_gB.unit

    # Print out relevant variables - test against IDL routine bgb.pro
    
    if verbose:
        print '\nn_net: %i' % n_net
        print 'gamma: %.2f' % gamma
        print 'da: %.2f' % da.value
        print 'theta: %.2e' % theta.value
        print 'a_theta: %.2e' % a_theta.value
        print 'i_gamma: %.2f' % i_gamma
        print 'counting mag: %.2f\n' % counting_mag
        print 'phi_star: %.3e' % phi_star
        print 'Mr_star: %.3f' % Mr_star
        print 'alpha: %.3f\n' % alpha
        print 'h: %.4f\n' % c.h
        print 'psi: %.2e\n' % psi
    
    return B_gB, err_B_gB

def bstats(data):

    z = data['z'][0]

    # Angular size of the field to count neighbors
    field_size = c.arcsec_per_kpc_comoving(z) * 500. * u.kpc        # arcsec per 500 kpc at redshift z
    field_size_arcmin = field_size.to(u.arcmin)
    
    # Limit on magnitude down to which a galaxy is counted as a neighbor
    Mr_star = -21.1837		# Blanton et al. (2003) - shouldn't this change with distance?
    abs_mstar = Mr_star + 2
    
    sdss_lim = 22.2			# SDSS apparent limiting magnitude
    lumdist = c.luminosity_distance(z)
    sdss_lim_abs = sdss_lim - 5*np.log10(lumdist.to(u.pc).value) + 5
    
    counting_mag = min(sdss_lim_abs,abs_mstar)
        
    # Count number of galaxies within field_size down to counting mag
    nt = np.sum(
                (data['sep_arcmin'] <= field_size_arcmin.value) & 
                (data['absmagR'] <= counting_mag)
                )

    # For galaxies above a redshift cut, use counts from an equal-sized external radius
    # in the 10 arcmin SDSS search as background

    zlim_control = 0.0613
    if z > zlim_control:

        # Results should check that I'm not near the edge of SDSS for any of these ...

        nb = np.sum(
                    (data['sep_arcmin'] > field_size_arcmin.value) & 
                    (data['sep_arcmin'] <= field_size_arcmin.value * np.sqrt(2.)) & 
                    (data['absmagR'] <= counting_mag)
                    )

    # Otherwise, use galaxies from separate control field centered 1 degree to north
    else:

        # Open separate data file and match to current blazar
        with fits.open('%s/newdata/output/sdss/blazars_unique_lowz_control_10arcmin_sdss_willettk.fit' % blazardir) as f:
            cdata = f[1].data

        cdata_match = cdata[cdata['bname']==data['bname'][0]]

        nb = np.sum(
                    (cdata_match['sep_arcmin'] <= field_size_arcmin.value) & 
                    (cdata_match['absmagR'] <= counting_mag)
                    )

        '''
        print 'Using control field as background'
        print 'nt = %i, nb = %i' % (nt,nb)
        '''
        
    return nt, nb, field_size, z, counting_mag

def load_blazar_data(zcut=True):

    with fits.open('%s/newdata/output/sdss/blazars_unique_all_10arcmin_sdss.fits' % blazardir) as f:
        data = f[1].data

    if zcut:
        data = data[(data['z'] < 0.75) & (data['z'] >= 0.043)]     # Reduces data from ~2M to ~1M rows

    return data
    
def bgb_blazars(data):

    blazars = set(data['bname'])

    col_bname = []
    col_btype = []
    col_ra = []
    col_dec = []
    col_catalogue = []
    col_nt = []
    col_nb = []
    col_z = []
    col_field_size = []
    col_counting_mag = []
    col_bgb = []
    col_bgb_err = []

    for b in blazars:
    
        name_match = (data['bname'] == b)
        ngals = np.sum(name_match)
    
        if ngals >= 10:

            nt,nb,field_size,z,counting_mag = bstats(data[name_match])
            b_val,b_err = get_bgb(nt, nb, field_size, z, counting_mag, lf = 'dr6ages')
            #bd[b] = {'nt':nt,'nb':nb,'field_size':field_size,'z':z,'counting_mag':counting_mag,'bgb':b_val,'bgb_err':b_err}
            
            col_bname.append(b)
            col_btype.append((data[name_match][0]['btype']).strip())
            col_ra.append((data[name_match][0]['ra']))
            col_dec.append((data[name_match][0]['dec']))
            col_catalogue.append((data[name_match][0]['catalogue']).strip())
            col_nt.append(nt)
            col_nb.append(nb)
            col_z.append(z)
            col_field_size.append(field_size.value)
            col_counting_mag.append(counting_mag)
            col_bgb.append(b_val.value)
            col_bgb_err.append(b_err.value)

            #print 'B_gB = (%.2f +- %.1f) Mpc^1.77' % (b_val.value,b_err.value)
        else:
            pass
            #print '%s had only %i galaxies within r_proj = 500 kpc' % (b,ngals)
    
    #print good,bad
    btable = Table([col_bname,col_btype,col_ra,col_dec,col_catalogue,col_nt,col_nb,col_z,col_field_size,col_counting_mag,col_bgb,col_bgb_err], names=('bname','btype','ra','dec','catalogue','nt','nb','z','field_size','counting_mag','bgb','bgb_err'), meta={'name': 'All blazars in SDSS'})

    return btable

def plot_zhist(btable,savefig=False):

    fig = plt.figure(1,figsize=(9,7))
    fig.clf()
    ax = fig.add_subplot(111)
    ax.hist(btable['z'],bins=10,range=(0,0.8),histtype='step',lw=3,color='black',label='All blazars')

    color_cycle = ['377EB8', 'E41A1C', '4DAF4A', '984EA3', 'FF7F00'][::-1]

    for bt in set(btable['btype']):
        clr = '#'+color_cycle.pop().lower()
        #ax.hist(btable['z'][btable['btype'] == bt],bins=10,range=(0,0.8),histtype='stepfilled',alpha=0.3,label=btype_dict[bt])
        ax.hist(btable['z'][btable['btype'] == bt],bins=10,range=(0,0.8),histtype='step',alpha=1.0,lw=2,label=btype_dict[bt],color=clr)
        ax.axvline(np.mean(btable['z'][btable['btype'] == bt]),linestyle='--',lw=1,color=clr)

    ax.set_xlabel('Redshift',fontsize=20)
    ax.set_ylabel('Count',fontsize=20)

    # Now add the legend with some customizations.
    legend = ax.legend(loc='upper right', shadow=True)
    
    # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    frame = legend.get_frame()
    frame.set_facecolor('0.90')
    
    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize('large')
    
    for label in legend.get_lines():
        label.set_linewidth(1.5)  # the legend line width
        plt.show()

    if savefig:
        fig.savefig('%s/paper/figures/blazars_zhist.pdf' % blazardir)
    else:
        plt.show()

    return None

def plot_bgb_hist(btable,savefig=False):

    fig,axarr = plt.subplots(2,3,sharex=True,sharey=True,figsize=(14,8),num=1)
    btypes = set(btable['btype'])
    print '\nAll blazars: median B_gB = %i' % np.median(btable['bgb'])

    color_cycle = ['377EB8', 'E41A1C', '4DAF4A', '984EA3', 'FF7F00'][::-1]

    for ax,bt in zip(axarr.flatten(),btypes):
        clr = '#'+color_cycle.pop().lower()
        ax.hist(btable['bgb'],bins=20,range=(-1000,1000),histtype='step',lw=2,color='black')
        ax.hist(btable['bgb'][btable['btype'] == bt],bins=20,range=(-1000,1000),histtype='stepfilled',color=clr)
        ax.set_ylim(0,200)
        ax.set_title(btype_dict[bt])
        print '%s blazars: median B_gB = %i' % (bt,np.median(btable['bgb'][btable['btype'] == bt]))

        ax.vlines(np.median(btable['bgb']),ax.get_ylim()[0],ax.get_ylim()[1],lw=2,color='k')
        ax.vlines(np.median(btable['bgb'][btable['btype']==bt]),ax.get_ylim()[0],ax.get_ylim()[1],lw=2,linestyle='--',color=clr)

    for ax in axarr[1,:]:
        ax.set_xlabel(r'$B_{gB}$',fontsize=28)

    for ax in axarr[:,0]:
        ax.set_ylabel('Count',fontsize=28)


    fig.tight_layout()
    if savefig:
        fig.savefig('%s/paper/figures/bgb_hist.pdf' % blazardir)
    else:
        plt.show()


    return None

def plot_bgb_redshift(btable,savefig=False):

    fig = plt.figure(2,figsize=(9,7))
    fig.clf()
    ax = fig.add_subplot(111)

    ax.errorbar(btable['z'],btable['bgb'],btable['bgb_err'],color='grey',alpha=0.2,fmt='.')

    # Bin by type and overplot

    bllac = (btable['btype'] == 'bllac')
    fsrq = (btable['btype'] == 'fsrq')

    binwidth = 0.05
    lowz = 0.043
    highz = 0.75

    zbins = np.arange(lowz,highz,binwidth)
    z1 = zbins[:-1]
    z2 = zbins[1:]

    z_fsrq,z_bllac = [],[]
    z_fsrq_err,z_bllac_err = [],[]
    for zl,zh in zip(z1,z2):
        zindex = (btable['z'] >= zl) & (btable['z'] < zh)
        z_fsrq.append(np.mean(btable[zindex & fsrq]['bgb']))
        z_bllac.append(np.mean(btable[zindex & bllac]['bgb']))
        z_fsrq_err.append(np.std(btable[zindex & fsrq]['bgb']))
        z_bllac_err.append(np.std(btable[zindex & bllac]['bgb']))

        
    l1 = ax.errorbar(z1+binwidth,z_fsrq,z_fsrq_err,color='red',label='FSRQ',markersize=10,fmt='o')
    l2 = ax.errorbar(z1+binwidth,z_bllac,z_bllac_err,markerfacecolor='white',markeredgecolor='none',label='BL Lac',markersize=1,fmt='.',ecolor='blue')
    l3 = ax.scatter(z1+binwidth,z_bllac,facecolors='none',linewidth=10,edgecolors='blue',label='BL Lac',s=15)

    ax.set_xlabel('Redshift',fontsize=20)
    ax.set_ylabel(r'$B_{gB}$ [Mpc$^{-1.77}$]',fontsize=30)

    ax.set_ylim(-1000,1000)

    # Now add the legend with some customizations.
    legend = ax.legend((l1,l3),('FSRQ','BL Lac'),loc='upper right', shadow=True,scatterpoints=2)
    
    # The frame is matplotlib.patches.Rectangle instance surrounding the legend.
    #frame = legend.get_frame()
    #frame.set_facecolor('0.90')
    
    # Set the fontsize
    for label in legend.get_texts():
        label.set_fontsize(20)
    
    fig.tight_layout()
    if savefig:
        fig.savefig('%s/paper/figures/bgb_redshift.pdf' % blazardir)
    else:
        plt.show()


    return None

def check_circularity(rows):

    # Figure out whether the distribution of points around a center indicates that it's near the edge of a run

    from astropy import units as u
    from astropy.coordinates import SkyCoord,Angle

    ra_cen = rows[0]['ra']
    dec_cen = rows[0]['dec']

    position_angles = []
    c1 = SkyCoord(ra_cen,dec_cen,unit=u.degree,frame='fk5')
    c2 = SkyCoord([rows['neighbor_ra']],[rows['neighbor_dec']],unit=u.degree,frame='fk5')
    position_angles = c1.position_angle(c2).degree
    n,bins = np.histogram(position_angles[0,:],bins=np.arange(37)*10.,normed=True)

    if max(n) > 0.01:
        is_circular = False
    else:
        is_circular = True

    return is_circular

def plot_circularity(data):

    blazars = set(data['bname'])

    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(121)
    ax_wrapped = fig.add_subplot(122)

    fig2 = plt.figure(figsize=(10,10))

    for i,b in enumerate(blazars):
    
        name_match = (data['bname'] == b)
        rows = data[name_match]
        ra_cen = rows['ra'][0]
        dec_cen = rows['dec'][0]
        position_angles = check_circularity(ra_cen,dec_cen,rows,wrap=False)
        position_angles_wrapped = check_circularity(ra_cen,dec_cen,rows,wrap=True)
    
        bins = np.arange(37)*10.
        if name_match.sum() > 25:
            n,bins,edges = ax.hist(position_angles[0,:],bins=bins,histtype='step',color='k',alpha=0.2,normed=True)
            n_wrapped,bins_wrapped,edges_wrapped = ax_wrapped.hist(position_angles_wrapped[0,:],bins=bins,histtype='step',color='b',alpha=0.2,normed=True)
            if max(n) > 0.0075:
                pcolor = 'r'
            elif max(n_wrapped) > 0.010:
                pcolor = 'b'
            else:
                pcolor = 'k'

        axp = fig2.add_subplot(9,9,i+1)
        axp.scatter(rows['neighbor_ra'],rows['neighbor_dec'],color=pcolor)

    fig.tight_layout()
    plt.show()
    
    return None

def bzcat5(nmin=10):

    from astropy.io import ascii
    #data = ascii.read("../bzcat5/bzcat_sdss.csv",format='csv',guess=False)
    data = ascii.read("../newdata/all_blazars_sdss.csv",format='csv',guess=False)
    
    blazars = set(data['bname'])

    bad = 0
    good = 0
    noncirc = 0

    col_bname = []
    col_btype = []
    col_ra = []
    col_dec = []
    col_catalogue = []
    col_nt = []
    col_nb = []
    col_z = []
    col_field_size = []
    col_counting_mag = []
    col_bgb = []
    col_bgb_err = []

    from datetime import datetime as dt

    print dt.today().strftime("%H:%M:%S.%f")
    for i,b in enumerate(blazars):
    
        name_match = (data['bname'] == b)
        ngals = np.sum(name_match)

        rows = data[name_match]
    
        nt,nb,field_size,z,counting_mag = bstats(rows)
        b_val,b_err = get_bgb(nt, nb, field_size, z, counting_mag, lf = 'dr6ages')
            
        if nt+nb >= nmin:
            if check_circularity(rows):
                good += 1
                #print '{2:20} -- nt={3:4d}, nb={4:4d}, field={5:.1f}, z={6:.2f}, mag={7:.1f}, B_gB = ({0:.2f} +- {1:.2f}) Mpc^1.77'.format(b_val.value,b_err.value,b,nt,nb,field_size,z,counting_mag)
                col_bname.append(b)
                col_btype.append((data[name_match][0]['btype']).strip())
                col_ra.append((data[name_match][0]['ra']))
                col_dec.append((data[name_match][0]['dec']))
                col_catalogue.append((data[name_match][0]['catalog']).strip())
                col_nt.append(nt)
                col_nb.append(nb)
                col_z.append(z)
                col_field_size.append(field_size.value)
                col_counting_mag.append(counting_mag)
                col_bgb.append(b_val.value)
                col_bgb_err.append(b_err.value)

            else:
                noncirc += 1
                #print '{0} is non-circular'.format(b)
        else:
            bad += 1
            #print '{0} had only {1:d} galaxies within r_proj = 500 kpc'.format(b,nt)
            #pass

        if not (i % 50):
            print "Completed {0}/{1} at {2}".format(i,len(blazars),dt.today().strftime("%H:%M:%S.%f"))
    
    print "\n{0} with BgB, {1} with < {3}, {2} non-circular".format(len(blazars)-bad-noncirc,bad,noncirc,nmin)
    print dt.today().strftime("%H:%M:%S.%f")

    btable = Table([col_bname,col_btype,col_ra,col_dec,col_catalogue,col_nt,col_nb,col_z,col_field_size,col_counting_mag,col_bgb,col_bgb_err], names=('bname','btype','ra','dec','catalogue','nt','nb','z','field_size','counting_mag','bgb','bgb_err'), meta={'name': 'All blazars in SDSS'})

    #btable.write('../bzcat5/bzcat5_bgb.fits',format='fits')

    return btable


