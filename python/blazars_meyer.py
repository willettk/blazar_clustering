from __future__ import division

from astropy.io import fits
from astropy import units as u
from astroML.plotting import hist as histML
from matplotlib import pyplot as plt
from matplotlib import cm

import numpy as np
import bgb

fits_path = "../fits"
paper_path = "../paper"

def get_data():

    # Plot the peak frequency vs. peak luminosity for the blazar sequence,
    # distinguishing points by both redshift and clustering strength (B_gB)

    # Based on Figure 4 in Meyer et al. (2011)

    all_blazars = fits.getdata('{0}/meyer_bz_allseds.fits'.format(fits_path),1)
    
    lowz = 0.043    # Redshift at which angular size of 10 arcmin = 500 kpc
    highz = 0.75    # Redshift at which an M* galaxy reaches the SDSS completeness limit
    neighbors = 10
    
    ztemp = []
    for blazar in all_blazars:
        try:
            z = np.float(blazar['used_redshift'])
        except ValueError:
            z = blazar['z']
        ztemp.append(z)
    _zarr = np.array(ztemp)
    
    # Include only blazars with redshifts in SDSS completeness limit
    redshift_limits = (_zarr > lowz) & (_zarr < highz)
    # Include only blazars in the "trusted extended" or "unknown extended" categories of Meyer+11
    blazar_type_limits = (all_blazars['sed_code'] == 'Uex') | (all_blazars['sed_code'] == 'Tex')
    # Include only blazars that had a sufficient number of galaxies in SDSS to estimate environment
    counting_limits = (all_blazars['n500_cmag'] >= neighbors) & (all_blazars['bg'] > 0)

    b = all_blazars[redshift_limits & blazar_type_limits]
    zarr = _zarr[redshift_limits & blazar_type_limits]

    
    # Compute spatial clustering amplitude for the sample

    bmeyer,bmeyer_err = [],[]
    for blazar,z in zip(b,zarr):
        bgb_val,bgb_err= bgb.get_bgb(blazar['n500_cmag'], blazar['bg'], blazar['fieldsize']*u.arcsec, z, blazar['cmag'], lf = 'dr6ages')
        bmeyer.append(bgb_val.value)
        bmeyer_err.append(bgb_err.value)

    return bmeyer,b,zarr

def plot_blazar_sequence(bgb_data,bdata,zarr,savefig=False):

    # Axis limits for the plot
    
    leftx = 0.12
    rightx = 0.20
    yheight = 0.80
    xwidth = 1. - leftx - rightx
    xsize = 12
    ysize = xsize * xwidth / yheight
    fig = plt.figure(figsize=(xsize,ysize))
    ax = fig.add_subplot(111,position=[0.12,0.10,xwidth,yheight])

    # Plot the blazars on the nu_peak vs. peak luminosity plot. Color = redshift, size = B_gB
    
    nupeak = np.array(bdata['nupeak'],dtype=float)
    lpeak = np.array(bdata['lpeak'],dtype=float)

    sizescale = 100. / (max(bgb_data) - min(bgb_data))
    minsize = 5.

    smin = 10
    smax = 140
    bmin = -500
    bmax = 1000

    def size_bgb(x):
        shifted = x - bmin
        stemp = (shifted * (smax - smin) / (bmax - bmin)) + smin
        return stemp

    sizearr = size_bgb(np.array(bgb_data))

    def bscatter(axis,index,marker='o',label='BL Lac',vmin=0,vmax=0.75):
        sc = axis.scatter(nupeak[index],lpeak[index], 
                        marker=marker, label=label, 
                        c = zarr[index], s=sizearr[index],  
                        cmap = cm.jet,
                        vmin=vmin, vmax=vmax)

        return sc

    # Categorize blazars by spectral type
    bllac = (bdata['btype'] == 'BLLac') | (bdata['btype'] == 'Plotkin_blazar') | (bdata['btype'] == 'HBL') | (bdata['btype'] == 'lBLLac')
    fsrq = (bdata['btype'] == 'FSRQ')
    uncertain = (bdata['btype'] == 'BLLac_candidate') | (bdata['btype'] == 'blazar_uncertain') | (bdata['btype'] == 'Blazar')
    
    sc_bllac = bscatter(ax,bllac,'o','BL Lac')
    sc_fsrq = bscatter(ax,fsrq,'s','FSRQ')
    sc_uncertain = bscatter(ax,uncertain,'+','Uncertain')

    # Add dashed lines indicating the blazar sequence from theoretical predictions of Meyer+11

    # Single-component jet
    track_a = np.loadtxt("../meyer/track_a.txt")
    x0,y0 = track_a[-1]
    dx = track_a[0][0] - track_a[1][0]
    dy = track_a[0][1] - track_a[1][1]
    ax.arrow(x0,y0,dx,dy,lw=2,fc='k',ec='k',head_width=0.1)

    # Decelerating jet
    track_b = np.loadtxt("../meyer/track_b.txt")
    seg_x = track_b[1:][:,0]
    seg_y = track_b[1:][:,1]
    ax.plot(seg_x,seg_y,lw=2,color='k')
    x0,y0 = track_b[1]
    dx = track_b[0][0] - track_b[1][0]
    dy = track_b[0][1] - track_b[1][1]
    ax.arrow(x0,y0,dx,dy,lw=2,fc='k',ec='k',head_width=0.1)

    ax.set_xlabel(r'$\log(\nu_{\rm peak})$ [Hz]',fontsize=22)
    ax.set_ylabel(r'$\log(\nu {\rm L}_\nu$) [erg s$^{-1}$]',fontsize=22)

    # Original limits on Figure 4 (Meyer+11)

    ax.set_xlim(12,18)
    ax.set_ylim(41,48)

    # More sensible limits for the range of our sample

    #ax.set_xlim(12,17)
    #ax.set_ylim(43.5,47)

    # Set up dummy axis to make an extra legend for the point sizes

    xdummy,ydummy = [0],[0]
    p1Artist=ax.scatter(xdummy,ydummy,marker='o',color='k',s=size_bgb(-500),label='-500')
    p2Artist=ax.scatter(xdummy,ydummy,marker='o',color='k',s=size_bgb(0),label='0')
    p3Artist=ax.scatter(xdummy,ydummy,marker='o',color='k',s=size_bgb(500),label='500')
    p0Artist=ax.scatter(xdummy,ydummy,marker='o',color='k',s=size_bgb(1000),label='1000')
    p4Artist=ax.scatter(xdummy,ydummy,marker='o',color='k',s=size_bgb(1500),label='1500')

    m2Artist=ax.scatter(xdummy,ydummy,marker='o',color='b',s=size_bgb(0),label='BL Lac')
    m3Artist=ax.scatter(xdummy,ydummy,marker='s',color='b',s=size_bgb(0),label='FSRQ')
    m4Artist=ax.scatter(xdummy,ydummy,marker='+',color='b',s=size_bgb(0),label='Uncertain')

    handles,labels = ax.get_legend_handles_labels()

    # Main legend
    legend1=ax.legend(handles[8:],labels[8:],scatterpoints=1)
    # Sizing legend
    ax.legend(handles[3:8],labels[3:8],scatterpoints=1,bbox_to_anchor=(1.02,0.35),loc=2)
    plt.gca().add_artist(legend1)

    # Add colorbar for the blazar redshift
    cb_axis=fig.add_axes([0.82,0.45,0.05,0.40])
    cb = plt.colorbar(sc_bllac,cax = cb_axis, orientation='vertical')
    cb.set_label('blazar redshift',fontsize=16)

    if savefig:
        plt.savefig('{0}/figures/bgb_blazarsequence_allseds.pdf'.format(paper_path))
    else:
        plt.show()

def plot_blazar_sequence_hist(bgb_data,savefig=False):

    # Plot just the distribution of B_gB for the Meyer+11 sample

    fig = plt.figure()
    ax = fig.add_subplot(111)
    histML(bgb_data, bins=25, ax=ax, histtype='step', color='b',weights=np.zeros_like(bgb_data) + 1./len(bgb_data))

    ax.set_title('Matched SDSS and (TEX+UEX) samples')
    ax.set_xlabel(r'$B_{gB}$',fontsize=24)
    ax.set_ylabel('Count',fontsize=16)
    if savefig:
        fig.savefig('{0}/figures/bgb_blazarsequence_hist.pdf'.format(paper_path))
    else:
        plt.show()


    return None

if __name__ == "__main__":
    bgb_data,bdata,zarr = get_data()
    plot_blazar_sequence(bgb_data,bdata,zarr,savefig=True)
    plot_blazar_sequence_hist(bgb_data,savefig=True)
