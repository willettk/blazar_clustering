from astropy.io import fits
from astropy import units as u

import numpy as np

from matplotlib import pyplot as plt

import bgb

def plot_blazar_sequence():

    blazar_dir = '/Users/willettk/Astronomy/Research/blazars'
    all_blazars = fits.getdata('%s/fits/meyer_bz_allseds.fits' % blazar_dir,ext=1)
    
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
    zarr = np.array(ztemp)
    
    # Limit blazars to redshift range and matches to Meyer+11
    redshift_limits = (zarr > lowz) & (zarr < highz)
    blazar_type_limits = (all_blazars['sed_code'] == 'Uex') | (all_blazars['sed_code'] == 'Tex')
    counting_limits = (all_blazars['n500_cmag'] >= neighbors) & (all_blazars['bg'] > 0)

    ind = redshift_limits & blazar_type_limits
    b = all_blazars[ind]
    
    # Sort blazars by overall type
    bllac = (b['btype'] == 'BLLac') | (b['btype'] == 'Plotkin_blazar') | (b['btype'] == 'HBL') | (b['btype'] == 'lBLLac')
    fsrq = (b['btype'] == 'FSRQ')
    uncertain = (b['btype'] == 'BLLac_candidate') | (b['btype'] == 'blazar_uncertain') | (b['btype'] == 'Blazar')
    
    # Run bgb.py on only the Meyer+11 sample
    bmeyer,bmeyer_err = [],[]
    for blazar,z in zip(b,zarr):
        bgb_val,bgb_err= bgb.bgb(blazar['n500_cmag'], blazar['bg'], blazar['fieldsize']*u.arcsec, z, blazar['cmag'], lf = 'dr6ages')
        bmeyer.append(bgb_val.value)
        bmeyer_err.append(bgb_err.value)

    # Axis limits for the plot
    bmin = -500.
    bmax = 1500.
    
    leftx = 0.12
    rightx = 0.20
    yheight = 0.80
    xwidth = 1. - leftx - rightx
    xsize = 8
    ysize = xsize * xwidth / yheight
    fig = plt.figure(1,(xsize,ysize))
    fig.clf()
    ax = fig.add_subplot(111,position=[0.12,0.10,xwidth,yheight])

    # Set colormap for the symbols

    cmap = plt.cm.get_cmap('jet')
    
    # Plot the blazars on the nu_peak vs. peak luminosity plot. Color = redshift, size = B_gB
    
    nupeak = np.array(b['nupeak'],dtype=float)
    lpeak = np.array(b['lpeak'],dtype=float)

    sizearr = (np.array(bmeyer) - np.min(bmeyer)) / 100.

    sc = ax.scatter(nupeak[bllac],lpeak[bllac],marker='o',label='BL Lac',c = zarr[bllac],s=sizearr[bllac],cmap=cmap,vmin=0,vmax=0.75)
    ax.scatter(nupeak[fsrq],lpeak[fsrq],marker='s',label='FSRQ',c = zarr[fsrq],s=sizearr[fsrq],cmap=cmap,vmin=0,vmax=0.75)
    ax.scatter(nupeak[uncertain],lpeak[uncertain],marker='+',label='Uncertain',c = zarr[uncertain],s=sizearr[uncertain],cmap=cmap,vmin=0,vmax=0.75)

    ax.set_xlabel(r'$\log(\nu_{\rm peak})$ [Hz]',fontsize=16)
    ax.set_ylabel(r'$\log(\nu {\rm L}_\nu$) [erg s$^{-1}$]',fontsize=16)

    ax.set_xlim(12,17)
    ax.set_ylim(43.5,47)

    ax.legend(scatterpoints=2)

    position=fig.add_axes([0.82,0.45,0.05,0.40])
    cb = plt.colorbar(sc,cax = position, orientation='vertical')
    cb.set_label('blazar redshift',fontsize=16)

    
    fig.savefig('%s/paper/figures/bgb_blazarsequence_allseds.pdf' % blazar_dir)
    #plt.show()

    return zarr

def idl_junk():

    '''
    #fig.savefig('%s/paper/figures/bgb_blazarsequence_allseds.pdf' % blazar_dir)
    
    cgloadct, 13
    minsize = 0.5
    sizescale = 5.
    
    # Overplot tracks from Meyer et al. (2011)
    
    readcol, '/Applications/Dexter/meyer_sequence.gif.tracka', xa, ya, format='f,f', /skipline, /silent
    readcol, '/Applications/Dexter/meyer_sequence.gif.trackb', xb, yb, format='f,f', /skipline, /silent
    
    slope = (ya[1] - ya[0]) / (xa[1] - xa[0])
    intercept = ya[0] - slope * xa[0]
    
    ;cgplots, xa, ya, color='black', linestyle=2
    ;cgplots, xb, yb, color='black', linestyle=2
    
    cgplot, b.nupeak, b.lpeak, $
        background='white', $
        position=[0.08,0.15,0.30,0.95], $
        charsize=cs, $
        charthick = ct, $
        thick=th, $
        /nodata, $
        yr=[43.95,47], /ystyle, $
        xr=[12,17], /xstyle, $
        title='BL Lacs', $
        xtitle='log ('+greek('nu')+'!Ipeak!N) [Hz]', $
        ytitle='log ('+greek('nu')+'L!I'+greek('nu')+'!N) [erg s!E-1!N]'
    
    for j=0,n_elements(bllac) - 1 do begin
        cgplot, b[bllac[j]].nupeak, b[bllac[j]].lpeak, $
            /over, $
            symsize=(bmeyer[bllac[j]] - bmin) / (bmax-bmin) * sizescale + minsize, $
            color=fix((z[bllac[j]] - min(z)) / (max(z)-min(z)) * 255.), $
            psym=16
    endfor
    
    cgarrow, (!Y.crange[1] - intercept) / slope, !y.crange[1],(!Y.crange[0] - intercept) / slope, !y.crange[0], color='black', linestyle=2, thick=th,/data,/solid,hsize=hsize
    cgplots, [15.80,16.47,17.00], [44.24,44.50,44.82], color='black', linestyle=2, thick=th
    cgarrow, 15.80,44.24, 15.05,44.00,color='black', linestyle=2, thick=th, /data, hsize=hsize,/solid
    
    cgplot, b.nupeak, b.lpeak, $
        background='white', $
        position=[0.37,0.15,0.58,0.95], $
        charsize=cs, $
        charthick = ct, $
        thick=th, $
        /nodata, $
        yr=[43.95,47], /ystyle, $
        xr=[12,17], /xstyle, $
        title='FSRQs', $
        xtitle='log ('+greek('nu')+'!Ipeak!N) [Hz]'
    
    for j=0,n_elements(fsrq) - 1 do begin
        cgplot, b[fsrq[j]].nupeak, b[fsrq[j]].lpeak, $
            /over, $
            symsize=(bmeyer[fsrq[j]] - bmin) / (bmax-bmin) * sizescale + minsize, $
            color=fix((z[fsrq[j]] - min(z)) / (max(z)-min(z)) * 255.), $
            psym=15
    endfor
    
    cgarrow, (!Y.crange[1] - intercept) / slope, !y.crange[1],(!Y.crange[0] - intercept) / slope, !y.crange[0], color='black', linestyle=2, thick=th,/data,/solid,hsize=hsize
    cgplots, [15.80,16.47,17.00], [44.24,44.50,44.82], color='black', linestyle=2, thick=th
    cgarrow, 15.80,44.24, 15.05,44.00,color='black', linestyle=2, thick=th, /data, hsize=hsize,/solid
    
    cgplot, b.nupeak, b.lpeak, $
        background='white', $
        position=[0.62,0.15,0.80,0.95], $
        charsize=cs, $
        charthick = ct, $
        thick=th, $
        /nodata, $
        yr=[43.95,47], /ystyle, $
        xr=[12,17], /xstyle, $
        title='Uncertain blazars', $
        xtitle='log ('+greek('nu')+'!Ipeak!N) [Hz]'
    
    for j=0,n_elements(uncertain) - 1 do begin
        cgplot, b[uncertain[j]].nupeak, b[uncertain[j]].lpeak, $
            /over, $
            symsize=(bmeyer[uncertain[j]] - bmin) / (bmax-bmin) * sizescale + minsize, $
            color=fix((z[uncertain[j]] - min(z)) / (max(z)-min(z)) * 255.), $
            psym=14
    
    
    
    al_legend, /top,/left, psym=[16,15,34], ['BL Lac','FSRQ','Uncertain'], charsize=labelsize*1.4, symsize=2, outline_color='black', textcolor='black', colors='black'
    cgcolorbar, position=[0.9,0.5,0.95,0.95], range=[min(z),max(z)], /vert, title='blazar redshift', color='black'
    
    s1 = (-500 - bmin) / (bmax-bmin) * sizescale + minsize
    s2 = (0 - bmin) / (bmax-bmin) * sizescale + minsize
    s3 = (500 - bmin) / (bmax-bmin) * sizescale + minsize
    s4 = (1000 - bmin) / (bmax-bmin) * sizescale + minsize
    al_legend, position=[0.84,0.35], psym=16, ['-500','0','500','1000'], charsize=cs, symsize=[s1,s2,s3,s4],/normal,spacing = 3, outline_color='black', textcolor='black', colors='black'
    cgtext, /normal, 0.895, 0.37, 'B!IgB!N',charsize=cs
    
    # Kluges to make sure lines don't run off the edge of my truncated plot
    
    cgarrow, (!Y.crange[1] - intercept) / slope, !y.crange[1],(!Y.crange[0] - intercept) / slope, !y.crange[0], color='black', linestyle=2, thick=th,/data,/solid,hsize=hsize
    cgplots, [15.80,16.47,17.00], [44.24,44.50,44.82], color='black', linestyle=2, thick=th
    cgarrow, 15.80,44.24, 15.05,44.00,color='black', linestyle=2, thick=th, /data, hsize=hsize,/solid
    
    #cgtext, 13.7, 46.7, /data, 'Single-component jet', charsize=labelsize
    #cgtext, 15.85, 44.2, /data, 'Decelerating jets', charsize=labelsize
    
    '''
    
    
    '''
    # Split between low and high-peaked environments?
    
    lowpeak = where(b.nupeak lt 14.5)
    highpeak = where(b.nupeak gt 14.5)
    highpower = where(b.lpeak gt 45.0)
    lowpower = where(b.lpeak le 45.0)
    
    print,''
    print,'B_gB lowpeak: ',mean(bmeyer[lowpeak]),stddev(bmeyer[lowpeak]), median(bmeyer[lowpeak])
    print,'B_gB highpeak: ',mean(bmeyer[highpeak]),stddev(bmeyer[highpeak]), median(bmeyer[highpeak])
    print,''
    print,'redshift lowpeak: ',mean(z[lowpeak]),stddev(z[lowpeak])
    print,'redshift highpeak: ',mean(z[highpeak]),stddev(z[highpeak])
    print,''
    kstwo, bmeyer[lowpeak],bmeyer[highpeak],d_bgb,prob_bgb
    print,'K-S BgB, split by peak:', prob_bgb, probgauss(prob_bgb)        # Probability of B_gB being drawn from the same distribution is 0.3% (3 sigma)
    kstwo, z[lowpeak],z[highpeak],d_z,prob_z            # Probability of redshift being drawn from the same distribution is 2.4% (2.3 sigma)
    print,'K-S redshift, split by peak:', prob_z, probgauss(prob_z)
    kstwo, b[lowpeak].lpeak,b[highpeak].lpeak,d_lpeak,prob_lpeak            # Probability of redshift being drawn from the same distribution is 2.4% (2.3 sigma)
    print,'K-S L_peak, split by peak:', prob_lpeak, probgauss(prob_lpeak)
    kstwo, bmeyer[lowpower],bmeyer[highpower],d_lpower,prob_lpower            # Probability of redshift being drawn from the same distribution is 2.4% (2.3 sigma)
    print,'K-S L_peak, split by power:', prob_lpower, probgauss(prob_lpower)
    
    # Blazar type in each
    #print,'High peaked'
    print,''
    
    #stop
    #
    ## Correlations
    #
    #mbh_flt = float(b.mass_bh)
    #mbhind = where(mbh_flt ne 0.)
    #lext_flt = float(b.lext)
    #lextind = where(lext_flt ne 0.)
    #
    #!p.multi=[0,2,2]
    #cgplot,bmeyer,b.lpeak,psym=16,xr=[-1000,1500],yr=[42,47],title='L!Ipeak!N', charsize=2
    #cgplot,bmeyer,b.nupeak,psym=16,xr=[-1000,1500],yr=[12,17],title=greek('nu')+'!Ipeak!N', charsize=2
    #cgplot,bmeyer,b.mass_bh,psym=16,xr=[-1000,1500],yr=[6,10],title='M!IBH!N', charsize=2
    #cgplot,bmeyer,b.lext,psym=16,xr=[-1000,1500],yr=[37,44],title='L!Iext!N', charsize=2
    #
    #print,'Bgb vs lpeak:   ',string(correlate(bmeyer,float(b.lpeak)),format='(f7.4)')
    #print,'Bgb vs nupeak:  ',string(correlate(bmeyer,float(b.nupeak)),format='(f7.4)')
    #print,'Bgb vs BH mass: ',string(correlate(bmeyer[mbhind],mbh_flt[mbhind]),format='(f7.4)')
    #print,'Bgb vs L_ext:   ',string(correlate(bmeyer[lextind],lext_flt[lextind]),format='(f7.4)')
    '''


if __name__ == "__main__":
    plt.ion()
    plot_blazar_sequence()
