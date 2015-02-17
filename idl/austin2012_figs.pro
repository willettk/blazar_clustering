;+
; NAME:
;       
;	
;
; PURPOSE:
;
;	
;
; INPUTS:
;
;	
;
; OUTPUTS:
;
;	
;
; KEYWORDS:
;
;	
;
; EXAMPLE:
;
;	IDL> 
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                
;-

pro austin2012_figs, stop=stop, noplot = noplot, ps=ps

	blazardir = '~/Astronomy/Research/blazars/'
	figdir = '~/Astronomy/meetings/austin2012/'

	restore,blazardir+'ab_structure.sav'

	; total N_gal, background N_gal, radius of field in arcsec, redshift of blazar

	ind = where(bz.z gt 0.162 and bz.z lt 0.75)
	bz=bz[ind]

	bllac = where(bz.type eq 'BLLac')
	fsrq = where(bz.type eq 'FSRQ')
	uncertain = where(bz.btype eq 'BLLac_candidate' or bz.btype eq 'Blazar')

	annulus_area, bz.cmag_app, ng_per_arcmin2, /quiet, /noplot

	fieldradius = 500	; kpc
	bzz = bz.z
	fieldsize = zang(float(fieldradius),bzz,/silent,/wmap7)
	nt = float(bz.n500_cmag)
	nb = !pi * (fieldsize/60.)^2 * ng_per_arcmin2
	countingmag = bz.cmag

	bgb, nt, nb, fieldsize, bzz, countingmag, bgb_blazar, bgb_blazar_err

		kstwo, double(bgb_blazar[fsrq]), double(bgb_blazar[bllac]), d, prob
		print,'BL Lac/FSRQ KS result:',d,prob,sqrt(2d) * inverf(1d - prob)

	print,'Blazars: ', mean(bgb_blazar),stddev(bgb_blazar), mean(bzz), n_elements(bgb_blazar)
	print,'BL Lac: ', mean(bgb_blazar[bllac]),stddev(bgb_blazar[bllac]), mean(bzz[bllac]), n_elements(bllac)
	print,'FSRQ: ', mean(bgb_blazar[fsrq]),stddev(bgb_blazar[fsrq]), mean(bzz[fsrq]), n_elements(fsrq)
	print,'Uncertain: ', mean(bgb_blazar[uncertain]),stddev(bgb_blazar[uncertain]), mean(bzz[uncertain]), n_elements(uncertain)

	; FRI and FRII galaxies

	restore,blazardir+'tcrr_structure.sav'
	tind = where(tcrr.z gt 0.162 and tcrr.z lt 0.75 and tcrr.frclass ge 1)
	tcrr=tcrr[tind]

	fr1 = where(tcrr.frclass eq 1)
	fr2 = where(tcrr.frclass eq 2)

	annulus_area, tcrr.cmag_app, ng_per_arcmin2, /quiet, /noplot

	fieldradius = 500	; kpc
	z = tcrr.z
	fieldsize = zang(float(fieldradius),z,/silent,/wmap7)
	nt = float(tcrr.n500_cmag)
	nb = !pi * (fieldsize/60.)^2 * ng_per_arcmin2
	countingmag = tcrr.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bgb_tcrr, bgb_tcrr_err

	print,'Radio galaxies: ', mean(bgb_tcrr),stddev(bgb_tcrr), mean(z), n_elements(bgb_tcrr)
	print,'FR I: ', mean(bgb_tcrr[fr1]),stddev(bgb_tcrr[fr1]), mean(z[fr1]), n_elements(fr1)
	print,'FR II: ', mean(bgb_tcrr[fr2]),stddev(bgb_tcrr[fr2]), mean(z[fr2]), n_elements(fr2)

	;	kstwo, bgb_tcrr[fr1], bgb_tcrr[fr2], d_tcrr, prob_tcrr
	;	print,'FRI/FRII KS result:',d_tcrr,prob_tcrr,sqrt(2d) * inverf(1d - prob_tcrr)

	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=figdir+'bgb_blazars.eps', /quiet, /color, /encap, ysize=10, xsize=12
		cs = 1.5
		labelsize=1.5
		th = 5
	endif else begin
		cs = 2
		labelsize = 1.5
	endelse
	
	; Plot the distribution of B_gB

	bgbbin = 200
	cghistoplot, bgb_blazar, $
		thick=th, ythick=th, xthick=th, $
		datacolor='black', $
		binsize=bgbbin, $
		;xr=[-500,2500], $
		xtitle='B!IgB!N', $
		ytitle='N (all blazars)', $
		charsize=cs

	cgplots, [mean(bgb_blazar),mean(bgb_blazar)], !y.crange, linestyle=2, thick=2, color='Black'

	cghistoplot, bgb_blazar[bllac], $
		thick=th, ythick=th, xthick=th, $
		binsize=bgbbin, $
		yr=[0,300], $
		xtitle='B!IgB!N', $
		ytitle='N (by type)', $
		datacolor="Blue",$ 
		charsize=cs

	cghistoplot,bgb_blazar[fsrq], /oplot, $
		thick=th, $
		binsize=bgbbin, $
		datacolor="Red"

	cgplots, [mean(bgb_blazar[fsrq]),mean(bgb_blazar[fsrq])], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, [mean(bgb_blazar[bllac]),mean(bgb_blazar[bllac])], !y.crange, linestyle=2, color="Blue", thick=2

	al_legend, /top, /right, $
		thick=th, bthick=th, $
		['FSRQ','BL Lac'], $
		color=['Red','Blue'], $
		;linestyle=0, $
		psym=28, $
		charsize=labelsize

	th = 5 

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	cghistoplot, bgb_tcrr[fr2], $
		thick=th, ythick=th, xthick=th, $
		binsize=100, $
		xtitle='B!IgB!N', $
		ytitle='N (FRII radio galaxies)', $
		datacolor="Dark Green",$ 
		charsize=cs

	cgplots, [median(bgb_tcrr[fr2]),median(bgb_tcrr[fr2])], !y.crange, linestyle=2, color="Dark Green", thick=2

	th = 5

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	; Plot binned Bgb as fn. of redshift

	binsize = 0.05
	bins = fillarr(binsize,0.162,0.8)
	nb = n_elements(bins)
	fsrq_binned = fltarr(nb)
	fsrq_binned_err = fltarr(nb)
	bllac_binned = fltarr(nb)
	bllac_binned_err = fltarr(nb)
	fsrq_bgb = bgb_blazar[fsrq]
	bllac_bgb = bgb_blazar[bllac]
	for i=0,nb-2 do begin
		junk_fsrq = where(bzz[fsrq] gt bins[i] and bzz[fsrq] lt bins[i+1])
		fsrq_binned[i] = mean(fsrq_bgb[junk_fsrq])
		fsrq_binned_err[i] = stddev(fsrq_bgb[junk_fsrq])
		junk_bllac = where(bzz[bllac] gt bins[i] and bzz[bllac] lt bins[i+1])
		bllac_binned[i] = mean(bllac_bgb[junk_bllac])
		bllac_binned_err[i] = stddev(bllac_bgb[junk_bllac])
	endfor

	zarr = bins[0:nb-2]+binsize/2.

	;ploterror, bzz[fsrq], bgb_blazar[fsrq], bgb_blazar_err[fsrq], $
	ploterror, zarr, fsrq_binned, fsrq_binned_err, $
		thick=th, ythick=th, xthick=th, $
		color="Red", $
		errcolor="red", $
		charsize=cs, $
		xrange=[0,0.8], /xstyle, $
		/nohat, $
;		yr=[-5d3,2d4], $
		xtitle='Blazar redshift', $
		ytitle='B!IgB!N', $
		psym = 16

	;oploterror, bzz[bllac], bgb_blazar[bllac], bgb_blazar_err[bllac], $
	oploterror, zarr, bllac_binned, bllac_binned_err, $
		thick=th, $
		color="Blue", $
		/nohat, $
		psym = 9

	al_legend, /top, /right, $
		thick=th, bthick=th, $
		['FSRQ','BL Lac'], $
		color=['Red','Blue'], $
		psym=[16,9], $
		charsize=labelsize

	th = 5

	if keyword_set(ps) then ps_end


	endif	; noplot

	!p.multi=[0,1,1]	

	if keyword_set(ps) then begin
		ps_start, filename=figdir+'bgb_multiwave.eps', /quiet, /color, /encap, ysize=10, xsize=14
		cs = 2.5
		labelsize=1.3
		th = 5
	endif else begin
		cs = 2
		labelsize = 1.5
	endelse
	
	!p.multi=[0,3,2]

	mpc2cm = 3.086d24
	mpc2m = 3.086d22
	mjy2cgs = 1d-26
	mjy2si = 1d-29
	xray2cgs = 1d-12

	ind = where(bz.flux_radio gt 0.)
	ind_bllac = where(bz.flux_radio gt 0. and bz.btype eq 'BLLac')
	ind_fsrq = where(bz.flux_radio gt 0. and bz.btype eq 'FSRQ')
	cgplot, bz[ind].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind].z,/silent,/wmap7)*mpc2m)^2, bgb_blazar[ind], psym=9, xtitle='1.4 GHz L!I'+greek('nu')+'!N [W/Hz]', ytitle='B!IgB!N',/xlog, thick=th, xthick=th, ythick=th, charsize=cs, /nodata
	cgplot, /over, bz[ind_bllac].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2m)^2, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue'
	cgplot, /over, bz[ind_fsrq].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2m)^2, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red'
	cgtext, 2d23, -150, /data, color='blue', greek('rho')+' = '+string(correlate(alog10(bz[ind_bllac].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2m)^2),bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
	cgtext, 2d23, -300, /data, color='red', greek('rho')+' = '+string(correlate(alog10(bz[ind_fsrq].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2m)^2),bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

	al_legend, /top, /right, $
		thick=th, bthick=th, $
		['FSRQ','BL Lac'], $
		color=['Red','Blue'], $
		psym=[7,16], $
		charsize=labelsize

	th = 5 

	ind = where(bz.mag_r gt 0.)
	ind_bllac = where(bz.mag_r gt 0. and bz.btype eq 'BLLac')
	ind_fsrq = where(bz.mag_r gt 0. and bz.btype eq 'FSRQ')
	cgplot, bz[ind].mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind], psym=9, xtitle='M!IR!N', ytitle='B!IgB!N', thick=th, xthick=th, ythick=th, charsize=cs, /nodata
	cgplot, /over, bz[ind_bllac].mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind_bllac], psym=16, thick=th, color='blue'
	cgplot, /over, bz[ind_fsrq].mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind_fsrq], psym=7, thick=th, color='red'
	cgtext, -29.5, -150, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].mag_r - 5*alog10(lumdist(bz[ind_bllac].z,/silent,/wmap7)*1e6) + 5.,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
	cgtext, -29.5, -300, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].mag_r - 5*alog10(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*1e6) + 5.,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

	ind = where(bz.flux_xray gt 0.)
	ind_bllac = where(bz.flux_xray gt 0. and bz.btype eq 'BLLac')
	ind_fsrq = where(bz.flux_xray gt 0. and bz.btype eq 'FSRQ')
	cgplot, bz[ind].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind], psym=9, xtitle='(0.1-2.4) keV '+greek('nu')+'L!I'+greek('nu')+'!N [erg/s]', ytitle='B!IgB!N',/xlog, thick=th, xthick=th, ythick=th, charsize=cs, /nodata
	cgplot, /over, bz[ind_bllac].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue'
	cgplot, /over, bz[ind_fsrq].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red'
	cgtext, 2d42, -150, /data, color='blue', greek('rho')+' = '+string(correlate(alog10(bz[ind_bllac].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2cm)^2),bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
	cgtext, 2d42, -300, /data, color='red', greek('rho')+' = '+string(correlate(alog10(bz[ind_fsrq].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2cm)^2),bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

	ind = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0.)
	ind_bllac = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0. and bz.type eq 'BLLac')
	ind_fsrq = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0. and bz.btype eq 'FSRQ')
	cgplot, bz[ind].spindex_ro, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (radio-optical)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata
	cgplot, bz[ind_bllac].spindex_ro, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
	cgplot, bz[ind_fsrq].spindex_ro, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
	cgtext, 0.05, -150, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_ro,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
	cgtext, 0.05, -300, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_ro,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

	ind = where(bz.spindex_ox gt -9.9 and bz.spindex_ox ne 0.)
	ind_bllac = where(bz.spindex_ox gt -9.9 and bz.type eq 'BLLac')
	ind_fsrq = where(bz.spindex_ox gt -9.9 and bz.spindex_ox ne 0. and bz.btype eq 'FSRQ')
	cgplot, bz[ind].spindex_ox, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (optical-X-ray)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata
	cgplot, bz[ind_bllac].spindex_ox, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
	cgplot, bz[ind_fsrq].spindex_ox, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
	cgtext, 0.60, -150, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_rx,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
	cgtext, 0.60, -300, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_rx,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

	ind = where(bz.spindex_rx gt -9.9 and bz.spindex_rx ne 0.)
	ind_bllac = where(bz.spindex_rx gt -9.9 and bz.type eq 'BLLac')
	ind_fsrq = where(bz.spindex_rx gt -9.9 and bz.spindex_rx ne 0. and bz.btype eq 'FSRQ')
	cgplot, bz[ind].spindex_rx, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (radio-X-ray)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata
	cgplot, bz[ind_bllac].spindex_rx, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
	cgplot, bz[ind_fsrq].spindex_rx, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
	cgtext, 0.45, -150, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_rx,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
	cgtext, 0.45, -300, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_rx,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

	if keyword_set(ps) then ps_end

	; Plot number of neighbors for each blazar within 500 kpc

	!p.multi=[0,1,1]	

	if keyword_set(ps) then begin
		ps_start, filename=figdir+'blazar_neighbors.eps', /quiet, /color, /encap, ysize=5, xsize=7
		cs = 2.0
		labelsize=1.3
		th = 5
	endif else begin
		cs = 2
		labelsize = 1.5
	endelse

	cghistoplot, bz.n500_cmag, $
		datacolor='black', $
		xtitle='N!Igal!N (R < 500 kpc)', $
		ytitle='N!Iblazar!N', $
		charsize=cs, $
		thick=th, ythick=th, xthick=th
	

	if keyword_set(ps) then ps_end

	if keyword_set(stop) then stop

end
