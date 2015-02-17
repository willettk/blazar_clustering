
;+
; NAME:
;       
;	HEIDELBERG_FIGURES
;
; PURPOSE:
;
;	Make figures for the Heidelberg gamma-ray poster
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
;       Written by K. Willett                Jul 12
;-

pro heidelberg_figures, stop=stop, noplot = noplot, ps=ps, readfiles=readfiles, multiwave=multiwave, $
	neighbors = neighbors, lowz=lowz, highz=highz, cosmo=cosmo

	blazarsavdir = '~/Astronomy/Research/blazars/sav/'
	figdir = '~/Astronomy/meetings/heidelberg2012/'

	restore,blazarsavdir+'ab_structure.sav'

	; total N_gal, background N_gal, radius of field in arcsec, redshift of blazar

	if n_elements(lowz) eq 0 then lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	if n_elements(highz) eq 0 then highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	if n_elements(neighbors) eq 0 then neighbors = 10
	ind = where(bz.z lt highz and bz.z gt lowz and bz.n500_cmag ge neighbors and bz.bg gt 0)

	bz=bz[ind]

	bllac = where(bz.btype eq 'BLLac' or bz.btype eq 'Plotkin_blazar' or bz.btype eq 'HBL' or bz.btype eq 'lBLLac')
	fsrq = where(bz.btype eq 'FSRQ')
	uncertain = where(bz.btype eq 'BLLac_candidate' or bz.btype eq 'blazar_uncertain' or bz.btype eq 'Blazar')

	z = bz.z
	fieldsize = bz.fieldsize
	nt = float(bz.n500_cmag)
	nb = float(bz.bg)
	countingmag = bz.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bgb_blazar, bgb_blazar_err, cosmo=cosmo

		kstwo, bgb_blazar[fsrq], bgb_blazar[bllac], d, prob
		print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)

	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=figdir+'heidelberg_bgb.eps', /quiet, /color, /encap, ysize=10, xsize=12
		cs = 1.7
		labelsize=1.3
		th=5
	endif else begin
		cs = 2
		th=1
		labelsize = 1.5
	endelse
	
	; Plot the distribution of B_gB

	; Plot 1

	bgbbin = 200
	cghistoplot, bgb_blazar, $
		datacolor='black', $
		xr=[-2000,2000], $
		binsize=bgbbin, $
		xtitle='B!IgB!N [Mpc!E1.77!N]', $
		ytitle='Number of blazars', $
		thick=th, $
		charsize=cs

	cgplots, [mean(bgb_blazar),mean(bgb_blazar)], !y.crange, linestyle=2, thick=2, color='Black'
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.05 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		'<B!IgB!N>='+string(mean(bgb_blazar),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar),format='(i3)'), $
		charsize=labelsize

	; Plot 2: redshift distribution of the galaxies

	cghistoplot, bz.z, xtitle='Blazar redshift', ytitle='Number of blazars', charsize=cs, datacolor='black', binsize=0.05, thick=th, yr=[0,100]

	; Plot 3

	cghistoplot, bgb_blazar[fsrq], $
		binsize=bgbbin, $
		xr=[-2000,2000], $
		yr=[0,250], $
		xtitle='B!IgB!N [Mpc!E1.77!N]', $
		ytitle='Number of blazars', $
		datacolor="Red",$ 
		thick=th, $
		charsize=cs

	cghistoplot,bgb_blazar[bllac], /oplot, $
		binsize=bgbbin, $
		thick=th, $
		datacolor="Blue"

	cgplots, [mean(bgb_blazar[fsrq]),mean(bgb_blazar[fsrq])], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, [mean(bgb_blazar[bllac]),mean(bgb_blazar[bllac])], !y.crange, linestyle=2, color="Blue", thick=2
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.05 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
		color="red", $
		/data, $
		'<B!IgB!N>='+string(mean(bgb_blazar[fsrq]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar[fsrq]),format='(i3)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.05 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		color="blue", $
		'<B!IgB!N>='+string(mean(bgb_blazar[bllac]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar[bllac]),format='(i3)'), $
		charsize=labelsize

	al_legend, /top, /right, $
		bthick=th, $
		['FSRQ','BL Lac'], $
		color=['Red','Blue'], $
		psym=28 , $
		;linestyle=0, linsize=0.2, $
		;thick=th, $
		charsize=labelsize

	; Plot binned Bgb as fn. of redshift

	binsize = 0.05
	bins = fillarr(binsize,lowz,highz+binsize)
	nbins = n_elements(bins)
	ab_binned = fltarr(nbins-1)
	ab_binned_err = fltarr(nbins-1)
	fsrq_binned = fltarr(nbins-1)
	fsrq_binned_err = fltarr(nbins-1)
	bllac_binned = fltarr(nbins-1)
	bllac_binned_err = fltarr(nbins-1)
	fsrq_bgb = bgb_blazar[fsrq]
	bllac_bgb = bgb_blazar[bllac]
	fsrq_bgb_err = bgb_blazar_err[fsrq]
	bllac_bgb_err = bgb_blazar_err[bllac]
	for i=0,nbins-2 do begin
		junk_ab = where(bz.z gt bins[i] and bz.z lt bins[i+1],ja)
		if ja gt 0 then begin
			ab_binned[i] = wmean(bgb_blazar[junk_ab],bgb_blazar_err[junk_ab],error=wmean_err_ab,/nan)
			ab_binned_err[i] = stddev(bgb_blazar[junk_ab])
			ab_binned_err[i] = wmean_err_ab
		endif else print,'No FSRQs between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_fsrq = where(bz[fsrq].z gt bins[i] and bz[fsrq].z lt bins[i+1],jf)
		if jf gt 0 then begin
			fsrq_binned[i] = wmean(fsrq_bgb[junk_fsrq],fsrq_bgb_err[junk_fsrq],error=wmean_err_fsrq,/nan)
			fsrq_binned_err[i] = stddev(fsrq_bgb[junk_fsrq])
			fsrq_binned_err[i] = wmean_err_fsrq
		endif else print,'No FSRQs between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_bllac = where(bz[bllac].z gt bins[i] and bz[bllac].z lt bins[i+1],jb)
		if jb gt 0 then begin
			bllac_binned[i] = wmean(bllac_bgb[junk_bllac],bllac_bgb_err[junk_bllac],error=wmean_err_bllac,/nan)
			bllac_binned_err[i] = stddev(bllac_bgb[junk_bllac])
			bllac_binned_err[i] = wmean_err_bllac
		endif else print,'No BL Lacs between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
	endfor

	zarr = bins[0:nbins-2]+binsize/2.

	cgplot, indgen(10), $
		/nodata, $
		thick=th, ythick=th, xthick=th, $
		charsize=cs, $
		xrange=[0,round(highz*10.)/10.], /xstyle, $
;		yrange=[min(bgb_blazar),max(bgb_blazar)], /ystyle, $
		yr=[-2000,2000], $
		xtitle='Blazar redshift', $
		ytitle='B!IgB!N'

	oploterror, z, bgb_blazar, bgb_blazar_err, $
		psym=3, $
		color='grey', $
		thick=th, $
		/nohat

	;ab = where(ab_binned ne 0.)
	;oploterror, zarr[ab], ab_binned[ab], ab_binned_err[ab], $
	;	thick=th, $
	;	color="black", $
	;	/nohat, $
	;	psym = -16

	fb = where(fsrq_binned ne 0.)
	oploterror, zarr[fb], fsrq_binned[fb], fsrq_binned_err[fb], $
		thick=th, $
		color="red", $
		/nohat, $
		psym = -16

	bb = where(bllac_binned ne 0.)
	oploterror, zarr[bb], bllac_binned[bb], bllac_binned_err[bb], $
		thick=th, $
		color="Blue", $
		/nohat, $
		psym = -9

	al_legend, /top, /right, $
		bthick=th, $
		['FSRQ','BL Lac'], $
		color=['Red','Blue'], $
		psym=[16,9], $
		charsize=labelsize

		; Fitting the raw and binned distributions with a simple linear model


		linexp = 'p[0] + x*p[1]'
		xarr = fillarr(0.01,0,1)
		bgerrind = where(finite(bgb_blazar_err))
		blerrind = where(finite(bllac_binned_err) and bllac_binned_err ne 0.)
		fqerrind = where(finite(fsrq_binned_err) and fsrq_binned_err ne 0.)
		result_all = mpfitexpr(linexp,z[bgerrind],bgb_blazar[bgerrind],bgb_blazar_err[bgerrind],perr=perr_all,/quiet)
		result_bllac = mpfitexpr(linexp,zarr[blerrind],bllac_binned[blerrind],bllac_binned_err[blerrind],perr=perr_bllac,/quiet)
		result_fsrq = mpfitexpr(linexp,zarr[fqerrind],fsrq_binned[fqerrind],fsrq_binned_err[fqerrind],perr=perr_fsrq,/quiet)

		;cgplot, /over, xarr, result_all[0] + xarr*result_all[1], linestyle=2, color='black'
		;cgplot, /over, xarr, result_bllac[0] + xarr*result_bllac[1], linestyle=2, color='blue'
		;cgplot, /over, xarr, result_fsrq[0] + xarr*result_fsrq[1], linestyle=2, color='red'


		print,'All: ',result_all, perr_all
		print,'BL Lac: ',result_bllac, perr_bllac
		print,'FSRQ: ',result_fsrq, perr_fsrq


		; Does clustering amplitude correlate with redshift?

		print,'Correlation: ',correlate(bgb_blazar,z)


	;ploterror, z[fsrq], bgb_blazar[fsrq], bgb_blazar_err[fsrq], $
	;	color="Red", $
	;	errcolor="red", $
	;	charsize=cs, $
	;	xrange=[0,0.8], /xstyle, $
	;	/nohat, $
	;	xtitle='Redshift', $
	;	ytitle='B!IgB!N', $
	;	psym = 16

	;oploterror, z[bllac], bgb_blazar[bllac], bgb_blazar_err[bllac], $
	;	color="Blue", $
	;	/nohat, $
	;	psym = 9

	if keyword_set(ps) then ps_end

	!p.multi=[0,1,1]	


	endif	; noplot keyword

	; Compute how many of the neighbors might be associated with the blazar

	; Photo-z

	;nbz = n_elements(bz)
	;pzarr = lonarr(nbz)
	;meanz = fltarr(nbz)
	;meanzerr = fltarr(nbz)
	;mylim = 0.2
	;for i=0,nbz-1 do begin
	;	zblazar = bz[i].z
	;	tempz = *(bz[i].n_redshift)
	;	tempzerr = *(bz[i].n_redshift_err)
	;	tempzerr[where(tempzerr lt 0)] = 0.
	;	meanz[i] = mean(tempz)
	;	meanzerr[i] = mean(tempzerr[where(tempzerr gt 0)])
	;	junk = where(abs(zblazar - tempz) lt (tempzerr > mylim),tcount)
	;	pzarr[i] = tcount
	;endfor

	;!p.multi=[0,2,2]	
	;cghistoplot, pzarr, xtitle='N of possible spectral neighbors per blazar',charsize=1
	;cghistoplot, bz.n500, xtitle='N of neighbors per blazar',charsize=1
	;cghistoplot, bz.n500_cmag, xtitle='N of neighbors per blazar (counting mag)',charsize=1
	;cgplot, meanz, meanzerr, psym=16, xtitle='Mean redshift', ytitle='Mean redshift err',charsize=1

	if keyword_set(multiwave) then begin

		!p.multi=[0,3,2]

		if keyword_set(ps) then begin
			ps_start, filename=figdir+'heidelberg_bgb_multiwave.eps', /quiet, /color, /encap, ysize=10, xsize=14
			cs = 2.5
			labelsize=1.3
			th = 5
		endif else begin
			cs = 2
			labelsize = 1.5
		endelse
	
		mpc2cm = 3.086d24
		mpc2m = 3.086d22
		mjy2cgs = 1d-26
		mjy2si = 1d-29
		xray2cgs = 1d-12

		ind = where(bz.flux_radio gt 0.)
		ind_bllac = where(bz.flux_radio gt 0. and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.flux_radio gt 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind].z,/silent,/wmap7)*mpc2m)^2, bgb_blazar[ind], psym=9, xtitle='1.4 GHz L!I'+greek('nu')+'!N [W/Hz]', ytitle='B!IgB!N',/xlog, thick=th, xthick=th, ythick=th, charsize=cs, /nodata, yr=[-1200,2000], /ystyle
		cgplot, /over, bz[ind_bllac].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2m)^2, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue'
		cgplot, /over, bz[ind_fsrq].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2m)^2, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red'
		cgtext, 1d26, 1700, /data, color='blue', greek('rho')+' = '+string(correlate(alog10(bz[ind_bllac].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2m)^2),bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 1d26, 1400, /data, color='red', greek('rho')+' = '+string(correlate(alog10(bz[ind_fsrq].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2m)^2),bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		al_legend, /bottom, /right, $
			bthick=th, $
			['FSRQ','BL Lac'], $
			color=['Red','Blue'], $
			psym=[7,16], $
			charsize=labelsize

		ind = where(bz.sp_mag_r gt 0.)
		ind_bllac = where(bz.sp_mag_r gt 0. and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.sp_mag_r gt 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].sp_mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind], psym=9, xtitle='M!IR!N', ytitle='B!IgB!N', thick=th, xthick=th, ythick=th, charsize=cs, /nodata, yr=[-1200,2000],/ystyle
		cgplot, /over, bz[ind_bllac].sp_mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind_bllac], psym=16, thick=th, color='blue'
		cgplot, /over, bz[ind_fsrq].sp_mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind_fsrq], psym=7, thick=th, color='red'
		cgtext, -20, 1700, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].sp_mag_r - 5*alog10(lumdist(bz[ind_bllac].z,/silent,/wmap7)*1e6) + 5.,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, -20, 1400, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].sp_mag_r - 5*alog10(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*1e6) + 5.,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		ind = where(bz.flux_xray gt 0.)
		ind_bllac = where(bz.flux_xray gt 0. and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.flux_xray gt 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind], psym=9, xtitle='(0.1-2.4) keV '+greek('nu')+'L!I'+greek('nu')+'!N [erg/s]', ytitle='B!IgB!N',/xlog, thick=th, xthick=th, ythick=th, charsize=cs, /nodata, yr=[-1200,2000], /ystyle, xr=[1d42,1d47], /xstyle
		cgplot, /over, bz[ind_bllac].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue'
		cgplot, /over, bz[ind_fsrq].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red'
		cgtext, 1d45, 1700, /data, color='blue', greek('rho')+' = '+string(correlate(alog10(bz[ind_bllac].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2cm)^2),bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 1d45, 1400, /data, color='red', greek('rho')+' = '+string(correlate(alog10(bz[ind_fsrq].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2cm)^2),bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		ind = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0.)
		ind_bllac = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0. and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].spindex_ro, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (radio-optical)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata, yr=[-1200,2000], /ystyle
		cgplot, bz[ind_bllac].spindex_ro, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
		cgplot, bz[ind_fsrq].spindex_ro, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
		cgtext, 0.60, 1700, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_ro,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 0.60, 1400, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_ro,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		ind = where(bz.spindex_ox gt -9.9 and bz.spindex_ox ne 0.)
		ind_bllac = where(bz.spindex_ox gt -9.9 and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.spindex_ox gt -9.9 and bz.spindex_ox ne 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].spindex_ox, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (optical-X-ray)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata, yr=[-1200,2000], /ystyle, xr=[0.5,2.5], /xstyle
		cgplot, bz[ind_bllac].spindex_ox, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
		cgplot, bz[ind_fsrq].spindex_ox, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
		cgtext, 1.8, 1700, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_rx,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 1.8, 1400, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_rx,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		ind = where(bz.spindex_rx gt -9.9 and bz.spindex_rx ne 0.)
		ind_bllac = where(bz.spindex_rx gt -9.9 and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.spindex_rx gt -9.9 and bz.spindex_rx ne 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].spindex_rx, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (radio-X-ray)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata, yr=[-1200,2000], /ystyle, xr=[0.4,1.0], /xstyle
		cgplot, bz[ind_bllac].spindex_rx, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
		cgplot, bz[ind_fsrq].spindex_rx, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
		cgtext, 0.80, 1700, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_rx,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 0.80, 1400, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_rx,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		if keyword_set(ps) then ps_end

		!p.multi=[0,1,1]

	endif

	if keyword_set(stop) then stop

end
