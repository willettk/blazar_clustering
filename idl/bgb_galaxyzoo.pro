
;+
; NAME:
;       
;	BGB_GALAXYZOO
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
;       Written by K. Willett               Dec 11 
;-

pro bgb_galaxyzoo, stop=stop, noplot = noplot, ps=ps, readfiles=readfiles, multiwave=multiwave, $
	neighbors = neighbors, lowz=lowz, highz=highz, cosmo=cosmo

	blazardir = '~/Astronomy/Research/blazars/sav/'
	paperdir = '~/Astronomy/Research/blazars/paper/'

	restore,blazardir+'el_modz_structure.sav'
	restore,blazardir+'sp_modz_structure.sav'

	; Plot the redshift distribution of the galaxies

	!p.multi=[0,1,1]
	ps_start, filename=paperdir+'galaxyzoo_zhist.eps',/encap,/color,/quiet
	cghistoplot, el.z, xtitle='Galaxy redshift', charsize=2, datacolor='red'
	cghistoplot, sp.z, datacolor='blue', /oplot
	al_legend,/top,/right,['Elliptical','Spiral'],psym=[28,28],color=['red','blue'], charsize=1.5
	ps_end

	; total N_gal, background N_gal, radius of field in arcsec, redshift of blazar

	if n_elements(lowz) eq 0 then lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	if n_elements(highz) eq 0 then highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	if n_elements(neighbors) eq 0 then neighbors = 10
	;ind = where(bz.z lt highz and bz.z gt lowz and bz.n500_cmag ge neighbors and bz.bg gt 0)
	spind = where(sp.z lt highz and sp.z gt lowz and sp.bg gt 0 and sp.n500_cmag ge neighbors)
	elind = where(el.z lt highz and el.z gt lowz and el.bg gt 0 and el.n500_cmag ge neighbors)

	sp=sp[spind]
	el=el[elind]

	z_sp = sp.z
	fieldsize_sp = sp.fieldsize
	nt_sp = float(sp.n500_cmag)
	nb_sp = float(sp.bg)
	countingmag_sp = sp.cmag

	z_el = el.z
	fieldsize_el = el.fieldsize
	nt_el = float(el.n500_cmag)
	nb_el = float(el.bg)
	countingmag_el = el.cmag

	z = [z_el,z_sp]
	fieldsize = [fieldsize_el,fieldsize_sp]
	nt = [nt_el,nt_sp]
	nb = [nb_el,nb_sp]
	countingmag = [countingmag_el,countingmag_sp]

	bgb, nt_el, nb_el, fieldsize_el, z_el, countingmag_el, bgb_el, bgb_el_err, cosmo=cosmo
	bgb, nt_sp, nb_sp, fieldsize_sp, z_sp, countingmag_sp, bgb_sp, bgb_sp_err, cosmo=cosmo

		kstwo, bgb_el, bgb_sp, d, prob
		print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)

	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=paperdir+'bgb_galaxyzoo.eps', /quiet, /color, /encap
		cs = 1.5
		labelsize=1
		th=2
	endif else begin
		cs = 2
		th=1
		labelsize = 1.5
	endelse
	
	; Plot the distribution of B_gB

	; Plot 1

	bgbbin = 200
	cghistoplot, [bgb_sp,bgb_el], $
		datacolor='black', $
		binsize=bgbbin, $
		xtitle='B!IgB!N', $
		ytitle='N!Igal!N', $
		charsize=cs

	cgplots, [mean([bgb_sp,bgb_el]),mean([bgb_sp,bgb_el])], !y.crange, linestyle=2, thick=2, color='Black'
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		'<B!IgB!N>='+string(mean([bgb_sp,bgb_el]),format='(i4)')+cgsymbol('+-')+string(stddev([bgb_sp,bgb_el]),format='(i4)'), $
		charsize=labelsize

	; Plot 2

	cghistoplot, bgb_el, $
		binsize=bgbbin, $
		yr=[0,150], $
		xtitle='B!IgB!N', $
		ytitle='N!Iblazars!N', $
		datacolor="Red",$ 
		charsize=cs

	cghistoplot,bgb_sp, /oplot, $
		binsize=bgbbin, $
		datacolor="Blue"

	cgplots, [mean(bgb_el),mean(bgb_el)], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, [mean(bgb_sp),mean(bgb_sp)], !y.crange, linestyle=2, color="Blue", thick=2
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
		color="red", $
		/data, $
		'<B!IgB!N>='+string(mean(bgb_el),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_el),format='(i4)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		color="blue", $
		'<B!IgB!N>='+string(mean(bgb_sp),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_sp),format='(i4)'), $
		charsize=labelsize

	; Plot 3

	cgplot, z, nt, $
		/nodata, $
		charsize=cs, $
		xrange=[0,round(highz*10.)/10.], /xstyle, $
		yrange=[(min(nb) < min(nt) < min(nt-nb)),(max(nb) > max(nt) > max(nt-nb))], $
		xtitle='z', $
		ytitle='Nt - Nb'

	cgplot, z, float(nt)-float(nb), psym=9, /over,color='black'
	cgplot, z, float(nb), psym=15, /over, color='pink'
	cgplot, z, float(nt), psym=7, /over, color='green'

	; Plot 4 - binned Bgb as fn. of redshift

	bgb_all = [bgb_el,bgb_sp]
	bgb_all_err = [bgb_el_err,bgb_sp_err]

	binsize = 0.05
	bins = fillarr(binsize,lowz,highz+binsize)
	nbins = n_elements(bins)
	ab_binned = fltarr(nbins-1)
	ab_binned_err = fltarr(nbins-1)
	el_binned = fltarr(nbins-1)
	el_binned_err = fltarr(nbins-1)
	sp_binned = fltarr(nbins-1)
	sp_binned_err = fltarr(nbins-1)
	el_bgb = bgb_el
	sp_bgb = bgb_sp
	el_bgb_err = bgb_el_err
	sp_bgb_err = bgb_sp_err
	for i=0,nbins-2 do begin
		junk_ab = where(z gt bins[i] and z lt bins[i+1],ja)
		if ja gt 0 then begin
			ab_binned[i] = wmean(bgb_all[junk_ab],bgb_all_err[junk_ab],error=wmean_err_ab,/nan)
			ab_binned_err[i] = stddev(bgb_all[junk_ab])
			ab_binned_err[i] = wmean_err_ab
		endif else print,'No galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_el = where(el.z gt bins[i] and el.z lt bins[i+1],jf)
		if jf gt 0 then begin
			el_binned[i] = wmean(el_bgb[junk_el],el_bgb_err[junk_el],error=wmean_err_el,/nan)
			el_binned_err[i] = stddev(el_bgb[junk_el])
			el_binned_err[i] = wmean_err_el
		endif else print,'No ellipticals between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_sp = where(sp.z gt bins[i] and sp.z lt bins[i+1],jb)
		if jb gt 0 then begin
			sp_binned[i] = wmean(sp_bgb[junk_sp],sp_bgb_err[junk_sp],error=wmean_err_sp,/nan)
			sp_binned_err[i] = stddev(sp_bgb[junk_sp])
			sp_binned_err[i] = wmean_err_sp
		endif else print,'No spirals between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
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

	oploterror, z, bgb_all, bgb_all_err, $
		psym=3, $
		color='grey', $
		thick=th, $
		/nohat

	ab = where(ab_binned ne 0.)
	oploterror, zarr[ab], ab_binned[ab], ab_binned_err[ab], $
		thick=th, $
		color="black", $
		/nohat, $
		psym = -16

	fb = where(el_binned ne 0.)
	oploterror, zarr[fb], el_binned[fb], el_binned_err[fb], $
		thick=th, $
		color="red", $
		/nohat, $
		psym = -9

	bb = where(sp_binned ne 0.)
	oploterror, zarr[bb], sp_binned[bb], sp_binned_err[bb], $
		thick=th, $
		color="Blue", $
		/nohat, $
		psym = -9

	al_legend, /top, /right, $
		bthick=th, $
		['Elliptical','Spiral'], $
		color=['Red','Blue'], $
		psym=[16,9], $
		charsize=labelsize

		; Fitting the raw and binned distributions with a simple linear model


		linexp = 'p[0] + x*p[1]'
		xarr = fillarr(0.01,0,1)
		bgerrind = where(finite(bgb_all_err))
		blerrind = where(finite(sp_binned_err) and sp_binned_err ne 0.)
		fqerrind = where(finite(el_binned_err) and el_binned_err ne 0.)
		result_all = mpfitexpr(linexp,z[bgerrind],bgb_all[bgerrind],bgb_all_err[bgerrind],perr=perr_all,/quiet)
		result_sp = mpfitexpr(linexp,zarr[blerrind],sp_binned[blerrind],sp_binned_err[blerrind],perr=perr_sp,/quiet)
		result_el = mpfitexpr(linexp,zarr[fqerrind],el_binned[fqerrind],el_binned_err[fqerrind],perr=perr_el,/quiet)

		;cgplot, /over, xarr, result_all[0] + xarr*result_all[1], linestyle=2, color='black'
		;cgplot, /over, xarr, result_bllac[0] + xarr*result_bllac[1], linestyle=2, color='blue'
		;cgplot, /over, xarr, result_fsrq[0] + xarr*result_fsrq[1], linestyle=2, color='red'


		print,'All galaxies: ',result_all, perr_all
		print,'Spiral: ',result_sp, perr_sp
		print,'Elliptical: ',result_el, perr_el


		; Does clustering amplitude correlate with redshift?

		print,'Correlation: ',correlate(bgb_all,z)


	;ploterror, z[fsrq], bgb_el, bgb_blazar_err[fsrq], $
	;	color="Red", $
	;	errcolor="red", $
	;	charsize=cs, $
	;	xrange=[0,0.8], /xstyle, $
	;	/nohat, $
	;	xtitle='Redshift', $
	;	ytitle='B!IgB!N', $
	;	psym = 16

	;oploterror, z[bllac], bgb_sp, bgb_blazar_err[bllac], $
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
		cgtext, 2d23, 1300, /data, color='blue', greek('rho')+' = '+string(correlate(alog10(bz[ind_bllac].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2m)^2),bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 2d23, 1200, /data, color='red', greek('rho')+' = '+string(correlate(alog10(bz[ind_fsrq].flux_radio*mjy2si*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2m)^2),bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		al_legend, /top, /right, $
			bthick=th, $
			['FSRQ','BL Lac'], $
			color=['Red','Blue'], $
			psym=[7,16], $
			charsize=labelsize

		ind = where(bz.sp_mag_r gt 0.)
		ind_bllac = where(bz.sp_mag_r gt 0. and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.sp_mag_r gt 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].sp_mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind], psym=9, xtitle='M!IR!N', ytitle='B!IgB!N', thick=th, xthick=th, ythick=th, charsize=cs, /nodata
		cgplot, /over, bz[ind_bllac].sp_mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind_bllac], psym=16, thick=th, color='blue'
		cgplot, /over, bz[ind_fsrq].sp_mag_r - 5*alog10(lumdist(bz[ind].z,/silent,/wmap7)*1e6) + 5., bgb_blazar[ind_fsrq], psym=7, thick=th, color='red'
		cgtext, -27.5, 1300, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].sp_mag_r - 5*alog10(lumdist(bz[ind_bllac].z,/silent,/wmap7)*1e6) + 5.,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, -27.5, 1200, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].sp_mag_r - 5*alog10(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*1e6) + 5.,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		ind = where(bz.flux_xray gt 0.)
		ind_bllac = where(bz.flux_xray gt 0. and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.flux_xray gt 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind], psym=9, xtitle='(0.1-2.4) keV '+greek('nu')+'L!I'+greek('nu')+'!N [erg/s]', ytitle='B!IgB!N',/xlog, thick=th, xthick=th, ythick=th, charsize=cs, /nodata
		cgplot, /over, bz[ind_bllac].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue'
		cgplot, /over, bz[ind_fsrq].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red'
		cgtext, 2d42, 1300, /data, color='blue', greek('rho')+' = '+string(correlate(alog10(bz[ind_bllac].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_bllac].z,/silent,/wmap7)*mpc2cm)^2),bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 2d42, 1200, /data, color='red', greek('rho')+' = '+string(correlate(alog10(bz[ind_fsrq].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind_fsrq].z,/silent,/wmap7)*mpc2cm)^2),bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		ind = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0.)
		ind_bllac = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0. and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].spindex_ro, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (radio-optical)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata
		cgplot, bz[ind_bllac].spindex_ro, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
		cgplot, bz[ind_fsrq].spindex_ro, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
		cgtext, 0.05, 1300, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_ro,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 0.05, 1200, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_ro,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		ind = where(bz.spindex_ox gt -9.9 and bz.spindex_ox ne 0.)
		ind_bllac = where(bz.spindex_ox gt -9.9 and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.spindex_ox gt -9.9 and bz.spindex_ox ne 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].spindex_ox, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (optical-X-ray)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata
		cgplot, bz[ind_bllac].spindex_ox, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
		cgplot, bz[ind_fsrq].spindex_ox, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
		cgtext, 0.60, 1300, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_rx,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 0.60, 1200, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_rx,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		ind = where(bz.spindex_rx gt -9.9 and bz.spindex_rx ne 0.)
		ind_bllac = where(bz.spindex_rx gt -9.9 and bz.btype eq 'BLLac')
		ind_fsrq = where(bz.spindex_rx gt -9.9 and bz.spindex_rx ne 0. and bz.btype eq 'FSRQ')
		cgplot, bz[ind].spindex_rx, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+' (radio-X-ray)', ytitle='B!IgB!N', charsize=cs, thick=th, xthick=th, ythick=th, /nodata
		cgplot, bz[ind_bllac].spindex_rx, bgb_blazar[ind_bllac], psym=16, thick=th, color='blue', /over
		cgplot, bz[ind_fsrq].spindex_rx, bgb_blazar[ind_fsrq], psym=7, thick=th, color='red', /over
		cgtext, 0.45, 1300, /data, color='blue', greek('rho')+' = '+string(correlate(bz[ind_bllac].spindex_rx,bgb_blazar[ind_bllac]),format='(f5.2)'),charsize=labelsize
		cgtext, 0.45, 1200, /data, color='red', greek('rho')+' = '+string(correlate(bz[ind_fsrq].spindex_rx,bgb_blazar[ind_fsrq]),format='(f5.2)'),charsize=labelsize

		!p.multi=[0,1,1]

	endif

	if keyword_set(stop) then stop

end
