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
; REVISION HISTORY
;       Written by K. Willett                Apr 12
;-

pro bgb_ki, stop=stop, noplot = noplot, ps=ps, readfiles=readfiles, multiwave=multiwave, $
	neighbors = neighbors, lowz=lowz, highz=highz, cosmo=cosmo

	blazardir = '~/Astronomy/Research/blazars/sav/'
	paperdir = '~/Astronomy/Research/blazars/paper/'

	restore,blazardir+'ki_modz_structure.sav'

	; Plot the redshift distribution of the galaxies

	ps_start, filename=paperdir+'blazar_zhist.eps',/encap,/color,/quiet
	cghistoplot, ki.z, xtitle='Blazar redshift', charsize=2
	ps_end

	; total N_gal, background N_gal, radius of field in arcsec, redshift of blazar

	if n_elements(lowz) eq 0 then lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	if n_elements(highz) eq 0 then highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	if n_elements(neighbors) eq 0 then neighbors = 10
	ind = where(ki.z lt highz and ki.z gt lowz and ki.n500_cmag ge neighbors and ki.bg gt 0)

	ki=ki[ind]

	mcompact = where(ki.radio_morph eq 'COMPACT')
	mcomplex = where(ki.radio_morph eq 'COMPLEX')
	mresolved = where(ki.radio_morph eq 'RESOLVED')

	z = ki.z
	fieldsize = ki.fieldsize
	nt = float(ki.n500_cmag)
	nb = float(ki.bg)
	countingmag = ki.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bgb_blazar, bgb_blazar_err, cosmo=cosmo

		kstwo, bgb_blazar[mcomplex], bgb_blazar[mcompact], d, prob
		print,'KS result (complex-compact):',d,prob,sqrt(2d) * inverf(1d - prob)
		kstwo, bgb_blazar[mcomplex], bgb_blazar[mresolved], d, prob
		print,'KS result (complex-resolved):',d,prob,sqrt(2d) * inverf(1d - prob)
		kstwo, bgb_blazar[mresolved], bgb_blazar[mcompact], d, prob
		print,'KS result (resolved-compact):',d,prob,sqrt(2d) * inverf(1d - prob)

	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=paperdir+'bgb_blazars.eps', /quiet, /color, /encap
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
	cghistoplot, bgb_blazar, $
		datacolor='black', $
		binsize=bgbbin, $
		xtitle='B!IgB!N', $
		ytitle='N!Iblazars!N', $
		charsize=cs

	cgplots, [mean(bgb_blazar),mean(bgb_blazar)], !y.crange, linestyle=2, thick=2, color='Black'
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		'<B!IgB!N>='+string(mean(bgb_blazar),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar),format='(i4)'), $
		charsize=labelsize

	; Plot 2

	cghistoplot, bgb_blazar[mcomplex], $
		binsize=bgbbin, $
		yr=[0,150], $
		xtitle='B!IgB!N', $
		ytitle='N!Iblazars!N', $
		datacolor="Red",$ 
		charsize=cs

	cghistoplot,bgb_blazar[mcompact], /oplot, $
		binsize=bgbbin, $
		datacolor="Blue"

	cghistoplot,bgb_blazar[mresolved], /oplot, $
		binsize=bgbbin, $
		datacolor="Green"

	cgplots, [mean(bgb_blazar[mcomplex]),mean(bgb_blazar[mcomplex])], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, [mean(bgb_blazar[mcompact]),mean(bgb_blazar[mcompact])], !y.crange, linestyle=2, color="Blue", thick=2
	cgplots, [mean(bgb_blazar[mresolved]),mean(bgb_blazar[mresolved])], !y.crange, linestyle=2, color="Green", thick=2
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
		color="red", $
		/data, $
		'<B!IgB!N>='+string(mean(bgb_blazar[mcomplex]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar[mcomplex]),format='(i4)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		color="blue", $
		'<B!IgB!N>='+string(mean(bgb_blazar[mcompact]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar[mcompact]),format='(i4)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.5 + !y.crange[0], $
		/data, $
		color="green", $
		'<B!IgB!N>='+string(mean(bgb_blazar[mresolved]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar[mresolved]),format='(i4)'), $
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

	; Plot binned Bgb as fn. of redshift

	binsize = 0.05
	bins = fillarr(binsize,lowz,highz+binsize)
	nbins = n_elements(bins)
	ki_binned = fltarr(nbins-1)
	ki_binned_err = fltarr(nbins-1)
	mcomplex_binned = fltarr(nbins-1)
	mcomplex_binned_err = fltarr(nbins-1)
	mcomplex_bgb = bgb_blazar[mcomplex]
	mcomplex_bgb_err = bgb_blazar_err[mcomplex]
	mcompact_binned = fltarr(nbins-1)
	mcompact_binned_err = fltarr(nbins-1)
	mcompact_bgb = bgb_blazar[mcompact]
	mcompact_bgb_err = bgb_blazar_err[mcompact]
	mresolved_binned = fltarr(nbins-1)
	mresolved_binned_err = fltarr(nbins-1)
	mresolved_bgb = bgb_blazar[mresolved]
	mresolved_bgb_err = bgb_blazar_err[mresolved]
	for i=0,nbins-2 do begin
		junk_ab = where(ki.z gt bins[i] and ki.z lt bins[i+1],ja)
		if ja gt 0 then begin
			ki_binned[i] = wmean(bgb_blazar[junk_ab],bgb_blazar_err[junk_ab],error=wmean_err_ab,/nan)
			ki_binned_err[i] = stddev(bgb_blazar[junk_ab])
			ki_binned_err[i] = wmean_err_ab
		endif else print,'No radio galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_mcomplex = where(ki[mcomplex].z gt bins[i] and ki[mcomplex].z lt bins[i+1],jf)
		if jf gt 0 then begin
			mcomplex_binned[i] = wmean(mcomplex_bgb[junk_mcomplex],mcomplex_bgb_err[junk_mcomplex],error=wmean_err_mcomplex,/nan)
			mcomplex_binned_err[i] = stddev(mcomplex_bgb[junk_mcomplex])
			mcomplex_binned_err[i] = wmean_err_mcomplex
		endif else print,'No complex galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_mcompact = where(ki[mcompact].z gt bins[i] and ki[mcompact].z lt bins[i+1],jb)
		if jb gt 0 then begin
			mcompact_binned[i] = wmean(mcompact_bgb[junk_mcompact],mcompact_bgb_err[junk_mcompact],error=wmean_err_mcompact,/nan)
			mcompact_binned_err[i] = stddev(mcompact_bgb[junk_mcompact])
			mcompact_binned_err[i] = wmean_err_mcompact
		endif else print,'No compact galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_mresolved = where(ki[mresolved].z gt bins[i] and ki[mresolved].z lt bins[i+1],jr)
		if jr gt 0 then begin
			mresolved_binned[i] = wmean(mresolved_bgb[junk_mresolved],mresolved_bgb_err[junk_mresolved],error=wmean_err_mresolved,/nan)
			mresolved_binned_err[i] = stddev(mresolved_bgb[junk_mresolved])
			mresolved_binned_err[i] = wmean_err_mresolved
		endif else print,'No resolved galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
	endfor

	zarr = bins[0:nbins-2]+binsize/2.

	cgplot, indgen(10), $
		/nodata, $
		thick=th, ythick=th, xthick=th, $
		charsize=cs, $
		xrange=[0,round(highz*10.)/10.], /xstyle, $
		;yrange=[min(bgb_blazar),max(bgb_blazar)], /ystyle, $
		yrange=[-2000,2000], $
		xtitle='Blazar redshift', $
		ytitle='B!IgB!N'

	oploterror, z, bgb_blazar, bgb_blazar_err, $
		psym=3, $
		color='grey', $
		thick=th, $
		/nohat

	ab = where(ki_binned ne 0.)
	oploterror, zarr[ab], ki_binned[ab], ki_binned_err[ab], $
		thick=th, $
		color="black", $
		/nohat, $
		psym = -16

	fb = where(mcomplex_binned ne 0.)
	oploterror, zarr[fb], mcomplex_binned[fb], mcomplex_binned_err[fb], $
		thick=th, $
		color="red", $
		/nohat, $
		psym = -9

	bb = where(mcompact_binned ne 0.)
	oploterror, zarr[bb], mcompact_binned[bb], mcompact_binned_err[bb], $
		thick=th, $
		color="Blue", $
		/nohat, $
		psym = -9

	rb = where(mresolved_binned ne 0.)
	oploterror, zarr[rb], mcompact_binned[rb], mcompact_binned_err[rb], $
		thick=th, $
		color="Green", $
		/nohat, $
		psym = -9

	al_legend, /top, /right, $
		bthick=th, $
		['Complex','Compact','Resolved'], $
		color=['Red','Blue','Green'], $
		psym=[16,9,9], $
		charsize=labelsize

		; Fitting the raw and binned distributions with a simple linear model


		linexp = 'p[0] + x*p[1]'
		xarr = fillarr(0.01,0,1)
		bgerrind = where(finite(bgb_blazar_err))
		blerrind = where(finite(mcompact_binned_err) and mcompact_binned_err ne 0.)
		fqerrind = where(finite(mcomplex_binned_err) and mcomplex_binned_err ne 0.)
		result_all = mpfitexpr(linexp,z[bgerrind],bgb_blazar[bgerrind],bgb_blazar_err[bgerrind],perr=perr_all,/quiet)
		result_mcompact = mpfitexpr(linexp,zarr[blerrind],mcompact_binned[blerrind],mcompact_binned_err[blerrind],perr=perr_mcompact,/quiet)
		result_mcomplex = mpfitexpr(linexp,zarr[fqerrind],mcomplex_binned[fqerrind],mcomplex_binned_err[fqerrind],perr=perr_mcomplex,/quiet)

		;cgplot, /over, xarr, result_all[0] + xarr*result_all[1], linestyle=2, color='black'
		;cgplot, /over, xarr, result_mcompact[0] + xarr*result_mcompact[1], linestyle=2, color='blue'
		;cgplot, /over, xarr, result_mcomplex[0] + xarr*result_mcomplex[1], linestyle=2, color='red'


		print,'All: ',result_all, perr_all
		print,'Compact: ',result_mcompact, perr_mcompact
		print,'Complex: ',result_mcomplex, perr_mcomplex

		; Does clustering amplitude correlate with redshift?

		print,'Correlation: ',correlate(bgb_blazar,z)

	if keyword_set(ps) then ps_end

	!p.multi=[0,1,1]	


	endif	; noplot keyword


	if keyword_set(stop) then stop

end
