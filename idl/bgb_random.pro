
;+
; NAME:
;       
;	BGB_RANDOM
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

pro bgb_random, stop=stop, noplot = noplot, ps=ps, readfiles=readfiles, multiwave=multiwave, $
	neighbors = neighbors, lowz=lowz, highz=highz, cosmo=cosmo

	blazardir = '~/Astronomy/Research/blazars/sav/'
	paperdir = '~/Astronomy/Research/blazars/paper/'

	restore,blazardir+'random_structure.sav'

	; total N_gal, background N_gal, radius of field in arcsec, redshift of blazar (matched to random location)

	if n_elements(lowz) eq 0 then lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	if n_elements(highz) eq 0 then highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	if n_elements(neighbors) eq 0 then neighbors = 10
	ind = where(random.z lt highz and random.z gt lowz and random.n500_cmag ge neighbors and random.bg gt 0)

	random=random[ind]

	z = random.z
	fieldsize = random.fieldsize
	nt = float(random.n500_cmag)
	nb = float(random.bg)
	countingmag = random.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bgb_random, bgb_random_err, cosmo=cosmo

		;kstwo, bgb_random[fsrq], bgb_random[bllac], d, prob
		;print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)

	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=paperdir+'bgb_random.eps', /quiet, /color, /encap
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
	cghistoplot, bgb_random, $
		datacolor='black', $
		binsize=bgbbin, $
		xtitle='B!IgB!N', $
		ytitle='N!Igal!N', $
		charsize=cs

	cgplots, [mean(bgb_random),mean(bgb_random)], !y.crange, linestyle=2, thick=2, color='Black'
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		'<B!IgB!N>='+string(mean(bgb_random),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_random),format='(i4)'), $
		charsize=labelsize

	; Plot 2 - split randomly

	nrandom = n_elements(random)
	sample1 = cgrandomindices(nrandom,nrandom/2)
	sample2 = setdifference(indgen(nrandom),sample1)

	cghistoplot, bgb_random[sample1], $
		binsize=bgbbin, $
		yr=[0,150], $
		xtitle='B!IgB!N', $
		ytitle='N!Igal!N', $
		datacolor="Red",$ 
		charsize=cs

	cghistoplot,bgb_random[sample2], /oplot, $
		binsize=bgbbin, $
		datacolor="Blue"

	cgplots, [mean(bgb_random[sample1]),mean(bgb_random[sample1])], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, [mean(bgb_random[sample2]),mean(bgb_random[sample2])], !y.crange, linestyle=2, color="Blue", thick=2
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
		color="red", $
		/data, $
		'<B!IgB!N>='+string(mean(bgb_random[sample1]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_random[sample1]),format='(i4)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		color="blue", $
		'<B!IgB!N>='+string(mean(bgb_random[sample2]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_random[sample2]),format='(i4)'), $
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
	ab_binned = fltarr(nbins-1)
	ab_binned_err = fltarr(nbins-1)
	;fsrq_binned = fltarr(nbins-1)
	;fsrq_binned_err = fltarr(nbins-1)
	;bllac_binned = fltarr(nbins-1)
	;bllac_binned_err = fltarr(nbins-1)
	;fsrq_bgb = bgb_random[fsrq]
	;bllac_bgb = bgb_random[bllac]
	;fsrq_bgb_err = bgb_random_err[fsrq]
	;bllac_bgb_err = bgb_random_err[bllac]
	for i=0,nbins-2 do begin
		junk_ab = where(random.z gt bins[i] and random.z lt bins[i+1],ja)
		if ja gt 0 then begin
			ab_binned[i] = wmean(bgb_random[junk_ab],bgb_random_err[junk_ab],error=wmean_err_ab,/nan)
			ab_binned_err[i] = stddev(bgb_random[junk_ab])
			ab_binned_err[i] = wmean_err_ab
		endif else print,'No galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		;junk_fsrq = where(random[fsrq].z gt bins[i] and random[fsrq].z lt bins[i+1],jf)
		;if jf gt 0 then begin
		;	fsrq_binned[i] = wmean(fsrq_bgb[junk_fsrq],fsrq_bgb_err[junk_fsrq],error=wmean_err_fsrq,/nan)
		;	fsrq_binned_err[i] = stddev(fsrq_bgb[junk_fsrq])
		;	fsrq_binned_err[i] = wmean_err_fsrq
		;endif else print,'No FSRQs between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		;junk_bllac = where(random[bllac].z gt bins[i] and random[bllac].z lt bins[i+1],jb)
		;if jb gt 0 then begin
		;	bllac_binned[i] = wmean(bllac_bgb[junk_bllac],bllac_bgb_err[junk_bllac],error=wmean_err_bllac,/nan)
		;	bllac_binned_err[i] = stddev(bllac_bgb[junk_bllac])
		;	bllac_binned_err[i] = wmean_err_bllac
		;endif else print,'No BL Lacs between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
	endfor

	zarr = bins[0:nbins-2]+binsize/2.

	cgplot, indgen(10), $
		/nodata, $
		thick=th, ythick=th, xthick=th, $
		charsize=cs, $
		xrange=[0,round(highz*10.)/10.], /xstyle, $
;		yrange=[min(bgb_random),max(bgb_random)], /ystyle, $
		yr=[-2000,2000], $
		xtitle='Redshift', $
		ytitle='B!IgB!N'

	oploterror, z, bgb_random, bgb_random_err, $
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

	;fb = where(fsrq_binned ne 0.)
	;oploterror, zarr[fb], fsrq_binned[fb], fsrq_binned_err[fb], $
	;	thick=th, $
	;	color="red", $
	;	/nohat, $
	;	psym = -9

	;bb = where(bllac_binned ne 0.)
	;oploterror, zarr[bb], bllac_binned[bb], bllac_binned_err[bb], $
	;	thick=th, $
	;	color="Blue", $
	;	/nohat, $
	;	psym = -9

	;al_legend, /top, /right, $
	;	bthick=th, $
	;	['FSRQ','BL Lac'], $
	;	color=['Red','Blue'], $
	;	psym=[16,9], $
	;	charsize=labelsize

		; Fitting the raw and binned distributions with a simple linear model


		linexp = 'p[0] + x*p[1]'
		xarr = fillarr(0.01,0,1)
		bgerrind = where(finite(bgb_random_err))
		;blerrind = where(finite(bllac_binned_err) and bllac_binned_err ne 0.)
		;fqerrind = where(finite(fsrq_binned_err) and fsrq_binned_err ne 0.)
		result_all = mpfitexpr(linexp,z[bgerrind],bgb_random[bgerrind],bgb_random_err[bgerrind],perr=perr_all,/quiet)
		;result_bllac = mpfitexpr(linexp,zarr[blerrind],bllac_binned[blerrind],bllac_binned_err[blerrind],perr=perr_bllac,/quiet)
		;result_fsrq = mpfitexpr(linexp,zarr[fqerrind],fsrq_binned[fqerrind],fsrq_binned_err[fqerrind],perr=perr_fsrq,/quiet)

		;cgplot, /over, xarr, result_all[0] + xarr*result_all[1], linestyle=2, color='black'
		;cgplot, /over, xarr, result_bllac[0] + xarr*result_bllac[1], linestyle=2, color='blue'
		;cgplot, /over, xarr, result_fsrq[0] + xarr*result_fsrq[1], linestyle=2, color='red'


		print,'All: ',result_all, perr_all
		;print,'BL Lac: ',result_bllac, perr_bllac
		;print,'FSRQ: ',result_fsrq, perr_fsrq


		; Does clustering amplitude correlate with redshift?

		print,'Correlation: ',correlate(bgb_random,z)

	;ploterror, z[fsrq], bgb_random[fsrq], bgb_random_err[fsrq], $
	;	color="Red", $
	;	errcolor="red", $
	;	charsize=cs, $
	;	xrange=[0,0.8], /xstyle, $
	;	/nohat, $
	;	xtitle='Redshift', $
	;	ytitle='B!IgB!N', $
	;	psym = 16

	;oploterror, z[bllac], bgb_random[bllac], bgb_random_err[bllac], $
	;	color="Blue", $
	;	/nohat, $
	;	psym = 9

	if keyword_set(ps) then ps_end

	!p.multi=[0,1,1]	


	endif	; noplot keyword

	if keyword_set(stop) then stop

end
