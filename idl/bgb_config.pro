;+
; NAME:
;       
;	BGB_CONFIG
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

pro bgb_config, stop=stop, noplot = noplot, ps=ps, readfiles=readfiles, multiwave=multiwave, $
	neighbors = neighbors, lowz=lowz, highz=highz, cosmo=cosmo

	blazardir = '~/Astronomy/Research/blazars/sav/'
	paperdir = '~/Astronomy/Research/blazars/paper/'

	restore,blazardir+'config_structure.sav'

	; Plot the redshift distribution of the galaxies

	ps_start, filename=paperdir+'blazar_zhist.eps',/encap,/color,/quiet
	cghistoplot, config.z, xtitle='Blazar redshift', charsize=2
	ps_end

	; total N_gal, background N_gal, radius of field in arcsec, redshift of blazar

	if n_elements(lowz) eq 0 then lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	if n_elements(highz) eq 0 then highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	if n_elements(neighbors) eq 0 then neighbors = 10
	ind = where(config.z lt highz and config.z gt lowz and config.n500_cmag ge neighbors and config.bg gt 0)

	config=config[ind]

	fr1 = where(strtrim(config.mt,2) eq 'I', fr1count)
	fr2 = where(strtrim(config.mt,2) eq 'II' or strtrim(config.mt,2) eq 'U',fr2count)
	compact = where(strtrim(config.mt,2) eq 'C' or strtrim(config.mt,2) eq 'C*' or strtrim(config.mt,2) eq 'S' or strtrim(config.mt,2) eq 'S*',compactcount)

	print,''
	print,'FR 1 galaxies: ',fr1count
	print,'FR 2 galaxies: ',fr2count
	print,'Compact galaxies: ',compactcount
	print,''

	z = config.z
	fieldsize = config.fieldsize
	nt = float(config.n500_cmag)
	nb = float(config.bg)
	countingmag = config.cmag

	bgb, nt, nb, fieldsize, z, countingmag, config_bgb, config_bgb_err, cosmo=cosmo

		kstwo, config_bgb[fr1], config_bgb[fr2], d, prob
		print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)

	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=paperdir+'bgb_config.eps', /quiet, /color, /encap
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
	cghistoplot, config_bgb, $
		datacolor='black', $
		binsize=bgbbin, $
		xtitle='B!IgB!N', $
		ytitle='N!Iblazars!N', $
		charsize=cs

	cgplots, [mean(config_bgb),mean(config_bgb)], !y.crange, linestyle=2, thick=2, color='Black'
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		'<B!IgB!N>='+string(mean(config_bgb),format='(i4)')+cgsymbol('+-')+string(stddev(config_bgb),format='(i4)'), $
		charsize=labelsize

	; Plot 2

	cghistoplot, config_bgb[fr1], $
		binsize=bgbbin, $
		yr=[0,150], $
		xtitle='B!IgB!N', $
		ytitle='N!Iblazars!N', $
		datacolor="Red",$ 
		charsize=cs

	cghistoplot,config_bgb[fr2], /oplot, $
		binsize=bgbbin, $
		datacolor="Blue"

	cghistoplot,config_bgb[compact], /oplot, $
		binsize=bgbbin, $
		datacolor="green"

	cgplots, [mean(config_bgb[fr1]),mean(config_bgb[fr1])], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, [mean(config_bgb[fr2]),mean(config_bgb[fr2])], !y.crange, linestyle=2, color="Blue", thick=2
	cgplots, [mean(config_bgb[compact]),mean(config_bgb[compact])], !y.crange, linestyle=2, color="green", thick=2
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
		color="red", $
		/data, $
		'<B!IgB!N>='+string(mean(config_bgb[fr1]),format='(i4)')+cgsymbol('+-')+string(stddev(config_bgb[fr1]),format='(i4)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		color="blue", $
		'<B!IgB!N>='+string(mean(config_bgb[fr2]),format='(i4)')+cgsymbol('+-')+string(stddev(config_bgb[fr2]),format='(i4)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.5 + !y.crange[0], $
		/data, $
		color="green", $
		'<B!IgB!N>='+string(mean(config_bgb[compact]),format='(i4)')+cgsymbol('+-')+string(stddev(config_bgb[compact]),format='(i4)'), $
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
	config_binned = fltarr(nbins-1)
	config_binned_err = fltarr(nbins-1)
	fr1_binned = fltarr(nbins-1)
	fr1_binned_err = fltarr(nbins-1)
	fr1_bgb = config_bgb[fr1]
	fr1_bgb_err = config_bgb_err[fr1]
	fr2_binned = fltarr(nbins-1)
	fr2_binned_err = fltarr(nbins-1)
	fr2_bgb = config_bgb[fr2]
	fr2_bgb_err = config_bgb_err[fr2]
	for i=0,nbins-2 do begin
		junk_ab = where(config.z gt bins[i] and config.z lt bins[i+1],ja)
		if ja gt 0 then begin
			config_binned[i] = wmean(config_bgb[junk_ab],config_bgb_err[junk_ab],error=wmean_err_ab,/nan)
			config_binned_err[i] = stddev(config_bgb[junk_ab])
			config_binned_err[i] = wmean_err_ab
		endif else print,'No radio galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_fr1 = where(config[fr1].z gt bins[i] and config[fr1].z lt bins[i+1],jf)
		if jf gt 0 then begin
			fr1_binned[i] = wmean(fr1_bgb[junk_fr1],fr1_bgb_err[junk_fr1],error=wmean_err_fr1,/nan)
			fr1_binned_err[i] = stddev(fr1_bgb[junk_fr1])
			fr1_binned_err[i] = wmean_err_fr1
		endif else print,'No FR1 galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_fr2 = where(config[fr2].z gt bins[i] and config[fr2].z lt bins[i+1],jb)
		if jb gt 0 then begin
			fr2_binned[i] = wmean(fr2_bgb[junk_fr2],fr2_bgb_err[junk_fr2],error=wmean_err_fr2,/nan)
			fr2_binned_err[i] = stddev(fr2_bgb[junk_fr2])
			fr2_binned_err[i] = wmean_err_fr2
		endif else print,'No FR2 galaxies between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
	endfor

	zarr = bins[0:nbins-2]+binsize/2.

	cgplot, indgen(10), $
		/nodata, $
		thick=th, ythick=th, xthick=th, $
		charsize=cs, $
		xrange=[0,round(highz*10.)/10.], /xstyle, $
		;yrange=[min(config_bgb),max(config_bgb)], /ystyle, $
		yrange=[-2000,2000], $
		xtitle='Blazar redshift', $
		ytitle='B!IgB!N'

	oploterror, z, config_bgb, config_bgb_err, $
		psym=3, $
		color='grey', $
		thick=th, $
		/nohat

	ab = where(config_binned ne 0.)
	oploterror, zarr[ab], config_binned[ab], config_binned_err[ab], $
		thick=th, $
		color="black", $
		/nohat, $
		psym = -16

	fb = where(fr1_binned ne 0.)
	oploterror, zarr[fb], fr1_binned[fb], fr1_binned_err[fb], $
		thick=th, $
		color="red", $
		/nohat, $
		psym = -9

	bb = where(fr2_binned ne 0.)
	oploterror, zarr[bb], fr2_binned[bb], fr2_binned_err[bb], $
		thick=th, $
		color="Blue", $
		/nohat, $
		psym = -9

	al_legend, /top, /right, $
		bthick=th, $
		['FR 1','FR 2'], $
		color=['Red','Blue'], $
		psym=[16,9], $
		charsize=labelsize

		; Fitting the raw and binned distributions with a simple linear model

		linexp = 'p[0] + x*p[1]'
		xarr = fillarr(0.01,0,1)
		bgerrind = where(finite(config_bgb_err))
		blerrind = where(finite(fr2_binned_err) and fr2_binned_err ne 0.)
		fqerrind = where(finite(fr1_binned_err) and fr1_binned_err ne 0.)
		result_all = mpfitexpr(linexp,z[bgerrind],config_bgb[bgerrind],config_bgb_err[bgerrind],perr=perr_all,/quiet)
		result_fr2 = mpfitexpr(linexp,zarr[blerrind],fr2_binned[blerrind],fr2_binned_err[blerrind],perr=perr_fr2,/quiet)
		result_fr1 = mpfitexpr(linexp,zarr[fqerrind],fr1_binned[fqerrind],fr1_binned_err[fqerrind],perr=perr_fr1,/quiet)

		;cgplot, /over, xarr, result_all[0] + xarr*result_all[1], linestyle=2, color='black'
		;cgplot, /over, xarr, result_fr2[0] + xarr*result_fr2[1], linestyle=2, color='blue'
		;cgplot, /over, xarr, result_fr1[0] + xarr*result_fr1[1], linestyle=2, color='red'


		print,'All: ',result_all, perr_all
		print,'FR 1: ',result_fr1, perr_fr1
		print,'FR 2: ',result_fr2, perr_fr2

		; Does clustering amplitude correlate with redshift?

		print,'Correlation: ',correlate(config_bgb,z)

	if keyword_set(ps) then ps_end

	!p.multi=[0,1,1]	


	endif	; noplot keyword


	if keyword_set(stop) then stop

end
