
;+
; NAME:
;       
;	BGB_TCRR
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

pro bgb_tcrr, stop=stop, noplot = noplot, ps=ps

	blazardir = '~/Astronomy/Research/blazars/sav/'
	paperdir = '~/Astronomy/Research/blazars/paper/'

	restore,blazardir+'tcrr_structure.sav'
	annulus_area, 22.2, nbg_per_arcmin2, /quiet, /noplot

	; Plot the redshift distribution of the galaxies

	cgPS_open, filename=paperdir+'tcrr_zhist.eps',/color,/quiet,/encap
	cghistoplot, tcrr.z, xtitle='3CRR redshift', charsize=2
	cgPS_close

	ind = where(tcrr.z lt 1.0)
	tcrr=tcrr[ind]

	fr1 = where(tcrr.frclass eq 1)
	fr2 = where(tcrr.frclass eq 2)

	; total N_gal, background N_gal, radius of field in arcsec, redshift of radio galaxy

	annulus_area, tcrr.cmag_app, ng_per_arcmin2, /quiet, /noplot

	fieldradius = 500	; kpc
	z = tcrr.z
	fieldsize = zang(float(fieldradius),z,/silent,/wmap7)
	nt = float(tcrr.n500_cmag)
	nb = !pi * (fieldsize/60.)^2 * ng_per_arcmin2
	countingmag = tcrr.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bgb_tcrr, bgb_tcrr_err

		kstwo, bgb_tcrr[fr1], bgb_tcrr[fr2], d, prob
		print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)

	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		cgPS_open, filename=paperdir+'bgb_tcrr.eps', /quiet, /color,/encap
		cs = 1.5
		labelsize=1
	endif else begin
		cs = 2
		labelsize = 1.5
	endelse
	
	; Plot the distribution of B_gB

	cghistoplot, bgb_tcrr, $
		binsize=500, $
		xr=[-3000,6000], $
		xtitle='B!IgB!N', $
		ytitle='Number of radio galaxies', $
		charsize=cs

	cgplots, /data, [median(bgb_tcrr),median(bgb_tcrr)], !y.crange, linestyle=2, thick=2, color='Black'
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.0 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		'<B!IgB!N>='+string(median(bgb_tcrr),format='(i5)')+cgsymbol('+-')+string(stddev(bgb_tcrr),format='(i5)'), $
		charsize=labelsize

	bgbbin = 200
	cghistoplot, bgb_tcrr[fr1], $
		binsize=bgbbin, $
		xr=[-300,2000], $
		yr=[0,50], $
		xtitle='B!IgB!N', $
		ytitle='Number of radio galaxies', $
		datacolor="Red",$ 
		charsize=cs

	cghistoplot,bgb_tcrr[fr2], /oplot, $
		binsize=bgbbin, $
		datacolor="Blue"

	cgplots, /data, [median(bgb_tcrr[fr1]),median(bgb_tcrr[fr1])], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, /data, [median(bgb_tcrr[fr2]),median(bgb_tcrr[fr2])], !y.crange, linestyle=2, color="Blue", thick=2
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
		color="red", $
		/data, $
		'<B!IgB!N>='+string(median(bgb_tcrr[fr1]),format='(i5)')+cgsymbol('+-')+string(stddev(bgb_tcrr[fr1]),format='(i5)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		color="blue", $
		'<B!IgB!N>='+string(median(bgb_tcrr[fr2]),format='(i5)')+cgsymbol('+-')+string(stddev(bgb_tcrr[fr2]),format='(i5)'), $
		charsize=labelsize

	; Plot Bgb as fn. of redshift

;	ploterror, z, bgb_tcrr, bgb_tcrr_err, $
;		charsize=cs, $
;		/nohat, $
;		yr=[-1d3,3d3], $
;		xrange=[0,1], $
;		xtitle='Redshift', $
;		ytitle='B!IgB!N', $
;		psym = 9
;
	cgplot, z, float(nt)-float(nb), psym=9, $
		charsize=cs, $
		xtitle='z', $
		ytitle='Nt - Nb'

	ploterror, z[fr1], bgb_tcrr[fr1], bgb_tcrr_err[fr1], $
		color="Red", $
		errcolor="red", $
		charsize=cs, $
		xrange=[0,1], $
		/nohat, $
		yr=[-1d3,3d3], $
		xtitle='Redshift', $
		ytitle='B!IgB!N', $
		psym = 16

	oploterror, z[fr2], bgb_tcrr[fr2], bgb_tcrr_err[fr2], $
		color="Blue", $
		/nohat, $
		psym = 9

	al_legend, /top, /left, $
		['FR I','FR II'], $
		color=['Red','Blue'], $
		psym=[16,9], $
		charsize=labelsize

	if keyword_set(ps) then cgPS_close

	!p.multi=[0,1,1]	

	endif

	if keyword_set(stop) then stop

end

