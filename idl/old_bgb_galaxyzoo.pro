
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

pro bgb_galaxyzoo, stop=stop, noplot = noplot, ps=ps

	blazardir = '~/Astronomy/Research/blazars/'
	paperdir = '~/Astronomy/Research/blazars/paper/'

	restore,blazardir+'elliptical_structure.sav'
	restore,blazardir+'spiral_structure.sav'
	annulus_area, nbg_per_arcmin2, /quiet, /noplot

	; Establish a limit on redshift due to angular size of the SDSS search,
	; while also eliminating objects in Stripe 82

	anglim_e = where(elliptical.z gt 0.162 and $
		(elliptical.ra gt 59 or $
		elliptical.ra lt -50) and $
		(elliptical.dec gt 1.25 or $
		elliptical.dec lt -1.25))
	anglim_s = where(spiral.z gt 0.162 and $
		(spiral.ra gt 59 or $
		spiral.ra lt -50) and $
		(spiral.dec gt 1.25 or $
		spiral.dec lt -1.25))

	elliptical = elliptical[anglim_e]
	spiral = spiral[anglim_s]

	eind = indgen(n_elements(elliptical.z))
	sind = indgen(n_elements(spiral.z))+n_elements(eind)

	; Plot the redshift distribution of the galaxies

	ps_start, filename=paperdir+'galaxyzoo_zhist.eps',/encap,/color,/quiet
	cghistoplot, spiral.z, xtitle='Redshift', charsize=2, binsize=0.02
	cghistoplot, elliptical.z, /oplot, datacolor="Blue", binsize=0.02
	al_legend, /top, /right, ['Elliptical','Spiral'], linestyle=0, color=['Red','Blue'], charsize=1
	ps_end

	; total N_gal, background N_gal, radius of field in arcsec, redshift of GZ object

	fieldradius = 500	; kpc
	z = [elliptical.z,spiral.z]
	fieldsize = zang(float(fieldradius),z,/silent,/wmap7)
	nt = [elliptical.n500,spiral.n500]
	nbg = !pi * (fieldsize/60.)^2 * nbg_per_arcmin2

	bgb, nt, nbg, fieldsize, z, bgb_galaxyzoo, bgb_galaxyzoo_err

		kstwo, bgb_galaxyzoo[eind], bgb_galaxyzoo[sind], d, prob
		print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)
		
	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=paperdir+'bgb_galaxyzoo.eps', /quiet, /color, /encap
		cs = 1.5
		labelsize=1
	endif else begin
		cs = 2
		labelsize = 1.5
	endelse
	
	; Plot the distribution of B_gB

	bgbbin = 5000
	cghistoplot, bgb_galaxyzoo, $
		binsize=bgbbin, $
		;xr=[-3000,6000], $
		xtitle='B!IgB!N', $
		ytitle='Number of GZ objects', $
		charsize=cs

	vline, /data, median(bgb_galaxyzoo), linestyle=2, thick=2, color='Black'
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.0 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		'<B!IgB!N>='+string(median(bgb_galaxyzoo),format='(i6)')+cgsymbol('+-')+string(stddev(bgb_galaxyzoo),format='(i6)'), $
		charsize=labelsize

	cghistoplot, bgb_galaxyzoo[eind], $
		binsize=bgbbin, $
		;xr=[-1000,9000], $
		xr=[min(bgb_galaxyzoo),max(bgb_galaxyzoo)], $
		;yr=[0,50], $
		xtitle='B!IgB!N', $
		ytitle='Number of GZ objects', $
		datacolor="Red",$ 
		charsize=cs

	cghistoplot,bgb_galaxyzoo[sind], /oplot, $
		binsize=bgbbin, $
		datacolor="Blue"

	vline, /data, median(bgb_galaxyzoo[eind]), linestyle=2, color="Red", thick=2, /noerase
	vline, /data, median(bgb_galaxyzoo[sind]), linestyle=2, color="Blue", thick=2, /noerase
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
		color="red", $
		/data, $
		'<B!IgB!N>='+string(median(bgb_galaxyzoo[eind]),format='(i6)')+cgsymbol('+-')+string(stddev(bgb_galaxyzoo[eind]),format='(i6)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		color="blue", $
		'<B!IgB!N>='+string(median(bgb_galaxyzoo[sind]),format='(i6)')+cgsymbol('+-')+string(stddev(bgb_galaxyzoo[sind]),format='(i6)'), $
		charsize=labelsize

	; Plot Bgb as fn. of redshift

	ploterror, z, bgb_galaxyzoo, bgb_galaxyzoo_err, $
		charsize=cs, $
		/nohat, $
;		yr=[-1d4,2d4], $
		xrange=[0,0.5], $
		xtitle='Redshift', $
		ytitle='B!IgB!N', $
		psym = 9

	ploterror, z[eind], bgb_galaxyzoo[eind], bgb_galaxyzoo_err[eind], $
		color="Red", $
		errcolor="red", $
		charsize=cs, $
		xrange=[0,0.5], $
		/nohat, $
;		yr=[-5d3,2d4], $
		xtitle='Redshift', $
		ytitle='B!IgB!N', $
		psym = 16

	oploterror, z[sind], bgb_galaxyzoo[sind], bgb_galaxyzoo_err[sind], $
		color="Blue", $
		/nohat, $
		psym = 9

	al_legend, /top, /left, $
		['Elliptical','Spiral'], $
		color=['Red','Blue'], $
		psym=[16,9], $
		charsize=labelsize

	if keyword_set(ps) then ps_end

	!p.multi=[0,1,1]	

	endif

	if keyword_set(stop) then stop

end


