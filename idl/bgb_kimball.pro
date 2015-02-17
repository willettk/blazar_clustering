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
;	Big question: what could I change in BGB to suppress the turnup at low z?
;
;		ngal, nbackground, LF, field size, gamma index?
;
; REVISION HISTORY
;       Written by K. Willett                
;-

pro bgb_kimball, stop=stop, noplot = noplot, ps=ps

	blazardir = '~/Astronomy/Research/blazars/sav/'
	paperdir = '~/Astronomy/Research/blazars/paper/'

	restore,blazardir+'kimball_structure.sav'

	; Plot the redshift distribution of the galaxies

	ps_start, filename=paperdir+'kimball_zhist.eps',/encap,/color,/quiet
	cghistoplot, kimball.z, xtitle='Radio galaxy redshift', charsize=2
	ps_end

	; total N_gal, background N_gal, radius of field in arcsec, redshift of radio galaxy

	ind = where(kimball.z lt 0.75 and kimball.z gt 0.4)
	kimball=kimball[ind]

	cx = where(kimball.radio_morph eq 'complex')
	cp = where(kimball.radio_morph eq 'compact')
	rs = where(kimball.radio_morph eq 'resolved')

	annulus_area, kimball.cmag_app, ng_per_arcmin2, /quiet, /noplot

	fieldradius = 500	; kpc
	z = kimball.z
	fieldsize = zang(float(fieldradius),z,/silent,/wmap7)
	nt = float(kimball.n500_cmag)
	nb = !pi * (fieldsize/60.)^2 * ng_per_arcmin2
	countingmag = kimball.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bgb_kimball, bgb_kimball_err

		kstwo, bgb_kimball[cp], bgb_kimball[cx], d_cp_cx, prob_cp_cx
		kstwo, bgb_kimball[cp], bgb_kimball[rs], d_cp_rs, prob_cp_rs
		kstwo, bgb_kimball[rs], bgb_kimball[cx], d_rs_cx, prob_rs_cx
		print,''
		print,'Compact-Complex KS result:',d_cp_cx,prob_cp_cx,sqrt(2d) * inverf(1d - prob_cp_cx)
		print,'Compact-Resolved KS result:',d_cp_rs,prob_cp_rs,sqrt(2d) * inverf(1d - prob_cp_rs)
		print,'Resolved-Complex KS result:',d_rs_cx,prob_rs_cx,sqrt(2d) * inverf(1d - prob_rs_cx)
		print,''

	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=paperdir+'bgb_kimball.eps', /quiet, /color, /encap
		cs = 1.5
		labelsize=1
	endif else begin
		cs = 2
		labelsize = 1.5
	endelse
	
	; Plot the distribution of B_gB

	bgbbin = 200
	cghistoplot, bgb_kimball, $
		;xr=[-500,3000], $
		datacolor='black', $
		binsize=bgbbin, $
		xtitle='B!IgB!N', $
		ytitle='N!Igalaxies!N', $
		charsize=cs

	cgplots, [mean(bgb_kimball),mean(bgb_kimball)], !y.crange, linestyle=2, thick=2, color='Black'
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		'<B!IgB!N>='+string(median(bgb_kimball),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_kimball),format='(i4)'), $
		charsize=labelsize

	cghistoplot, bgb_kimball[rs], $
		xr=[-500,3000], $
		binsize=bgbbin, $
		yr=[0,150], $
		xtitle='B!IgB!N', $
		ytitle='N!Igalaxies!N', $
		datacolor="Red",$ 
		charsize=cs

	cghistoplot,bgb_kimball[cx], /oplot, $
		binsize=bgbbin, $
		datacolor="Blue"

	cghistoplot,bgb_kimball[cp], /oplot, $
		binsize=bgbbin, $
		datacolor="forest green"

	cgplots, [mean(bgb_kimball[rs]),mean(bgb_kimball[rs])], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, [mean(bgb_kimball[cx]),mean(bgb_kimball[cx])], !y.crange, linestyle=2, color="Blue", thick=2
	cgplots, [mean(bgb_kimball[cp]),mean(bgb_kimball[cp])], !y.crange, linestyle=2, color="forest green", thick=2
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
		color="red", $
		/data, $
		'<B!IgB!N>='+string(median(bgb_kimball[rs]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_kimball[rs]),format='(i4)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
		/data, $
		color="blue", $
		'<B!IgB!N>='+string(median(bgb_kimball[cx]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_kimball[cx]),format='(i4)'), $
		charsize=labelsize
	cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
		(!y.crange[1] - !y.crange[0]) * 0.5 + !y.crange[0], $
		/data, $
		color="forest green", $
		'<B!IgB!N>='+string(median(bgb_kimball[cp]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_kimball[cp]),format='(i4)'), $
		charsize=labelsize

	cgplot, z, float(nt)-float(nb), psym=9, $
		charsize=cs, $
		xtitle='z', $
		ytitle='Nt - Nb'

	; Plot binned Bgb as fn. of redshift

	binsize = 0.05
	bins = fillarr(binsize,0.4,0.8)
	nb = n_elements(bins)
	compact_binned = fltarr(nb)
	compact_binned_err = fltarr(nb)
	complex_binned = fltarr(nb)
	complex_binned_err = fltarr(nb)
	resolved_binned = fltarr(nb)
	resolved_binned_err = fltarr(nb)
	complex_bgb = bgb_kimball[cx]
	compact_bgb = bgb_kimball[cp]
	resolved_bgb = bgb_kimball[rs]
	for i=0,nb-2 do begin
		junk_resolved = where(kimball[rs].z gt bins[i] and kimball[rs].z lt bins[i+1])
		resolved_binned[i] = mean(resolved_bgb[junk_resolved])
		resolved_binned_err[i] = stddev(resolved_bgb[junk_resolved])
		junk_complex = where(kimball[cx].z gt bins[i] and kimball[cx].z lt bins[i+1])
		complex_binned[i] = mean(complex_bgb[junk_complex])
		complex_binned_err[i] = stddev(complex_bgb[junk_complex])
		junk_compact = where(kimball[cp].z gt bins[i] and kimball[cp].z lt bins[i+1])
		compact_binned[i] = mean(compact_bgb[junk_compact])
		compact_binned_err[i] = stddev(compact_bgb[junk_compact])
	endfor

	zarr = bins[0:nb-2]+binsize/2.

	ploterror, zarr, resolved_binned, resolved_binned_err, $
		thick=th, ythick=th, xthick=th, $
		color="Red", $
		errcolor="red", $
		charsize=cs, $
		xrange=[0,1.0], /xstyle, $
		yr=[0,1500],$
		/nohat, $
;		yr=[-5d3,2d4], $
		xtitle='Radio galaxy redshift', $
		ytitle='B!IgB!N', $
		psym = 16

	oploterror, zarr, compact_binned, compact_binned_err, $
		thick=th, $
		color="Blue", $
		/nohat, $
		psym = 9

	oploterror, zarr, complex_binned, complex_binned_err, $
		thick=th, $
		color="forest green", $
		/nohat, $
		psym = 15

	al_legend, /top, /right, $
		thick=th, bthick=th, $
		['Resolved','Compact','Complex'], $
		color=['Red','Blue','forest green'], $
		psym=[16,9,15], $
		charsize=labelsize

	;ploterror, z[fsrq], bgb_kimball[fsrq], bgb_kimball_err[fsrq], $
	;	color="Red", $
	;	errcolor="red", $
	;	charsize=cs, $
	;	xrange=[0,0.8], /xstyle, $
	;	/nohat, $
	;	xtitle='Redshift', $
	;	ytitle='B!IgB!N', $
	;	psym = 16

	;oploterror, z[cx], bgb_kimball[cx], bgb_kimball_err[cx], $
	;	color="Blue", $
	;	/nohat, $
	;	psym = 9

	;al_legend, /bottom, /left, $
	;	['FSRQs','BL Lac'], $
	;	color=['Red','Blue'], $
	;	psym=[16,9], $
	;	charsize=labelsize

	if keyword_set(ps) then ps_end


	endif

	; Compute how many of the neighbors might be associated with the radio galaxy

	; Photo-z

	;nb = n_elements(kimball)
	;pzarr = lonarr(nb)
	;meanz = fltarr(nb)
	;meanzerr = fltarr(nb)
	;mylim = 0.2
	;for i=0,nb-1 do begin
	;	zgal = kimball[i].z
	;	tempz = *(kimball[i].n_redshift)
	;	tempzerr = *(kimball[i].n_redshift_err)
	;	tempzerr[where(tempzerr lt 0)] = 0.
	;	meanz[i] = mean(tempz)
	;	meanzerr[i] = mean(tempzerr[where(tempzerr gt 0)])
	;	junk = where(abs(zgal - tempz) lt (tempzerr > mylim),tcount)
	;	pzarr[i] = tcount
	;endfor

	;!p.multi=[0,2,2]	
	;cghistoplot, pzarr, xtitle='N of possible spectral neighbors per radio galaxy',charsize=1
	;cghistoplot, kimball.n500, xtitle='N of neighbors per radio galaxy',charsize=1
	;cghistoplot, kimball.n500_cmag, xtitle='N of neighbors per radio galaxy (counting mag)',charsize=1
	;cgplot, meanz, meanzerr, psym=16, xtitle='Mean redshift', ytitle='Mean redshift err',charsize=1

	!p.multi=[0,1,1]	

	if keyword_set(stop) then stop

end
