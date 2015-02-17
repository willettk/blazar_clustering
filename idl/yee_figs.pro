

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

pro yee_figs, stop=stop, noplot = noplot, ps=ps

	blazardir = '~/Astronomy/Research/blazars/'
	paperdir = '~/Astronomy/Research/blazars/paper/'

	restore,blazardir+'ab_structure.sav'

	; total N_gal, background N_gal, radius of field in arcsec, redshift of blazar

	ind = where(bz.z lt 0.80)
	bz=bz[ind]

	bllac = where(bz.btype eq 'BLLac')
	fsrq = where(bz.btype eq 'FSRQ')
	uncertain = where(bz.btype eq 'BLLac_candidate' or bz.btype eq 'blazar_uncertain')

	annulus_area, bz.cmag_app, ng_per_arcmin2, /quiet, /noplot

	fieldradius = 500	; kpc
	z = bz.z
	fieldsize = zang(float(fieldradius),z,/silent,/wmap7)
	nt = float(bz.n500_cmag)
	nb = !pi * (fieldsize/60.)^2 * ng_per_arcmin2
	countingmag = bz.cmag

	;nb = nt/2. * abs(randomn(seed,n_elements(nt))+1)

	bgb, nt, nb, fieldsize, z, countingmag, bgb_blazar, bgb_blazar_err

		kstwo, bgb_blazar[fsrq], bgb_blazar[bllac], d, prob
		print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)

	;;;;;;; CSV ;;;;;;

	; Print the data to a CSV file to send to Howard Yee

	writefile = '~/Astronomy/Research/blazars/yee_bgb.cat'
	openw, lun1, writefile, /get_lun
	printf, lun1, 'z        N_t         N_b        Fieldsize          m0           B_gB      B_gB_err                                                 x'
	printf, lun1, [transpose(string(z,format='(f5.3)')), transpose(string(nt,format='(i4)')), transpose(string(nb,format='(f4.2)')), transpose(string(fieldsize/60.,format='(f5.1)')), transpose(string(countingmag,format='(f5.1)')), transpose(string(bgb_blazar,format='(i6)')), transpose(string(bgb_blazar_err,format='(i6)'))]
;	printf, lun1, [transpose(z), nt[i], nb[i], fieldsize[i], countingmag[i], bgb_blazar[i], bgb_blazar_err[i]]
	close, lun1 & free_lun, lun1


	if ~keyword_set(noplot) then begin

	!p.multi=[0,2,2]	

	if keyword_set(ps) then begin
		ps_start, filename=blazardir+'yee_bgb.eps', /quiet, /color, /encap
		cs = 1.5
		labelsize=1
	endif else begin
		cs = 2
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
		'<B!IgB!N>='+string(median(bgb_blazar),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar),format='(i4)'), $
		charsize=labelsize

	; Plot 2

	cghistoplot, bgb_blazar[fsrq], $
		binsize=bgbbin, $
		yr=[0,150], $
		xtitle='B!IgB!N', $
		ytitle='N!Iblazars!N', $
		datacolor="Red",$ 
		charsize=cs

	cghistoplot,bgb_blazar[bllac], /oplot, $
		binsize=bgbbin, $
		datacolor="Blue"

	cgplots, [mean(bgb_blazar[fsrq]),mean(bgb_blazar[fsrq])], !y.crange, linestyle=2, color="Red", thick=2
	cgplots, [mean(bgb_blazar[bllac]),mean(bgb_blazar[bllac])], !y.crange, linestyle=2, color="Blue", thick=2
	;cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
	;	(!y.crange[1] - !y.crange[0]) * 0.7 + !y.crange[0], $
	;	color="red", $
	;	/data, $
	;	'<B!IgB!N>='+string(median(bgb_blazar[fsrq]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar[fsrq]),format='(i4)'), $
	;	charsize=labelsize
	;cgtext, (!x.crange[1] - !x.crange[0]) * 0.5 + !x.crange[0], $
	;	(!y.crange[1] - !y.crange[0]) * 0.6 + !y.crange[0], $
	;	/data, $
	;	color="blue", $
	;	'<B!IgB!N>='+string(median(bgb_blazar[bllac]),format='(i4)')+cgsymbol('+-')+string(stddev(bgb_blazar[bllac]),format='(i4)'), $
	;	charsize=labelsize

	al_legend, /top, /right, $
		thick=th, bthick=th, $
		['FSRQ','BL Lac'], $
		color=['Red','Blue'], $
		linestyle=[0,0], $
		charsize=labelsize

	; Plot 3

	;cgplot, z, float(nt)-float(nb), psym=9, $
	;	charsize=cs, $
	;	xtitle='z', $
	;	ytitle='Nt - Nb'

	;cgplot, z, float(nt), psym=7, /over, color='green'
	;cgplot, z, float(nb), psym=15, /over, color='pink'

	cgplot, z, bgb_blazar, $
		charsize = cs, $
		psym = 1, $
		xtitle='Blazar redshift', $
		ytitle='B!IgB!N'


	; Plot binned Bgb as fn. of redshift

	binsize = 0.05
	bins = fillarr(binsize,0.05,0.8)
	nbins = n_elements(bins)-1
	fsrq_binned = fltarr(nbins)
	fsrq_binned_err = fltarr(nbins)
	bllac_binned = fltarr(nbins)
	bllac_binned_err = fltarr(nbins)
	allblazars_binned = fltarr(nbins)
	allblazars_binned_err = fltarr(nbins)
	fsrq_bgb = bgb_blazar[fsrq]
	bllac_bgb = bgb_blazar[bllac]
	for i=0,nbins-1 do begin
		junk_fsrq = where(bz[fsrq].z gt bins[i] and bz[fsrq].z lt bins[i+1])
		fsrq_binned[i] = mean(fsrq_bgb[junk_fsrq])
		fsrq_binned_err[i] = stddev(fsrq_bgb[junk_fsrq])
		junk_bllac = where(bz[bllac].z gt bins[i] and bz[bllac].z lt bins[i+1])
		bllac_binned[i] = mean(bllac_bgb[junk_bllac])
		bllac_binned_err[i] = stddev(bllac_bgb[junk_bllac])
		junk_allblazars = where(bz.z gt bins[i] and bz.z lt bins[i+1])
		allblazars_binned[i] = mean(bgb_blazar[junk_allblazars])
		allblazars_binned_err[i] = stddev(bgb_blazar[junk_allblazars])
	endfor

	zarr = bins[0:nbins-1]+binsize/2.

	th = 2.0

	;ploterror, zarr, fsrq_binned, fsrq_binned_err, $
	ploterror, zarr, allblazars_binned, allblazars_binned_err, $
		thick=th, ythick=th, xthick=th, $
		color="Red", $
		errcolor="red", $
		charsize=cs, $
		xrange=[0,0.8], /xstyle, $
		/nohat, $
;		yr=[-5d3,2d4], $
		xtitle='Blazar redshift', $
		ytitle='B!IgB!N (binned)', $
		psym = 16

;	oploterror, zarr, bllac_binned, bllac_binned_err, $
;		thick=th, $
;		color="Blue", $
;		/nohat, $
;		psym = 9

;	al_legend, /top, /right, $
;		thick=th, bthick=th, $
;		['FSRQ','BL Lac'], $
;		color=['Red','Blue'], $
;		psym=[16,9], $
;		charsize=labelsize

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

	;al_legend, /bottom, /left, $
	;	['FSRQs','BL Lac'], $
	;	color=['Red','Blue'], $
	;	psym=[16,9], $
	;	charsize=labelsize

	if keyword_set(ps) then ps_end

	!p.multi=[0,1,1]	


	endif

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

	;!p.multi=[0,3,2]
	;mpc2cm = 3.086d24
	;mjy2cgs = 1d-29
	;xray2cgs = 1d-12

	;ind = where(bz.spindex_ro gt -9.9 and bz.spindex_ro ne 0.)
	;cgplot, bz[ind].spindex_ro, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+'!Iro!N', ytitle='B!IgB!N'

	;ind = where(bz.spindex_ox gt -9.9 and bz.spindex_ox ne 0.)
	;cgplot, bz[ind].spindex_ox, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+'!Iox!N', ytitle='B!IgB!N'

	;ind = where(bz.spindex_rx gt -9.9 and bz.spindex_rx ne 0.)
	;cgplot, bz[ind].spindex_rx, bgb_blazar[ind], psym=9, xtitle=greek('alpha')+'!Irx!N', ytitle='B!IgB!N'

	;ind = where(bz.flux_radio gt 0.)
	;cgplot, bz[ind].flux_radio*mjy2cgs*4*!dpi*(lumdist(bz[ind].z,/silent,/wmap7)*mpc2cm)^2*1420.*1e6, bgb_blazar[ind], psym=9, xtitle='1.4 GHz '+greek('nu')+'L!I'+greek('nu')+'!N [erg/s]', ytitle='B!IgB!N',/xlog

	;ind = where(bz.flux_xray gt 0.)
	;cgplot, bz[ind].flux_xray*xray2cgs*4*!dpi*(lumdist(bz[ind].z,/silent,/wmap7)*mpc2cm)^2, bgb_blazar[ind], psym=9, xtitle='0.1-2.4 keV '+greek('nu')+'L!I'+greek('nu')+'!N [erg/s]', ytitle='B!IgB!N',/xlog

	if keyword_set(stop) then stop

end
