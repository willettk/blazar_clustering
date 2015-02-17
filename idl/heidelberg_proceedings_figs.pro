
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

;pro heidelberg_proceedings_figs, stop=stop, noplot = noplot, ps=ps, readfiles=readfiles, multiwave=multiwave, $
	;neighbors = neighbors, lowz=lowz, highz=highz, cosmo=cosmo

	blazarsavdir = '~/Astronomy/Research/blazars/sav/'
	fitsdir = '~/Astronomy/Research/blazars/fits/'
	figdir = '~/Astronomy/meetings/heidelberg2012/aipproc_6s/'

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
		;print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)

	!p.multi=[0,2,1]	

	ps_start, filename=figdir+'heidelberg_bgb.eps', /quiet, /color, /encap, ysize=6, xsize=12
	cs = 1.7
	labelsize=1.3
	th=5
	
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
		endif ;else print,'No FSRQs between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_fsrq = where(bz[fsrq].z gt bins[i] and bz[fsrq].z lt bins[i+1],jf)
		if jf gt 0 then begin
			fsrq_binned[i] = wmean(fsrq_bgb[junk_fsrq],fsrq_bgb_err[junk_fsrq],error=wmean_err_fsrq,/nan)
			fsrq_binned_err[i] = stddev(fsrq_bgb[junk_fsrq])
			fsrq_binned_err[i] = wmean_err_fsrq
		endif ;else print,'No FSRQs between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
		junk_bllac = where(bz[bllac].z gt bins[i] and bz[bllac].z lt bins[i+1],jb)
		if jb gt 0 then begin
			bllac_binned[i] = wmean(bllac_bgb[junk_bllac],bllac_bgb_err[junk_bllac],error=wmean_err_bllac,/nan)
			bllac_binned_err[i] = stddev(bllac_bgb[junk_bllac])
			bllac_binned_err[i] = wmean_err_bllac
		endif ;else print,'No BL Lacs between '+string(bins[i],format='(f4.2)')+' and '+string(bins[i+1],format='(f4.2)')
	endfor

	zarr = bins[0:nbins-2]+binsize/2.

	cgplot, indgen(10), $
		/nodata, $
		position = [0.10,0.15,0.45,0.90], $
		thick=th, ythick=th, xthick=th, $
		charsize=cs, $
		xrange=[0,round(highz*10.)/10.], /xstyle, $
;		yrange=[min(bgb_blazar),max(bgb_blazar)], /ystyle, $
		yr=[-1000,1000], $
		xtitle='Blazar redshift', $
		ytitle='B!IgB!N [Mpc!E1.77!N]'

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

	bb = where(bllac_binned ne 0.)
	oploterror, zarr[bb], bllac_binned[bb], bllac_binned_err[bb], $
		thick=th, $
		color="Blue", $
		/nohat, $
		psym = -16

	fb = where(fsrq_binned ne 0.)
	oploterror, zarr[fb], fsrq_binned[fb], fsrq_binned_err[fb], $
		thick=th, $
		color="red", $
		/nohat, $
		psym = -6

	al_legend, /top, /right, $
		bthick=th, $
		['BL Lac','FSRQ'], $
		color=['Blue','Red'], $
		psym=[16,6], $
		charsize=labelsize

paperdir='~/Astronomy/Research/blazars/paper/'
b = mrdfits('~/Astronomy/Research/blazars/fits/meyer_bz_allseds.fits',1,/silent)

n = n_elements(b)

	if n_elements(lowz) eq 0 then lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	if n_elements(highz) eq 0 then highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	if n_elements(neighbors) eq 0 then neighbors = 10

	ztemp = fltarr(n)
	for i = 0,n-1 do begin
		znum = strnumber(b[i].used_redshift,val)
		if znum eq 1 then ztemp[i] = val else ztemp[i] = b[i].z
	endfor

;	ind = where(ztemp lt highz and ztemp gt lowz and b.n500_cmag ge neighbors and b.bg gt 0 and (b.sed_code eq 'Uex' or b.sed_code eq 'Tex'))
	ind = where(ztemp gt lowz and ztemp lt highz and (b.sed_code eq 'Uex' or b.sed_code eq 'Tex'),nb)
	b = b[ind]

	b.btype = strtrim(b.btype,2)
	bllac = where(b.btype eq 'BLLac' or b.btype eq 'Plotkin_blazar' or b.btype eq 'HBL' or b.btype eq 'lBLLac')
	fsrq = where(b.btype eq 'FSRQ')
	uncertain = where(b.btype eq 'BLLac_candidate' or b.btype eq 'blazar_uncertain' or b.btype eq 'Blazar')

	z = ztemp[ind]
	fieldsize = b.fieldsize
	nt = float(b.n500_cmag)
	nb = float(b.bg)
	countingmag = b.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bmeyer, bmeyer_err

bmin = min(bmeyer)
bmax = max(bmeyer)
	bmin = -500.
	bmax = 1500.

	cs = 1.7
	labelsize=1.3
	th = 5
	hsize=600
	
cgplot, b.nupeak, b.lpeak, $
	background='white', $
	position=[0.58,0.15,0.95,0.90], $
	charsize=cs, $
	thick=th, ythick=th, xthick=th, $
	/nodata, $
	yr=[43.95,47], /ystyle, $
	xr=[12,17], /xstyle, $
;	title='Updated Meyer plot', $
	xtitle='log ('+greek('nu')+'!Ipeak!N) [Hz]', $
	ytitle='log ('+greek('nu')+'L!I'+greek('nu')+'!N) [erg s!E-1!N]'

cgloadct, 13
minsize = 0.5
sizescale = 5.

for i=0,n_elements(bmeyer) - 1 do begin
	case strtrim(b[i].btype,2) of
		'BLLac': sym = 16
		'Plotkin_blazar': sym = 16
		'FSRQ': sym = 6
		else: sym=34
	endcase
	if sym ne 34 then begin
		cgplot, b[i].nupeak, b[i].lpeak, $
			/over, $
			symsize=(bmeyer[i] - bmin) / (bmax-bmin) * sizescale + minsize, $
			color=fix((z[i] - min(z)) / (max(z)-min(z)) * 255.), $
			psym=sym
	endif
endfor


al_legend, /top,/right, psym=[16,6], ['BL Lac','FSRQ'], charsize=labelsize, symsize=2, outline_color='black', textcolor='black', colors='black'
cgcolorbar, position=[0.60,0.95,0.90,0.99], range=[min(z),max(z)], /top, color='black'
cgtext, 0.92, 0.96, 'z', charsize=1.5, /normal, color='black'

s1 = (-500 - bmin) / (bmax-bmin) * sizescale + minsize
s2 = (0 - bmin) / (bmax-bmin) * sizescale + minsize
s3 = (500 - bmin) / (bmax-bmin) * sizescale + minsize
s4 = (1000 - bmin) / (bmax-bmin) * sizescale + minsize
;al_legend, position=[0.87,0.75], psym=16, ['-500','0','500','1000'], charsize=1.0, symsize=[s1,s2,s3,s4],/normal,spacing = 3, outline_color='black', textcolor='black', colors='black'
;cgtext, /normal, 0.91, 0.37, 'B!IgB!N',charsize=cs

; Overplot tracks from Meyer et al. (2011)

readcol, '/Applications/Dexter/meyer_sequence.gif.tracka', xa, ya, format='f,f', /skipline, /silent
readcol, '/Applications/Dexter/meyer_sequence.gif.trackb', xb, yb, format='f,f', /skipline, /silent

slope = (ya[1] - ya[0]) / (xa[1] - xa[0])
intercept = ya[0] - slope * xa[0]

;cgplots, xa, ya, color='black', linestyle=2
;cgplots, xb, yb, color='black', linestyle=2

; Kluges to make sure lines don't run off the edge of my truncated plot

;cgarrow, (!Y.crange[1] - intercept) / slope, !y.crange[1],(!Y.crange[0] - intercept) / slope, !y.crange[0], color='black', linestyle=2, thick=th,/data,/solid,hsize=hsize
;cgplots, [15.80,16.47,17.00], [44.24,44.50,44.82], color='black', linestyle=2, thick=th
;cgarrow, 15.80,44.24, 15.05,44.00,color='black', linestyle=2, thick=th, /data, hsize=hsize,/solid

;cgtext, 13.7, 46.7, /data, 'Single-component jet', charsize=labelsize
;cgtext, 15.85, 44.2, /data, 'Decelerating jets', charsize=labelsize


	ps_end

	; Make a FITS file with all blazar data

btemp = {$
   BNAME           :    'BZUJ0725-0054', $
   BTYPE           :    'Blazar', $
   RA              :          111.461, $
   DEC             :        -0.915560, $
   Z               :         0.128000, $
   SEARCHID        :               7, $
   BG              :         153, $
   BGB             :         1L, $
   BGB_ERR         :         1L, $
   FIELDSIZE       :          219.324, $
   CMAG            :         -19.1837, $
   CMAG_APP        :          19.7009, $
   N500            :          14, $
   N500_CMAG       :          10, $
   FLUX_RADIO      :          1399.00, $
   FLUX_XRAY       :          3.14000, $
   SP_MAG_R        :          16.6000, $
   SPINDEX_RO      :         -9.90000, $
   SPINDEX_RX      :         -9.90000, $
   SPINDEX_OX      :          1.37300 }

bz_all = replicate(btemp,n_elements(bgb_blazar))
	
bz_all.BNAME         = bz.BNAME       
bz_all.BTYPE         = bz.BTYPE       
bz_all.RA            = bz.RA          
bz_all.DEC           = bz.DEC         
bz_all.Z             = bz.Z           
bz_all.SEARCHID      = bz.SEARCHID    
bz_all.BG            = bz.BG          
bz_all.BGB           = bgb_blazar
bz_all.BGB_ERR       = bgb_blazar_err
bz_all.FIELDSIZE     = bz.FIELDSIZE   
bz_all.CMAG          = bz.CMAG        
bz_all.CMAG_APP      = bz.CMAG_APP    
bz_all.N500          = bz.N500        
bz_all.N500_CMAG     = bz.N500_CMAG   
bz_all.FLUX_RADIO    = bz.FLUX_RADIO  
bz_all.FLUX_XRAY     = bz.FLUX_XRAY   
bz_all.SP_MAG_R      = bz.SP_MAG_R    
bz_all.SPINDEX_RO    = bz.SPINDEX_RO  
bz_all.SPINDEX_RX    = bz.SPINDEX_RX  
bz_all.SPINDEX_OX    = bz.SPINDEX_OX  

mwrfits, bz_all, fitsdir+'bz_all.fits'

end
