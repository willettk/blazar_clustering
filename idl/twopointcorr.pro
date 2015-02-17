
;+
; NAME:
;       
;	TWOPOINTCORR
;
; PURPOSE:
;
;	Compute the two-point correlation function using the Landy & Szalay (1993) estimator
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
;	Plotting the one square degree field (33.85 arcmin radius) takes 9 minutes (551 sec).	
;
; REVISION HISTORY
;       Written by K. Willett                Oct 11
;-

; Load data

	; Choose a field of a given angular size (corresponds to 10 Mpc at some assumed redshift), download information for all galaxies therein
	
	; z=0.05 is 1 Mpc/arcsec
	; z=0.02 is 400 kpc/arcsec

; Try circle with radius of 3 arcmin; download all matching galaxies from CasJobs

dir = '~/Astronomy/Research/blazars/'

filelist = [$
'onearcmin_willettk.csv', $
'threearcmin_willettk.csv', $
'fivearcmin_willettk.csv', $
'tenarcmin_willettk.csv', $
'twentyarcmin_willettk.csv', $
'twentyfivearcmin_willettk.csv', $
'onesqdeg_willettk.csv']

nf = n_elements(filelist)
timearr = fltarr(nf)

nbins = 50
binsarr = dblarr(nbins-1,nf)
xiarr = dblarr(nbins-1,nf)

search_radii = [1., 3., 5., 10., 20., 25., 33.85]

for k=0, nf - 1 do begin

	time_start = systime(1)

	ra_cen = 185.0
	dec_cen = 40.0

	readcol, dir+filelist[k], $
		objID,ra,dec,type,r,distance, $
		format='a,f,f,i,f,f', $
		skipline=1, $
		/silent

	xarr = fillarr(1d-5,184d,186d)

	; Rough approximation of circular area searched on the sky

	disarr = fillarr(1d-5,0.04,0.06)
	darr = disarr*0.
	for i=0,n_elements(disarr) - 1 do begin
		gcirc, 2, ra_cen, dec_cen, ra_cen, dec_cen + disarr[i], dis
		darr[i] = dis
	endfor
	dismin = closeto(darr,search_radii[k])

	;cgplot, xarr, sqrt(disarr[dismin]^2 - (xarr - ra_cen)^2) + dec_cen, $
	;	/over,$
	;	color="blue"

	;cgplot, xarr, -1. * sqrt(disarr[dismin]^2 - (xarr - ra_cen)^2) + dec_cen, $
	;	/over,$
	;	color="blue"

	; Make believe I have five data points, arranged at a quincunx around the center at some specified angular offset

	angoffset = 50.
	offsetmag_dec = 0.05
	offsetmag_ra = 0.15
	offsetmag_dec = 1.0
	offsetmag_ra = 1.0

	disarr_n = fillarr(1d-4,0.00,offsetmag_dec)
	disarr_s = fillarr(1d-4,-1 * offsetmag_dec,0.00)
	darr = disarr_n*0.
	for i=0,n_elements(disarr_n) - 1 do begin
		gcirc, 2, ra_cen, dec_cen, ra_cen, dec_cen + disarr_n[i], dis
		darr[i] = dis
	endfor
	dec_n = dec_cen+disarr_n[closeto(darr,angoffset)]

	darr = disarr_s*0.
	for i=0,n_elements(disarr_s) - 1 do begin
		gcirc, 2, ra_cen, dec_cen, ra_cen, dec_cen + disarr_s[i], dis
		darr[i] = dis
	endfor
	dec_s = dec_cen+disarr_s[closeto(darr,angoffset)]

	disarr_e = fillarr(1d-4,0.00,offsetmag_ra)
	disarr_w = fillarr(1d-4,-1 * offsetmag_ra,0.00)
	darr = disarr_e*0.
	for i=0,n_elements(disarr_e) - 1 do begin
		gcirc, 2, ra_cen, dec_cen, ra_cen+disarr_e[i], dec_cen, dis
		darr[i] = dis
	endfor
	ra_e = ra_cen+disarr_e[closeto(darr,angoffset)]

	darr = disarr_w*0.
	for i=0,n_elements(disarr_w) - 1 do begin
		gcirc, 2, ra_cen, dec_cen, ra_cen+disarr_w[i], dec_cen, dis
		darr[i] = dis
	endfor
	ra_w = ra_cen+disarr_w[closeto(darr,angoffset)]

	ra_data = [ra_cen, ra_cen, ra_cen, ra_w, ra_e]
	dec_data= [dec_cen, dec_n, dec_s, dec_cen, dec_cen]

;	cgplot, /over, color="Green", ra_data, dec_data,psym=2

	; Measure distances between all pairs

	; Total number of unique pairs is N(N-1) / 2

	; DD pairs

	nd = long(n_elements(ra_data))
	DDarr = fltarr(nd*(nd-1)/2)
	i=0
	j=1
	DDind = 0L
	for i = 0L, nd-1L do begin
		j=i+1
		while j lt nd and j gt i do begin
			gcirc, 2, ra_data[i], dec_data[i], ra_data[j], dec_data[j], tempdist
			DDarr[DDind] = tempdist
			j+=1
			DDind+=1
		endwhile
	endfor

	; RR pairs

	nr = long(n_elements(ra))
	RRarr = fltarr(nr*(nr-1)/2)
	i=0
	j=1
	RRind = 0L
	for i = 0, nr-1L do begin
		j=i+1
		while j lt nr and j gt i do begin
			gcirc, 2, ra[i], dec[i], ra[j], dec[j], tempdist
			RRarr[RRind] = tempdist
			j+=1
			RRind+=1
		endwhile
	endfor

	; DR pairs

	DRarr = fltarr(nr*nd)
	i=0
	j=1
	DRind = 0L
	for i = 0L, nd-1L do begin
		for j=0L, nr-1L do begin
			gcirc, 2, ra_data[i], dec_data[i], ra[j], dec[j], tempdist
			DRarr[DRind] = tempdist
			DRind+=1
		endfor
	endfor

	; Compute the two-point correlation function
	
	; Bin the radii
	
	angbins = lindgen(nbins)/(nbins-1.)* max([DDarr,DRarr,RRarr])
	
	binarr_dd = intarr(nbins-1)
	binarr_dr = binarr_dd
	binarr_rr = binarr_dd
	
	for i=0,nbins-2 do begin
		tempbin = where(DDarr gt angbins[i] and DDarr le angbins[i+1], tempcount)
		binarr_dd[i] = tempcount
		tempbin = where(DRarr gt angbins[i] and DRarr le angbins[i+1], tempcount)
		binarr_dr[i] = tempcount
		tempbin = where(RRarr gt angbins[i] and RRarr le angbins[i+1], tempcount)
		binarr_rr[i] = tempcount
	endfor
	
	xi = (binarr_dd * float(nr^2) + binarr_rr * float(nd^2) - 2. * nr * nd * binarr_dr) / (float(nd^2) * binarr_rr)
	binsplot = (angbins[1:nbins-1] - angbins[0:nbins-2])/2. + angbins[0:nbins-2]

	binsarr[*,k] = binsplot
	xiarr[*,k] = xi

	time_end = systime(1)
	timearr[k] = time_end - time_start
	print, 'Elapsed time for '+string(search_radii[k],format='(f6.1)')+' arcmin search: ', string(time_end - time_start,format='(f6.2)')+' sec'

endfor

;save, timearr, search_radii, binsarr, xiarr, filename=dir+'twopointcorr.sav'

; Plot the sky distribution of galaxies

!p.multi=[0,2,2]
cs=1.5

cgplot, ra, dec, $
	charsize=cs, $
;	/isotropic, $
	/xstyle, /ystyle, $
	psym=3 

; Plot the two-point correlation function
 
 	; Big peak in xi at the value of the search radius -- likely caused by edge effects. Maybe use a subset of galaxies within a larger grid?

ploterror, binsplot, xi, sqrt(abs(xi)), $
	charsize = cs, $
	xtitle='Angular separation [arcsec]', $
	ytitle=greek('xi')+'['+greek('theta')+']', $
	psym=15

; Plot the time elapsed as function of search radius

cgplot, !dpi * (search_radii)^2, timearr, $
	/ylog, $
	charsize=cs, $
	psym = 9, $
	xtitle='Search area [sq. arcmin]', $
	ytitle='Time elapsed [sec]'

!p.multi=[0,1,1]

stop

end
