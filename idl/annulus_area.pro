
;+
; NAME:
;       
;	ANNULUS_AREA
;
; PURPOSE:
;
;	Determine the mean number of background field galaxies 
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
;       Written by K. Willett                Nov 11
;-

pro annulus_area, cmag, nbg_per_arcmin2, quiet=quiet, noplot=noplot, stop=stop

blazardir = '~/Astronomy/Research/blazars/'
paperdir = '~/Astronomy/Research/blazars/paper/'
;restore,blazardir+'ab_structure.sav'
restore,blazardir+'sav/random_structure.sav'

zarr = fillarr(0.01,0.1,5)
angsize_1mpc = zang(1000.,zarr,/silent,/wmap7)
angsize_15mpc = zang(1500.,zarr,/silent,/wmap7)
angsize_2mpc = zang(2000.,zarr,/silent,/wmap7)
area = !dpi/3600. * (angsize_2mpc^2 - angsize_1mpc^2)

	; Changing the upper redshift limit makes no significant difference,
	; indicating that the counting magnitude is mostly the SDSS completeness limit

	bz = random
	sdss_applim = 22.2
	
	;ind = where(bz.z gt 0.43)
	ind = where(bz.z gt 0.0,nlim)

	zlim = bz[ind].z
	;arealim = bz[ind].area

	bgmagarr = [0.]
	;for i=0,n_elements(zlim)-1 do bgmagarr=[bgmagarr,*(bz[ind[i]].bg_mag)] 
	for i=0,nlim-1 do bgmagarr=[bgmagarr,*(bz[ind[i]].n_mag_r)] 
	bgmagarr=bgmagarr[1:n_elements(bgmagarr)-1]
	bgmagarr = bgmagarr[where(bgmagarr gt -100)]
	sqarcmin_to_sqdeg = (60*180.)^2 / (206265.*!dpi)^2 
	sqarcmin_to_sr = (60.)^2 / (206265.)^2 
	;totalarea = total(arealim) * sqarcmin_to_sqdeg 
	totalarea = !dpi * (3.)^2 * nlim ;* sqarcmin_to_sqdeg

	n_with_mag = where(bgmagarr lt float(cmag), nwm)

	;nbg_per_arcmin2 = n_elements(bgmagarr)/(!dpi * (3.)^2 * nlim)
	nbg_per_arcmin2 = nwm/(!dpi * (3.)^2 * nlim)

	bgbin = 0.5
	absbins = fillarr(bgbin,-25,-15)
	appbins = fillarr(bgbin,15,25)
	bins = appbins
	n = n_elements(bins)

	nbgarr = lonarr(n)
	for i=0,n-1 do begin
		tempind = where((bgmagarr gt bins[i] - bgbin/2.) and (bgmagarr lt bins[i] + bgbin/2.), bgcount2)
		nbgarr[i] = bgcount2
	endfor
	nbg = nbgarr / totalarea
	nbg_density = total(nbg)

if ~keyword_set(noplot) then begin
	!p.multi=[0,2,2]	
	cs = 1.5
	cgplot, bins, nbgarr/totalarea, psym=-16, $
		charsize=cs, $
		xtitle='m!IR!N', $
		ytitle='n!Ibg!N/mag/arcsec!E2!N'
	cgplots, [sdss_applim, sdss_applim], !y.crange, linestyle=2
	cgplot, bins, nbgarr/(totalarea*sqarcmin_to_sqdeg), psym=-16, $
		/ylog, $
		charsize=cs, $
		xtitle='m!IR!N', $
		ytitle='n!Ibg!N/mag/deg!E2!N'
	cgplots, [sdss_applim, sdss_applim], 10^!y.crange, linestyle=2
	cgplot, bins, nbgarr/(totalarea*sqarcmin_to_sr), psym=-16, $
		/ylog, $
		charsize=cs, $
		xtitle='m!IR!N', $
		ytitle='n!Ibg!N/mag/sr'
	cgplots, [sdss_applim, sdss_applim], 10^!y.crange, linestyle=2
endif

	if ~keyword_set(quiet) then print,'Mean of the total galaxies:                          ', nbg_per_arcmin2
	if ~keyword_set(quiet) then print,'Mean of the total galaxies up to counting magnitude: ', nbg_density

	if keyword_set(stop) then stop

end
