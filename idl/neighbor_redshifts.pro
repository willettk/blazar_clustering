
;+
; NAME:
;       
;	NEIGHBOR_REDSHIFTS
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
;       Written by K. Willett                
;-


; 1. Bin the sample galaxies (outside of the clustering radius) by redshift

	; delta z \simeq 0.1? 0.2? Play with this value

; 2. Plot the number of galaxies per magnitude bin as a function of observed absolute magnitude

	; bin in sizes of about half M_r

; 3. Fit a scaling constant for each redshift bin, assuming a fixed Mstar and alpha

	; simple; experiment slightly with mrstar, alpha

; 4. Save the array of scaling constants 

; 5. Run BGB.pro using the new scaling; see what the absolute value of the LF yields
;restore, blazardir+'ab_neighbors3_willettk.sav'

;;;;;;;;;;;;;;;;;;;;
; Can I replicate Figure 1 in YG87? Read it in using one of the web tools. 
;;;;;;;;;;;;;;;;;;;;

blazardir = '~/Astronomy/Research/blazars/'
paperdir = '~/Astronomy/Research/blazars/paper/'

;!p.multi=[0,2,2]
restore,blazardir+'ab_structure.sav'

;zlow = 0.162
zlow = 0.5
zhigh = 1.0

; Try four bins between the redshift limits

nrb = 4
!p.multi=[0,2,2]

for j=0,nrb-1 do begin
	
	lowbin = zlow+(zhigh-zlow)/nrb*j
	highbin =lowbin+(zhigh-zlow)/nrb

	ind = where(bz.z gt lowbin and bz.z le highbin)

	zlim = bz[ind].z
	arealim = bz[ind].area

	bgmagarr = [0.]
	for i=0,n_elements(zlim)-1 do bgmagarr=[bgmagarr,*(bz[ind[i]].bg_mag)] 
	bgmagarr=bgmagarr[1:n_elements(bgmagarr)-1]
	; Remove galaxies with no magnitude (-9999)
	bgmagarr = bgmagarr[where(bgmagarr gt -100)]
	totalarea = total(arealim)

	bgbin = 0.5
	;bins = fillarr(bgbin,floor(min(bgmagarr)),ceil(max(bgmagarr)))
	bins = fillarr(bgbin,-25,ceil(max(bgmagarr)))
	n = n_elements(bins)
	nbgarr = lonarr(n)
	for i=0,n-1 do begin
		tempind = where((bgmagarr gt bins[i] - bgbin/2.) and (bgmagarr lt bins[i] + bgbin/2.), bgcount2)
		nbgarr[i] = bgcount2
	endfor
	nbg = nbgarr / totalarea
	nbg_per_arcmin2 = total(nbg)

	nbgerr = sqrt(nbgarr)
	nbgerr[where(nbgerr eq 0)] = 1

	ploterror, bins, nbgarr, nbgerr, $
		/ylog, $
		psym=16, $
		title='<z>='+string(mean([lowbin,highbin]),format='(f4.1)'), $
		xtitle='M!Ir!I', $
		ytitle='Number per magnitude', $
		charsize=1.5, $
		;xrange=[max(bins)+1,min(bins)],$
		xrange=[-19,-25], $
		/xstyle, $
		yrange=[1d-1,1d4]

		cmaglim = 22.5 - 5*alog10(lumdist(mean([lowbin,highbin]),/wmap7,/silent)*1e6) + 5
		f0 = closeto(bins,-25)		; Should be the completeness limit of survey
		f1 = closeto(bins,cmaglim)-1

		vline, cmaglim, linestyle=1, color='black', /data, /noerase, range=[1d-1,1d4], thick=3

	parinfo = replicate({value:0.,fixed:0,limited:intarr(2),limits:fltarr(2)},3)
	parinfo[*].value = [1.,-21.73,-1.2]
	parinfo[2].fixed = 1 
	parinfo[1].limited=[1,1]
	parinfo[1].limits=[-30.,-15.]

	sf = mpfitfun('schechter_diff',bins[f0:f1],nbgarr[f0:f1],nbgerr[f0:f1],parinfo=PARINFO,/quiet)

	sarr = fillarr(0.1,-25,-19)
	cgplot, sarr, schechter_diff(sarr, sf), /over
	cgplot, sarr, schechter_diff(sarr, [sf[0],sf[1]+2,sf[2]]), /over
	cgplot, sarr, schechter_diff(sarr, [sf[0],sf[1]-2,sf[2]]), /over

	;vline, bins[f0], color='red',/data,/noerase,range=[1d-1,1d4]
	vline, bins[f1], color='red',/data,/noerase,range=[1d-1,1d4]
	vline, sf[1], color='blue',/data,/noerase,range=[1d-1,1d4]

	cgtext, -24, 1d3, $
		/data, $
		greek('phi')+'!E*!N = '+string(sf[0],format='(i4)'), $
		charsize=1.5

		print,sf


endfor

end
