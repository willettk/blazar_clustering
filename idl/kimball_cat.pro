
;+
; NAME:
;       
;	KIMBALL_CAT
;
; PURPOSE:
;
;	Find SDSS-detected radio galaxies with RA, dec and classified by radio morphology
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

all4radio = '~/Astronomy/Research/blazars/kimball_all4radio_SDSS.csv'
first_nvss = '~/Astronomy/Research/blazars/kimball_FIRST_NVSS.csv'

readcol, all4radio, $
	/silent, $
	skipline=1, $
	UNIQ_ID,$
	RA,DEC,$
	FIRST_PEAK_FLUX,FIRST_FLUX,$
	NVSS_FLUX,WENSS_FLUX,GB6_FLUX, $
	NEAR_TYPE,NEAR_FLAGS,$
	NEAR_MODELMAG_U,NEAR_MODELMAG_G,NEAR_MODELMAG_R,NEAR_MODELMAG_I,NEAR_MODELMAG_Z,$
	NEAR_MODELMAGERR_U,NEAR_MODELMAGERR_G,NEAR_MODELMAGERR_R,NEAR_MODELMAGERR_I,NEAR_MODELMAGERR_Z,$
	DISTANCE,WENSS_DISTANCE,GB6_DISTANCE,NEAR_DISTANCE, $
	format=strjoin(replicate('a, ',24))+'a'

ngals = n_elements(uniq_id)

;readcol, first_nvss, $
;	/silent, $
;	skipline=1, $
;	FN_UNIQ_ID,$
;	FN_RA,FN_DEC,$
;	FN_FIRST_PEAK_FLUX,FN_FIRST_FLUX,$
;	FN_NVSS_FLUX,FN_GB6_FLUX,FN_WENSS_FLUX,$
;	FN_DISTANCE,FN_WENSS_DISTANCE,FN_GB6_DISTANCE, $
;	format=strjoin(replicate('a, ',10))+'a'
;
;match, uniq_id, fn_uniq_id, sdssind, fnind

; Compute the radio flux "AB magnitude"

first_flux = float(first_flux)
first_peak_flux = float(first_peak_flux)
nvss_flux = float(nvss_flux)

tfirst = -2.5 * alog10(first_flux * 1d-3 /3631.)
tnvss = -2.5 * alog10(nvss_flux * 1d-3 /3631.)

deltat = tfirst - tnvss

theta = sqrt(first_flux / first_peak_flux)

!p.multi=[0,1,2]

cghistoplot, deltat, xtitle='t!IFIRST!N - t!INVSS!N', xr=[-1,2]
cgplots, [0.35,0.35], !y.crange, linestyle=2, /data, thick=2
cgtext, 1.0, 1500, 'Complex'
cgtext, -0.8, 1500, 'Simple',charsize=2.0
cgtext, -0.8, 1500, '(Compact + Resolved)',charsize=1.0

cghistoplot, alog10(theta^2), xtitle='log '+greek('theta')+'!E2!N', xr=[-0.20,0.50]
cgplots, [0.05,0.05], !y.crange, linestyle=2, /data, thick=2
cgtext, 0.2, 1500, 'Resolved'
cgtext, -0.15, 1500, 'Compact'

ind_complex = where(deltat gt 0.35, ncomplex)
ind_compact = where(deltat lt 0.35 and alog10(theta^2) lt 0.05, ncompact)
ind_resolved = where(deltat lt 0.35 and alog10(theta^2) gt 0.05, nresolved)

print,'Number of complex  sources: ', ncomplex
print,'Number of compact  sources: ', ncompact
print,'Number of resolved sources: ', nresolved

; Write a CSV file similar to the blazar one on which SDSS can search for photometry and nearby neighbors

kimballid = 'KI_'+uniq_id
kimball_rmorph = strarr(ngals)
kimball_rmorph[ind_complex] = 'COMPLEX'
kimball_rmorph[ind_compact] = 'COMPACT'
kimball_rmorph[ind_resolved] = 'RESOLVED'
kimballarr = [transpose(kimballid), transpose(string(ra)), transpose(string(dec)), transpose(kimball_rmorph)]

openw, lun1, '~/Astronomy/Research/blazars/kimball_radiogals.cat', /get_lun
printf, lun1, 'kimball_id              ra              dec               radio_morph'
printf, lun1, kimballarr
close, lun1
free_lun, lun1

end
