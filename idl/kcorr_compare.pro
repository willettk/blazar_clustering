
;+
; NAME:
;       
;	KCORR_COMPARE
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
;       Written by K. Willett                Dec 11
;-

blazardir = '~/Astronomy/Research/blazars/'
restore, blazardir+'ab_neighbors3_willettk.sav'

match, objid_kcorr, mag_objid, ki, mi

kr = absmag_kcorr_r[ki]
distance = lumdist(mag_redshift[mi],/silent)
abs_mr = mag_r[mi] - 5 * alog10(distance*1e6) + 5

gi = where(kr gt -50 and finite(abs_mr))

zphot = mag_redshift[mi]

!p.multi=[0,2,1]

cgplot, kr[gi], abs_mr[gi], psym=2, $
	charsize=2, $
	xtitle='Sloan K-corrected M!Ir!N', $
	ytitle='Distance modulus M!Ir!N'

cgplot, /overplot, indgen(100)-50, indgen(100)-50, $
	linestyle=2, color='Red'

cgplot, zphot[gi], kr[gi] - abs_mr[gi], psym=2, $
	charsize=2, $
	xtitle='Redshift', $
	ytitle=greek('Delta')+'M!Ir!N'


end
