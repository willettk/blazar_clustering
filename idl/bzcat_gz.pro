
;+
; NAME:
;       
;	BZCAT_GZ
;
; PURPOSE:
;
;	Plot data on BL Lacs from BZCat with respect to number of GZ-classified neighbors
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
;	IDL> .r bzcat_gz
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Jul 11
;-

;;;;;;;;;;
; Plot the environment for a 3 arcminute radius around each BL Lac
;;;;;;;;;;

zs_3arcmin = '~/Astronomy/Research/blazars/bzcat_zs_3arcmin.csv'
zns_3arcmin = '~/Astronomy/Research/blazars/bzcat_zns_3arcmin.csv'

readcol, zs_3arcmin, $
	zs3_col0,zs3_objID,zs3_ra,zs3_dec,zs3_type,zs3_modelMag_u,zs3_modelMag_g,zs3_modelMag_r,zs3_modelMag_i,zs3_modelMag_z,zs3_p_el, $
	skipline=1, format='a,a,f,f,a,f,f,f,f,f,f,f', /quiet

readcol, zns_3arcmin, $
	zns3_col0,zns3_objID,zns3_ra,zns3_dec,zns3_type,zns3_modelMag_u,zns3_modelMag_g,zns3_modelMag_r,zns3_modelMag_i,zns3_modelMag_z,zns3_p_el, $
	skipline=1, format='a,a,f,f,a,f,f,f,f,f,f,f', /quiet

; Find host galaxies themselves; assumed to be any detection within 0.5 arcmin of coordinates 

zs_primary = '~/Astronomy/Research/blazars/bzcat_zs.csv'
zns_primary = '~/Astronomy/Research/blazars/bzcat_zns.csv'

readcol, zs_primary, $
	zs_col0,zs_objID,zs_ra,zs_dec,zs_type,zs_modelMag_u,zs_modelMag_g,zs_modelMag_r,zs_modelMag_i,zs_modelMag_z,zs_p_el, $
	skipline=1, format='a,a,f,f,a,f,f,f,f,f,f,f', /quiet

readcol, zns_primary, $
	zns_col0,zns_objID,zns_ra,zns_dec,zns_type,zns_modelMag_u,zns_modelMag_g,zns_modelMag_r,zns_modelMag_i,zns_modelMag_z,zns_p_el, $
	skipline=1, format='a,a,f,f,a,f,f,f,f,f,f,f', /quiet

; Find the object IDs for both samples

both_zs3objids = [zs3_objid, zns3_objid]
both_zsobjids = [zs_objid, zns_objid]

; Match
match, both_zs3objids, both_zsobjids, aa, bb

; Find the indices for objects that are not the parent BL Lac target
norepeats = setdifference(indgen(n_elements(both_zs3objids)),aa)

both_ra = [zns3_ra, zs3_ra]

norepeat_ra = both_ra[norepeats]

h_norepeats = histogram(both_ra[norepeats])
hh = where(h_norepeats gt 1, hcount)

ps_start, file='~/Astronomy/Research/blazars/bzcat_gz.ps', /color, /quiet
cghistoplot, h_norepeats, xtitle="Number of BZCat BL Lac neighbors w/in 3'", ytitle='N', thick=3, charthick=2
cgtext, 0.6, 0.6, charsize=1.5, 'N!Itot!N = '+string(n_elements(h_norepeats),format='(i4)')+' (1178)', /normal
cgtext, 0.6, 0.5, charsize=1.5, 'N!I>1!N = '+string(hcount,format='(i4)'), /normal
ps_end

end
