;+
; NAME:
;       
;	spiral_STRMAKE
;
; PURPOSE:
;
;	Create structure with SDSS information on neighbors for GZ-identified spiral galaxies
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
;	IDL> spiral_strmake
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Dec 11
;-

pro spiral_strmake, stop=stop, readbfiles=readbfiles, quiet=quiet

	heap_gc		; Clean up heap variables

	if n_elements(phi_star) eq 0 then phi_star = 0.5

	csvdir = '~/Astronomy/Research/blazars/csv/'
	savdir = '~/Astronomy/Research/blazars/sav/'

	spiralfile = csvdir+'spiral_willettk.csv'	; 3 arcmin radius
	neighborsfile = csvdir+'spiral_neighbors_willettk.csv'	; 3 arcmin radius
	sdssfile = csvdir+'spiral_sdss_willettk.csv'

	if keyword_set(readbfiles) then begin

		readcol, spiralfile, $
			objID,$
			ra,dec,$
			redshift,redshift_err,$
			p_el_debiased,p_cs_debiased,$
			search_id, $
			format='a,f,f,f,f,f,f,i', $
			skipline=1

		readcol, neighborsfile, $
			n_objID,$
			n_ra,n_dec,$
			n_redshift,n_redshift_err,$
			n_p_el_debiased,n_p_cs_debiased,$
			n_search_id,n_matched_id, $
			format='a,f,f,f,f,f,f,i,a', $
			skipline=1

		readcol, sdssfile, $
			n_mag_u,n_mag_g,n_mag_r,n_mag_i,n_mag_z,$
			n_mag_err_u,n_mag_err_g,n_mag_err_r,n_mag_err_i,n_mag_err_z,$
			sdss_ra, sdss_dec, $
			n_type,$
			sdss_redshift,sdss_redshift_err,$
			n_absmag_U,n_absmag_G,n_absmag_R,n_absmag_I,n_absmag_Z,$
			n_kcorrU,n_kcorrG,n_kcorrR,n_kcorrI,n_kcorrZ,$
			sdss_search_id,sdss_objid, $
			format='f,f,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f,f,f,f,f,f,f,f,f,i,a', $
			skipline=1

		save, filename=savdir+'spiral_neighbors_willettk.sav', $
			objID,ra,dec,redshift,redshift_err,p_el_debiased,p_cs_debiased,search_id, $
			n_objID,$
			n_ra,n_dec,$
			n_redshift,n_redshift_err,$
			n_p_el_debiased,n_p_cs_debiased,$
			n_search_id,n_matched_id, $
			n_mag_u,n_mag_g,n_mag_r,n_mag_i,n_mag_z,$
			n_mag_err_u,n_mag_err_g,n_mag_err_r,n_mag_err_i,n_mag_err_z,$
			sdss_ra, sdss_dec, $
			n_type,$
			sdss_redshift,sdss_redshift_err,$
			n_absmag_U,n_absmag_G,n_absmag_R,n_absmag_I,n_absmag_Z,$
			n_kcorrU,n_kcorrG,n_kcorrR,n_kcorrI,n_kcorrZ,$
			sdss_search_id,sdss_objid

	endif else restore, savdir+'spiral_neighbors_willettk.sav'

	r = {objid:'', $
		ra:0., $
		dec:0., $
		z:0., $
		zerr:0., $
		search_id:0L, $
		p_el_debiased:0., $
		p_cs_debiased:0., $
		n_p_el_debiased:ptr_new(), $
		n_p_cs_debiased:ptr_new(), $
		n_matched_id:ptr_new(), $
		n_mag_u:ptr_new(), $
		n_mag_g:ptr_new(), $
		n_mag_r:ptr_new(), $
		n_mag_i:ptr_new(), $
		n_mag_z:ptr_new(), $
		n_mag_err_u:ptr_new(), $
		n_mag_err_g:ptr_new(), $
		n_mag_err_r:ptr_new(), $
		n_mag_err_i:ptr_new(), $
		n_mag_err_z:ptr_new(), $
		n_absmag_u:ptr_new(), $
		n_absmag_g:ptr_new(), $
		n_absmag_r:ptr_new(), $
		n_absmag_i:ptr_new(), $
		n_absmag_z:ptr_new(), $
		n_redshift:ptr_new(), $
		n_redshift_err:ptr_new(), $
		n_ra:ptr_new(), $
		n_dec:ptr_new(),$
		bg:0,$
		area:0.,$
		bg_mag:ptr_new(), $
		cmag:0.,$
		n500:0, $
		n500_cmag:0}

	; Match the SDSS data to its galaxy parent ID

	unique_search_id1 = unique(search_id)

	; Remove galaxies that had no neighbors within 3 arcmin (kluge for now; should retain info)

	match,unique_search_id1,n_search_id,a,b, count=matchcount
	if matchcount ne n_elements(n_search_id) then begin
		nogals = setdifference(unique_search_id1, n_search_id)
		nogalsind = [0]
		for i=0, n_elements(nogals)-1 do nogalsind = [nogalsind,where(search_id eq nogals[i])]
		nogalsind = nogalsind[1:n_elements(nogalsind)-1]
		galsind = setdifference(lindgen(n_elements(ra)),nogalsind)
	
		ra = ra[galsind]
		dec = dec[galsind]
		search_id = search_id[galsind]
		redshift = redshift[galsind]
		redshift_err = redshift_err[galsind]
		p_el_debiased = p_el_debiased[galsind]
		p_cs_debiased = p_cs_debiased[galsind]

	endif

	; Isolate unique parent galaxies

	uobjid = unique(objid)
	nb = n_elements(uobjid)
	unique_search_id2 = unique(search_id)

	ptrarr_mag_u = ptrarr(nb,/allocate_heap)
	ptrarr_mag_g = ptrarr(nb,/allocate_heap)
	ptrarr_mag_r = ptrarr(nb,/allocate_heap)
	ptrarr_mag_i = ptrarr(nb,/allocate_heap)
	ptrarr_mag_z = ptrarr(nb,/allocate_heap)
	ptrarr_mag_err_u = ptrarr(nb,/allocate_heap)
	ptrarr_mag_err_g = ptrarr(nb,/allocate_heap)
	ptrarr_mag_err_r = ptrarr(nb,/allocate_heap)
	ptrarr_mag_err_i = ptrarr(nb,/allocate_heap)
	ptrarr_mag_err_z = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_u = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_g = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_r = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_i = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_z = ptrarr(nb,/allocate_heap)
	ptrarr_redshift = ptrarr(nb,/allocate_heap)
	ptrarr_redshift_err = ptrarr(nb,/allocate_heap)
	ptrarr_ra = ptrarr(nb,/allocate_heap)
	ptrarr_dec = ptrarr(nb,/allocate_heap)
	ptrarr_p_el_debiased = ptrarr(nb,/allocate_heap)
	ptrarr_p_cs_debiased = ptrarr(nb,/allocate_heap)
	ptrarr_matched_id = ptrarr(nb,/allocate_heap)

	for i=0, nb-1 do begin
		temp_sind = where(n_search_id eq unique_search_id2[i], ts_count)
		if ts_count gt 0 then begin
			*(ptrarr_mag_u[i]) = n_mag_u[temp_sind]
			*(ptrarr_mag_g[i]) = n_mag_g[temp_sind]
			*(ptrarr_mag_r[i]) = n_mag_r[temp_sind]
			*(ptrarr_mag_i[i]) = n_mag_i[temp_sind]
			*(ptrarr_mag_z[i]) = n_mag_z[temp_sind]
			*(ptrarr_mag_err_u[i]) = n_mag_err_u[temp_sind]
			*(ptrarr_mag_err_g[i]) = n_mag_err_g[temp_sind]
			*(ptrarr_mag_err_r[i]) = n_mag_err_r[temp_sind]
			*(ptrarr_mag_err_i[i]) = n_mag_err_i[temp_sind]
			*(ptrarr_mag_err_z[i]) = n_mag_err_z[temp_sind]
			*(ptrarr_absmag_u[i]) = n_absmag_u[temp_sind]
			*(ptrarr_absmag_g[i]) = n_absmag_g[temp_sind]
			*(ptrarr_absmag_r[i]) = n_absmag_r[temp_sind]
			*(ptrarr_absmag_i[i]) = n_absmag_i[temp_sind]
			*(ptrarr_absmag_z[i]) = n_absmag_z[temp_sind]
			*(ptrarr_redshift[i]) = n_redshift[temp_sind]
			*(ptrarr_redshift_err[i]) = n_redshift_err[temp_sind]
			*(ptrarr_ra[i]) = n_ra[temp_sind]
			*(ptrarr_dec[i]) = n_dec[temp_sind]
			*(ptrarr_p_el_debiased[i]) = n_p_el_debiased[temp_sind]
			*(ptrarr_p_cs_debiased[i]) = n_p_cs_debiased[temp_sind]
			*(ptrarr_matched_id[i]) = n_matched_id[temp_sind]
		endif else print,'No match for ID ',unique_search_id2[i], '#',i
	endfor 

	spiral = replicate(r,nb)

	spiral.n_mag_u = ptrarr_mag_u
	spiral.n_mag_g = ptrarr_mag_g
	spiral.n_mag_r = ptrarr_mag_r
	spiral.n_mag_i = ptrarr_mag_i
	spiral.n_mag_z = ptrarr_mag_z
	spiral.n_mag_err_u = ptrarr_mag_err_u
	spiral.n_mag_err_g = ptrarr_mag_err_g
	spiral.n_mag_err_r = ptrarr_mag_err_r
	spiral.n_mag_err_i = ptrarr_mag_err_i
	spiral.n_mag_err_z = ptrarr_mag_err_z
	spiral.n_absmag_u = ptrarr_absmag_u
	spiral.n_absmag_g = ptrarr_absmag_g
	spiral.n_absmag_r = ptrarr_absmag_r
	spiral.n_absmag_i = ptrarr_absmag_i
	spiral.n_absmag_z = ptrarr_absmag_z
	spiral.n_redshift = ptrarr_redshift
	spiral.n_redshift_err = ptrarr_redshift_err
	spiral.n_ra = ptrarr_ra
	spiral.n_dec = ptrarr_dec
	spiral.n_p_el_debiased = ptrarr_p_el_debiased
	spiral.n_p_cs_debiased = ptrarr_p_cs_debiased
	spiral.n_matched_id=ptrarr_matched_id

	eind = uniq(objid)
	bindarr = [uniq(objid),n_elements(objid)]

	spiral.objid=objid[eind]
	spiral.ra=ra[eind]
	spiral.dec=dec[eind]
	spiral.search_id=search_id[eind]
	spiral.z=redshift[eind]
	spiral.zerr=redshift_err[eind]
	spiral.p_el_debiased=p_el_debiased[eind]
	spiral.p_cs_debiased=p_cs_debiased[eind]

	;ptrarr_matched_id = ptrarr(nb,/allocate_heap)
	;for i=0,nb-1 do begin
	;	*(ptrarr_matched_id[i]) = n_matched_id[bindarr[i]:bindarr[i+1]-1]
	;endfor
	;spiral.matched_id=ptrarr_matched_id

	; For each galaxy, find the 500 kpc angular distance [arcsec]

	ang500 = zang(500., spiral.z,/silent,/wmap7)
	ang1000 = zang(1000., spiral.z,/silent,/wmap7)
	Mstar = -21.22
	sdss_lim = 22.2

	ptrarr_bg_mag = ptrarr(nb,/allocate_heap)

	for i=0,nb-1 do begin

		; Find number of neighbors within projected 500 Mpc

		gcirc, 2, spiral[i].ra, spiral[i].dec, *(spiral[i].n_ra), *(spiral[i].n_dec), angdist
		ind3 = where(angdist lt ang500[i], ind3count)
		spiral[i].n500 = ind3count

		; Find the number of neighbors within 500 Mpc and APPARENT magnitude brighter than counting mag

		distance = lumdist(spiral[i].z,/silent,/wmap7)
		app_mstar = (Mstar+2) + 5*alog10(distance*1e6) - 5
		countingmag_app = sdss_lim < app_mstar
		ind3_cmagapp = where((angdist lt ang500[i]) and (*(spiral[i].n_mag_r) lt countingmag_app), cmagcount_app)

			; Set the counting magnitude to ABSOLUTE, rather than apparent magnitude

			sdss_lim_abs = sdss_lim - 5*alog10(distance*1e6) + 5
			abs_mstar = Mstar + 2
			countingmag_abs = sdss_lim_abs < abs_mstar
			abs_rmag = *(spiral[i].n_mag_r) - 5*alog10(distance*1e6) + 5
			abs_rmag_kcorr = *(spiral[i].n_absmag_r) - 5*alog10(distance*1e6) + 5
			ind3_cmagabs = where((angdist lt ang500[i]) and $
				(abs_rmag_kcorr lt countingmag_abs), cmagcount_abs)

		n500_cmag = cmagcount_abs
		cmag = countingmag_abs

		spiral[i].n500_cmag = n500_cmag
		spiral[i].cmag = cmag

		; Determine the background galaxy density (Smith, O'Dea, & Baum 1995)

		bg_inds = where((angdist gt ang1000[i]) and $
			(*(spiral[i].n_mag_r) lt cmag), bgcount)
		bg_mag = (*(spiral[i].n_mag_r))[bg_inds]
		area = !dpi * (3.^2 - (ang1000[i]/60.)^2)	; arcmin
		spiral[i].bg = bgcount
		spiral[i].area = area
		*(ptrarr_bg_mag[i]) = bg_mag
	endfor

	spiral.bg_mag = ptrarr_bg_mag
		
	; Radio galaxy morphological classes

	save, filename=savdir+'spiral_structure.sav', spiral

	if ~keyword_set(quiet) then begin
		print,''
		print,'spiral galaxy data structure written to '+savdir+'spiral_structure.sav'
		print,''
	endif

	if keyword_set(stop) then stop

end

