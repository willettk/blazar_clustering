;+
; NAME:
;       
;	TCRR_STRMAKE
;
; PURPOSE:
;
;	Create structure with SDSS, 3CRR information on morphologically-classified radio galaxies
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
;	IDL> tcrr_strmake
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Dec 11
;-

pro tcrr_strmake, stop=stop, readbfiles=readbfiles

	heap_gc		; Clean up heap variables

	if n_elements(phi_star) eq 0 then phi_star = 0.5

	csvdir = '~/Astronomy/Research/blazars/csv/'
	savdir = '~/Astronomy/Research/blazars/sav/'

	tcrrfile = csvdir+'tcrr_neighbors3_willettk.csv'
	magfile = csvdir+'mags_3crr_willettk.csv'
	kcorrfile = csvdir+'kcorr_3crr_willettk.csv'

	if keyword_set(readbfiles) then begin

		readcol, tcrrfile, $
			name_tcrr,name_iau,redshift,frclass,ra,dec,search_id,matched_id, $
			format='a,a,f,i,f,f,i,a', $
			skipline=1

		readcol, magfile, $
			mag_u,mag_g,mag_r,mag_i,mag_z,$
			mag_err_u,mag_err_g,mag_err_r,mag_err_i,mag_err_z,$
			mag_ra,mag_dec,$
			mag_type,$
			mag_redshift,mag_redshift_err,$
			mag_search_id, $
			mag_objid, $
			format='f,f,f,f,f,f,f,f,f,f,f,f,i,f,f,i,a', $
			skipline=1

		readcol, kcorrfile, $
			absmag_kcorr_u, $ 
			absmag_kcorr_g, $ 
			absmag_kcorr_r, $ 
			absmag_kcorr_i, $ 
			absmag_kcorr_z, $ 
			kcorr_u, $ 
			kcorr_g, $ 
			kcorr_r, $ 
			kcorr_i, $ 
			kcorr_z, $ 
			searchid_kcorr, $
			objid_kcorr, $
			format='f,f,f,f,f,f,f,f,f,f,i,a', $
			skipline=1

		save, filename=savdir+'tcrr_neighbors3_willettk.sav', $
			name_tcrr,name_iau,redshift,frclass,ra,dec,search_id,matched_id, $
			mag_u,mag_g,mag_r,mag_i,mag_z,$
			mag_err_u,mag_err_g,mag_err_r,mag_err_i,mag_err_z,$
			mag_ra,mag_dec,$
			mag_type,$
			mag_redshift,mag_redshift_err,$
			mag_search_id, $
			mag_objid, $
			absmag_kcorr_u, $ 
			absmag_kcorr_g, $ 
			absmag_kcorr_r, $ 
			absmag_kcorr_i, $ 
			absmag_kcorr_z, $ 
			kcorr_u, $ 
			kcorr_g, $ 
			kcorr_r, $ 
			kcorr_i, $ 
			kcorr_z, $ 
			searchid_kcorr, $
			objid_kcorr

	endif else restore, savdir+'tcrr_neighbors3_willettk.sav'

	b = {name_tcrr:'', name_iau:'', $
		z:0., $
		frclass:0, $
		ra:0., dec:0., $
		search_id:0L, $
		matched_id:ptr_new(), $
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
		n_absmag_kcorr_u:ptr_new(), $
		n_absmag_kcorr_g:ptr_new(), $
		n_absmag_kcorr_r:ptr_new(), $
		n_absmag_kcorr_i:ptr_new(), $
		n_absmag_kcorr_z:ptr_new(), $
		n_redshift:ptr_new(), $
		n_redshift_err:ptr_new(), $
		n_ra:ptr_new(), $
		n_dec:ptr_new(),$
		bg:0,$
		area:0.,$
		bg_mag:ptr_new(), $
		cmag:0.,$
		cmag_app:0.,$
		n500:0, $
		n500_cmag:0}

	; Match the SDSS data to its blazar parent ID

	unique_searchid1 = unique(search_id,/sort)

	; Remove blazars that had no galaxies within 3 arcmin (kluge for now; should retain info)

	nogals = setdifference(unique_searchid1, mag_search_id)
	nogalsind = [0]
	for i=0, n_elements(nogals)-1 do nogalsind = [nogalsind,where(search_id eq nogals[i])]
	nogalsind = nogalsind[1:n_elements(nogalsind)-1]
	galsind = setdifference(lindgen(n_elements(name_tcrr)),nogalsind)

		name_tcrr = name_tcrr[galsind]
		name_iau = name_iau[galsind]
		ra = ra[galsind]
		dec = dec[galsind]
		z = redshift[galsind]
		frclass = frclass[galsind]
		matched_id = matched_id[galsind]
		search_id = search_id[galsind]

	uname = unique(name_iau,/sort)
	nb = n_elements(uname)
	unique_searchid2 = unique(search_id,/sort)

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
	ptrarr_absmag_kcorr_u = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_kcorr_g = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_kcorr_r = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_kcorr_i = ptrarr(nb,/allocate_heap)
	ptrarr_absmag_kcorr_z = ptrarr(nb,/allocate_heap)
	ptrarr_redshift = ptrarr(nb,/allocate_heap)
	ptrarr_redshift_err = ptrarr(nb,/allocate_heap)
	ptrarr_ra = ptrarr(nb,/allocate_heap)
	ptrarr_dec = ptrarr(nb,/allocate_heap)

	for i=0, nb-1 do begin
		temp_sind = where(mag_search_id eq unique_searchid2[i], ts_count)
		if ts_count gt 0 then begin
			*(ptrarr_mag_u[i]) = mag_u[temp_sind]
			*(ptrarr_mag_g[i]) = mag_g[temp_sind]
			*(ptrarr_mag_r[i]) = mag_r[temp_sind]
			*(ptrarr_mag_i[i]) = mag_i[temp_sind]
			*(ptrarr_mag_z[i]) = mag_z[temp_sind]
			*(ptrarr_mag_err_u[i]) = mag_err_u[temp_sind]
			*(ptrarr_mag_err_g[i]) = mag_err_g[temp_sind]
			*(ptrarr_mag_err_r[i]) = mag_err_r[temp_sind]
			*(ptrarr_mag_err_i[i]) = mag_err_i[temp_sind]
			*(ptrarr_mag_err_z[i]) = mag_err_z[temp_sind]
			*(ptrarr_redshift[i]) = mag_redshift[temp_sind]
			*(ptrarr_redshift_err[i]) = mag_redshift_err[temp_sind]
			*(ptrarr_ra[i]) = mag_ra[temp_sind]
			*(ptrarr_dec[i]) = mag_dec[temp_sind]
		endif else print,'No match for ID ',unique_searchid2[i], '#',i

		temp_kind = where(searchid_kcorr eq unique_searchid2[i], tk_count)
		if tk_count gt 0 then begin
			*(ptrarr_absmag_kcorr_u[i]) = absmag_kcorr_u[temp_kind]
			*(ptrarr_absmag_kcorr_g[i]) = absmag_kcorr_g[temp_kind]
			*(ptrarr_absmag_kcorr_r[i]) = absmag_kcorr_r[temp_kind]
			*(ptrarr_absmag_kcorr_i[i]) = absmag_kcorr_i[temp_kind]
			*(ptrarr_absmag_kcorr_z[i]) = absmag_kcorr_z[temp_kind]
		endif else print,'No match for K-corrected ID ',unique_searchid2[i], '#',i
	endfor 

	tcrr = replicate(b,nb)

	tcrr.n_mag_u = ptrarr_mag_u
	tcrr.n_mag_g = ptrarr_mag_g
	tcrr.n_mag_r = ptrarr_mag_r
	tcrr.n_mag_i = ptrarr_mag_i
	tcrr.n_mag_z = ptrarr_mag_z
	tcrr.n_mag_err_u = ptrarr_mag_err_u
	tcrr.n_mag_err_g = ptrarr_mag_err_g
	tcrr.n_mag_err_r = ptrarr_mag_err_r
	tcrr.n_mag_err_i = ptrarr_mag_err_i
	tcrr.n_mag_err_z = ptrarr_mag_err_z
	tcrr.n_absmag_kcorr_u = ptrarr_absmag_kcorr_u
	tcrr.n_absmag_kcorr_g = ptrarr_absmag_kcorr_g
	tcrr.n_absmag_kcorr_r = ptrarr_absmag_kcorr_r
	tcrr.n_absmag_kcorr_i = ptrarr_absmag_kcorr_i
	tcrr.n_absmag_kcorr_z = ptrarr_absmag_kcorr_z
	tcrr.n_redshift = ptrarr_redshift
	tcrr.n_redshift_err = ptrarr_redshift_err
	tcrr.n_ra = ptrarr_ra
	tcrr.n_dec = ptrarr_dec

	tcrrind = uniq(name_iau,sort(search_id))

	tcrr.name_tcrr=name_tcrr[tcrrind]
	tcrr.name_iau=name_iau[tcrrind]
	tcrr.frclass=frclass[tcrrind]
	tcrr.ra=ra[tcrrind]
	tcrr.dec=dec[tcrrind]
	tcrr.z=z[tcrrind]
	tcrr.search_id=search_id[tcrrind]

	ptrarr_matchedid = ptrarr(nb,/allocate_heap)
	for i=0,nb-1 do begin
		*(ptrarr_matchedid[i]) = matched_id[where(search_id eq search_id[tcrrind[i]])]
	endfor
	tcrr.matched_id=ptrarr_matchedid

	; For each galaxy, find the 500 kpc angular distance [arcsec]

	ang500 = zang(500., tcrr.z,/silent,/wmap7)
	ang1000 = zang(1000., tcrr.z,/silent,/wmap7)
	Mstar = -21.1837		; Blanton et al. (2003)
	sdss_lim = 22.2

	ptrarr_bg_mag = ptrarr(nb,/allocate_heap)

	for i=0,nb-1 do begin

		; Find number of neighbors within projected 500 Mpc

		gcirc, 2, tcrr[i].ra, tcrr[i].dec, *(tcrr[i].n_ra), *(tcrr[i].n_dec), angdist
		ind3 = where(angdist lt ang500[i], ind3count)
		tcrr[i].n500 = ind3count

		; Find the number of neighbors within 500 Mpc and APPARENT magnitude brighter than counting mag

		distance = lumdist(tcrr[i].z,/silent,/wmap7)
		app_mstar = (Mstar+2) + 5*alog10(distance*1e6) - 5
		countingmag_app = sdss_lim < app_mstar
		ind3_cmagapp = where((angdist lt ang500[i]) and (*(tcrr[i].n_mag_r) lt countingmag_app), cmagcount_app)
		tcrr[i].cmag_app = countingmag_app

			; Set the counting magnitude to ABSOLUTE, rather than apparent magnitude

			sdss_lim_abs = sdss_lim - 5*alog10(distance*1e6) + 5
			abs_mstar = Mstar + 2
			countingmag_abs = sdss_lim_abs < abs_mstar
			abs_rmag = *(tcrr[i].n_mag_r) - 5*alog10(distance*1e6) + 5
			abs_rmag_kcorr = *(tcrr[i].n_absmag_kcorr_r)
			ind3_cmagabs = where((angdist lt ang500[i]) and $
				(abs_rmag_kcorr lt countingmag_abs), cmagcount_abs)

		n500_cmag = cmagcount_abs
		cmag = countingmag_abs

		tcrr[i].n500_cmag = n500_cmag
		tcrr[i].cmag = cmag

		; Determine the background galaxy density (Smith, O'Dea, & Baum 1995)

		bg_inds = where((angdist gt ang1000[i]) and $
			(*(tcrr[i].n_mag_r) lt cmag), bgcount)
		bg_mag = (*(tcrr[i].n_mag_r))[bg_inds]
		area = !dpi * (3.^2 - (ang1000[i]/60.)^2)	; arcmin
		tcrr[i].bg = bgcount
		tcrr[i].area = area
		*(ptrarr_bg_mag[i]) = bg_mag
	endfor

	tcrr.bg_mag = ptrarr_bg_mag
		
	; Radio galaxy morphological classes

	fr1 = where(tcrr.frclass eq 1)
	fr2 = where(tcrr.frclass eq 2)

	save, filename=savdir+'tcrr_structure.sav', tcrr, $
		fr1, fr2

	print,''
	print,'Radio galaxy data structure written to '+savdir+'tcrr_structure.sav'
	print,''

	if keyword_set(stop) then stop

end
