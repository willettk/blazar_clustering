;+
; NAME:
;       
;	KIMBALL_STRMAKE
;
; PURPOSE:
;
;	Create structure with SDSS, radio morphology information on galaxies from the Kimball & Ivezic radio/SDSS overlap survey
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
;	IDL> bz_strmake
;
; NOTES:
;
;	Reference: http://adsabs.harvard.edu/abs/2008AJ....136..684K
;	Source: http://www.astro.washington.edu/users/akimball/radiocat/
;
;	This data is from catalog 5, consisting of the galaxies detected in all four radio surveys (FIRST, NVSS, WENSS, GB6) and 
;		have spectroscopic detections (with redshifts) from the SDSS.
;
; REVISION HISTORY
;       Written by K. Willett                Nov 11
;-

; Begin program

pro kimball_strmake, stop=stop, readbfiles=readbfiles

	heap_gc		; Clean up heap variables

	if n_elements(phi_star) eq 0 then phi_star = 0.5

	csvdir = '~/Astronomy/Research/blazars/csv/'
	savdir = '~/Astronomy/Research/blazars/sav/'

	neighbors_file = csvdir+'kimball_3arcmin_willettk.csv'
	sdss_file = csvdir+'kimball_sdss_willettk.csv'

	if keyword_set(readbfiles) then begin

		readcol, neighbors_file, $
			kimball_id,$
			ra,dec,$
			first_peak_flux,first_flux,nvss_flux,$
			kimball_rmag,$
			spec_type,spec_redshift,spec_redshifterr,spec_redshift_warning,$
			searchid,match_objid, $
			format='a,f,f,f,f,f,f,i,f,f,l,i,a', $
			skipline=1

		readcol, sdss_file, $
			n_mag_u,n_mag_g,n_mag_r,n_mag_i,n_mag_z,$
			n_mag_err_u,n_mag_err_g,n_mag_err_r,n_mag_err_i,n_mag_err_z,$
			n_ra, n_dec, $
			n_type,$
			n_redshift,n_redshift_err,$
			n_absmag_U,n_absmag_G,n_absmag_R,n_absmag_I,n_absmag_Z,$
			n_kcorr_U,n_kcorr_G,n_kcorr_R,n_kcorr_I,n_kcorr_Z,$
			n_searchid,n_objid, $
			format='f,f,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f,f,f,f,f,f,f,f,f,i,a', $
			skipline=1

		; Kluge to ensure unique search IDs

		nlines_neighbors = file_lines('~/Astronomy/Research/blazars/kimball_3arcmin_part1_willettk.csv')
		nlines_sdss = file_lines('~/Astronomy/Research/blazars/kimball_part1_sdss_willettk.csv')

		searchid[0:nlines_neighbors-2] += 5000
		n_searchid[0:nlines_sdss-2] += 5000

		save, filename=savdir+'kimball_neighbors_all.sav', $
			kimball_id, ra, dec, $
			searchid, match_objid, $
			first_peak_flux,first_flux,nvss_flux,$
			kimball_rmag,$
			spec_type,spec_redshift,spec_redshifterr,spec_redshift_warning,$
			n_mag_u,n_mag_g,n_mag_r,n_mag_i,n_mag_z,$
			n_mag_err_u,n_mag_err_g,n_mag_err_r,n_mag_err_i,n_mag_err_z,$
			n_ra,n_dec,$
			n_type,$
			n_redshift,n_redshift_err,$
			n_searchid, $
			n_objid, $
			n_absmag_u, $ 
			n_absmag_g, $ 
			n_absmag_r, $ 
			n_absmag_i, $ 
			n_absmag_z, $ 
			n_kcorr_u, $ 
			n_kcorr_g, $ 
			n_kcorr_r, $ 
			n_kcorr_i, $ 
			n_kcorr_z, $ 
			n_searchid, $
			n_objid

	endif else restore, savdir+'kimball_neighbors_all.sav'

	print,'BP 1'

	b = {kimball_id:'', $
		ra:0., dec:0., $
		z:0., zerr:0., $
		first_peak_flux: 0., $
		first_flux: 0., $
		nvss_flux: 0., $
		kimball_rmag: 0., $
		radio_morph:'', $
		searchid:0L, $
		match_objid:ptr_new(), $
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
		cmag_app:0.,$
		n500:0, $
		n500_cmag:0}

	print,'BP 2'

	; Match the SDSS data to its blazar parent ID

	unique_searchid1 = unique(searchid,/sort)
	unique_nsearchid1 = unique(n_searchid,/sort)

	print,'BP 3'

	; Remove blazars that had no galaxies within the search radii (kluge for now; should retain info)

	if n_elements(unique_searchid1) ne n_elements(unique_nsearchid1) then begin

		nogals = setdifference(unique_searchid1, unique_nsearchid1)
		nogalsind = [0]
		for i=0, n_elements(nogals)-1 do nogalsind = [nogalsind,where(searchid eq nogals[i])]
		nogalsind = nogalsind[1:n_elements(nogalsind)-1]
		galsind = setdifference(lindgen(n_elements(kimball_id)),nogalsind)
	
		print,'BP 4'
	
		kimball_id = kimball_id[galsind]
		ra = ra[galsind]
		dec = dec[galsind]
		z = spec_redshift[galsind]
		zerr = spec_redshifterr[galsind]
		first_peak_flux = first_peak_flux[galsind]
		first_flux = first_flux[galsind]
		nvss_flux = nvss_flux[galsind]
		match_objid = match_objid[galsind]
		searchid = searchid[galsind]
	
	endif else begin
		z = spec_redshift
		zerr = spec_redshifterr
	endelse

	print,'BP 5'

	; Find unique parent blazars

	name_uniq = unique(kimball_id,bcount,/sort)
	unique_searchid2 = unique(searchid,ucount,/sort)		

	nb = n_elements(name_uniq)
	nu = n_elements(unique_searchid2)

	if nb ne nu then message,'Number of target names is not the same as number of search IDs',/info

	print,'BP 6'

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

	print,'BP 7'


	print,'BP 8'

	for i=0, nb-1 do begin
		temp_sind = where(n_searchid eq unique_searchid2[i], ts_count)
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
			*(ptrarr_redshift[i]) = n_redshift[temp_sind]
			*(ptrarr_redshift_err[i]) = n_redshift_err[temp_sind]
			*(ptrarr_ra[i]) = n_ra[temp_sind]
			*(ptrarr_dec[i]) = n_dec[temp_sind]
			*(ptrarr_absmag_u[i]) = n_absmag_u[temp_sind]
			*(ptrarr_absmag_g[i]) = n_absmag_g[temp_sind]
			*(ptrarr_absmag_r[i]) = n_absmag_r[temp_sind]
			*(ptrarr_absmag_i[i]) = n_absmag_i[temp_sind]
			*(ptrarr_absmag_z[i]) = n_absmag_z[temp_sind]
		endif else print,'No match for search ID ',unique_searchid2[i], '#',i
	endfor 

	kimball = replicate(b,nb)

	kimball.n_mag_u = ptrarr_mag_u
	kimball.n_mag_g = ptrarr_mag_g
	kimball.n_mag_r = ptrarr_mag_r
	kimball.n_mag_i = ptrarr_mag_i
	kimball.n_mag_z = ptrarr_mag_z
	kimball.n_mag_err_u = ptrarr_mag_err_u
	kimball.n_mag_err_g = ptrarr_mag_err_g
	kimball.n_mag_err_r = ptrarr_mag_err_r
	kimball.n_mag_err_i = ptrarr_mag_err_i
	kimball.n_mag_err_z = ptrarr_mag_err_z
	kimball.n_absmag_u = ptrarr_absmag_u
	kimball.n_absmag_g = ptrarr_absmag_g
	kimball.n_absmag_r = ptrarr_absmag_r
	kimball.n_absmag_i = ptrarr_absmag_i
	kimball.n_absmag_z = ptrarr_absmag_z
	kimball.n_redshift = ptrarr_redshift
	kimball.n_redshift_err = ptrarr_redshift_err
	kimball.n_ra = ptrarr_ra
	kimball.n_dec = ptrarr_dec

	rgind = uniq(searchid,sort(searchid))

	kimball.kimball_id=kimball_id[rgind]
	kimball.ra=ra[rgind]
	kimball.dec=dec[rgind]
	kimball.z=z[rgind]
	kimball.zerr=zerr[rgind]
	kimball.first_peak_flux=first_peak_flux[rgind]
	kimball.first_flux=first_flux[rgind]
	kimball.nvss_flux=nvss_flux[rgind]
	kimball.searchid=searchid[rgind]

	ptrarr_matchedid = ptrarr(nb,/allocate_heap)
	for i=0,nb-1 do begin
		*(ptrarr_matchedid[i]) = match_objid[where(searchid eq searchid[rgind[i]])]
	endfor
	kimball.match_objid=ptrarr_matchedid

	; Radio morphology

	tfirst = -2.5 * alog10(kimball.first_flux * 1d-3 /3631.)
	tnvss = -2.5 * alog10(kimball.nvss_flux * 1d-3 /3631.)
	
	deltat = tfirst - tnvss
	theta = sqrt(kimball.first_flux / kimball.first_peak_flux)
	
	ind_complex = where(deltat gt 0.35, ncomplex)
	ind_compact = where(deltat lt 0.35 and alog10(theta^2) lt 0.05, ncompact)
	ind_resolved = where(deltat lt 0.35 and alog10(theta^2) gt 0.05, nresolved)

	kimball[ind_complex].radio_morph  = 'complex'
	kimball[ind_compact].radio_morph  = 'compact'
	kimball[ind_resolved].radio_morph = 'resolved'

	; For each galaxy, find the 500 kpc angular distance [arcsec]

	ang500 = zang(500., kimball.z,/silent,/wmap7)
	ang1000 = zang(1000., kimball.z,/silent,/wmap7)
	Mr_star = -21.1837		; Blanton et al. (2003)
	sdss_lim = 22.2			; SDSS apparent limiting magnitude

	ptrarr_bg_mag = ptrarr(nb,/allocate_heap)

	for i=0,nb-1 do begin

		; Find number of neighbors within projected 500 Mpc

		gcirc, 2, kimball[i].ra, kimball[i].dec, *(kimball[i].n_ra), *(kimball[i].n_dec), angdist
		ind3 = where(angdist lt ang500[i], ind3count)
		kimball[i].n500 = ind3count

		; Find the number of neighbors within 500 Mpc and APPARENT magnitude brighter than counting mag

		distance = lumdist(kimball[i].z,/silent,/wmap7)
		app_mstar = (Mr_star+2) + 5*alog10(distance*1e6) - 5
		countingmag_app = sdss_lim < app_mstar
		ind3_cmagapp = where((angdist lt ang500[i]) and (*(kimball[i].n_mag_r) lt countingmag_app), cmagcount_app)
		kimball[i].cmag_app = countingmag_app

			; Set the counting magnitude to ABSOLUTE, rather than apparent magnitude

			sdss_lim_abs = sdss_lim - 5*alog10(distance*1e6) + 5
			abs_mstar = Mr_star + 2
			countingmag_abs = sdss_lim_abs < abs_mstar
			abs_rmag = *(kimball[i].n_absmag_r)
			ind3_cmagabs = where((angdist lt ang500[i]) and $
				(abs_rmag lt countingmag_abs), cmagcount_abs)

		n500_cmag = cmagcount_abs
		cmag = countingmag_abs

		kimball[i].n500_cmag = n500_cmag
		kimball[i].cmag = cmag

		; Determine the background galaxy density in the vicinity (Smith, O'Dea, & Baum 1995)

		bg_inds = where((angdist gt ang1000[i]) and $
			(*(kimball[i].n_absmag_r) lt cmag), bgcount)
		bg_mag = (*(kimball[i].n_absmag_r))[bg_inds]
		if kimball[i].z gt 0.162 then outer_radius = 3. else outer_radius = 10.
		area = !dpi * (outer_radius^2 - (ang1000[i]/60.)^2)	; arcmin
		kimball[i].bg = bgcount
		kimball[i].area = area
		*(ptrarr_bg_mag[i]) = bg_mag
	endfor

	kimball.bg_mag = ptrarr_bg_mag
		
	save, filename=savdir+'kimball_structure.sav', kimball

	print,''
	print,'Structure written to '+savdir+'kimball_structure.sav'
	print,''

	if keyword_set(stop) then stop

end

