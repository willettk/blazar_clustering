
;+
; NAME:
;       
;	RANDOM_STRMAKE
;
; PURPOSE:
;
;	Create structure with SDSS information on neighbors for randomly-selected points
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
;	IDL> random_strmake
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Dec 11
;-

pro random_strmake, stop=stop, readbfiles=readbfiles, quiet=quiet

	heap_gc		; Clean up heap variables

	if n_elements(phi_star) eq 0 then phi_star = 0.5

	csvdir = '~/Astronomy/Research/blazars/csv/'
	savdir = '~/Astronomy/Research/blazars/sav/'

	random_file_3 = csvdir+'rand_3arcmin_n_willettk.csv'
	random_sdssfile_3 = csvdir+'rand_3arcmin_sdss_willettk.csv'
	random_file_10 = csvdir+'rand_10arcmin_n_willettk.csv'
	random_sdssfile_10 = csvdir+'rand_10arcmin_sdss_willettk.csv'

	random_neighbors_all = csvdir+'random_neighbors_all.csv'
	random_sdss_all = csvdir+'random_sdss_all.csv'

	spawn,'cat '+random_file_3+' '+random_file_10+' > '+random_neighbors_all
	spawn,'cat '+random_sdssfile_3+' '+random_sdssfile_10+' > '+random_sdss_all

	if keyword_set(readbfiles) then begin

		readcol, random_neighbors_all, $
			rname,ra,dec,z,searchid,match_objid, $
			format='a,f,f,f,i,a', $
			skipline=1
		bbadind=where(strmid(match_objid,19,4) eq 'name')
		match_objid[bbadind] = strmid(match_objid[bbadind],0,19)

		readcol, random_sdss_all, $
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

		sbadind=where(strmid(n_objid,19,1) eq 'u')
		n_objid[sbadind] = strmid(n_objid[sbadind],0,19)

		; Kluge to ensure unique search IDs

		nlines_bgb3 = file_lines(random_file_3)
		nlines_sdss3 = file_lines(random_sdssfile_3)

		searchid[0:nlines_bgb3-2] += 5000
		n_searchid[0:nlines_sdss3-2] += 5000

		; Load information for the control fields

		; North one degree from blazar
		file_connorth_3arcmin_neighbors = csvdir+'rand_3arcmin_connorth_n_willettk.csv'
		file_connorth_3arcmin_sdss = csvdir+'rand_3arcmin_connorth_sdss_willettk.csv'
		file_connorth_10arcmin_neighbors = csvdir+'rand_10arcmin_connorth_n_willettk.csv'
		file_connorth_10arcmin_sdss = csvdir+'rand_10arcmin_connorth_sdss_willettk.csv'

		; South one degree from blazar
		file_consouth_3arcmin_neighbors = csvdir+'rand_3arcmin_consouth_n_willettk.csv'
		file_consouth_3arcmin_sdss = csvdir+'rand_3arcmin_consouth_sdss_willettk.csv'
		file_consouth_10arcmin_neighbors = csvdir+'rand_10arcmin_consouth_n_willettk.csv'
		file_consouth_10arcmin_sdss = csvdir+'rand_10arcmin_consouth_sdss_willettk.csv'

		readcol, file_connorth_3arcmin_neighbors, $
			c3_rname,c3_ra,c3_dec,c3_z,c3_searchid,c3_matched_id, $
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_connorth_3arcmin_sdss, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_searchid,cn3_objid, $
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_connorth_10arcmin_neighbors, $
			c10_rname,c10_ra,c10_dec,c10_z,c10_searchid,c10_matched_id,$
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_connorth_10arcmin_sdss, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_searchid,cn10_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_consouth_3arcmin_neighbors, $
			c3d_rname,c3d_ra,c3d_dec,c3d_z,c3d_searchid,c3d_matched_id, $
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_consouth_3arcmin_sdss, $
			cn3d_u,cn3d_g,cn3d_r,cn3d_i,cn3d_z,cn3d_err_u,cn3d_err_g,cn3d_err_r,cn3d_err_i,cn3d_err_z,$
			cn3d_ra,cn3d_dec,cn3d_type,cn3d_redshift,cn3d_redshift_err,$
			cn3d_absmagU,cn3d_absmagG,cn3d_absmagR,cn3d_absmagI,cn3d_absmagZ,$
			cn3d_kcorrU,cn3d_kcorrG,cn3d_kcorrR,cn3d_kcorrI,cn3d_kcorrZ,$
			cn3d_searchid,cn3d_objid, $
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_consouth_10arcmin_neighbors, $
			c10d_rname,c10d_ra,c10d_dec,c10d_z,c10d_searchid,c10d_matched_id,$
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_consouth_10arcmin_sdss, $
			cn10d_u,cn10d_g,cn10d_r,cn10d_i,cn10d_z,cn10d_err_u,cn10d_err_g,cn10d_err_r,cn10d_err_i,cn10d_err_z,$
			cn10d_ra,cn10d_dec,cn10d_type,cn10d_redshift,cn10d_redshift_err,$
			cn10d_absmagU,cn10d_absmagG,cn10d_absmagR,cn10d_absmagI,cn10d_absmagZ,$
			cn10d_kcorrU,cn10d_kcorrG,cn10d_kcorrR,cn10d_kcorrI,cn10d_kcorrZ,$
			cn10d_searchid,cn10d_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		save, filename=savdir+'random_neighbors_all.sav', $
			rname, ra, dec, z, searchid, match_objid, $
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

		save, filename=savdir+'random_controlfields.sav', $
			c3_rname,c3_ra,c3_dec,c3_z,c3_searchid,c3_matched_id, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_searchid,cn3_objid, $
			c10_rname,c10_ra,c10_dec,c10_z,c10_searchid,c10_matched_id, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_searchid,cn10_objid, $
			c3d_rname,c3d_ra,c3d_dec,c3d_z,c3d_searchid,c3d_matched_id, $
			cn3d_u,cn3d_g,cn3d_r,cn3d_i,cn3d_z,cn3d_err_u,cn3d_err_g,cn3d_err_r,cn3d_err_i,cn3d_err_z,$
			cn3d_ra,cn3d_dec,cn3d_type,cn3d_redshift,cn3d_redshift_err,$
			cn3d_absmagU,cn3d_absmagG,cn3d_absmagR,cn3d_absmagI,cn3d_absmagZ,$
			cn3d_kcorrU,cn3d_kcorrG,cn3d_kcorrR,cn3d_kcorrI,cn3d_kcorrZ,$
			cn3d_searchid,cn3d_objid, $
			c10d_rname,c10d_ra,c10d_dec,c10d_z,c10d_searchid,c10d_matched_id, $
			cn10d_u,cn10d_g,cn10d_r,cn10d_i,cn10d_z,cn10d_err_u,cn10d_err_g,cn10d_err_r,cn10d_err_i,cn10d_err_z,$
			cn10d_ra,cn10d_dec,cn10d_type,cn10d_redshift,cn10d_redshift_err,$
			cn10d_absmagU,cn10d_absmagG,cn10d_absmagR,cn10d_absmagI,cn10d_absmagZ,$
			cn10d_kcorrU,cn10d_kcorrG,cn10d_kcorrR,cn10d_kcorrI,cn10d_kcorrZ,$
			cn10d_searchid,cn10d_objid
	
	endif else begin
		restore, savdir+'random_controlfields.sav'
		restore, savdir+'random_neighbors_all.sav'
	endelse

	r = {rname:'', $
		ra:0., dec:0., $
		z:0., $
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
		bg_mag:ptr_new(), $
		fieldsize:0.,$
		cmag:0.,$
		cmag_app:0.,$
		n500:0, $
		n500_cmag:0}

	; Match the SDSS data to its random point parent ID

	unique_searchid1 = unique(searchid,/sort)

	; Remove the random points that had no galaxies within 3 arcmin (kluge for now; should retain info)

	nogals = setdifference(unique_searchid1, n_searchid)
	nogalsind = [0]
	for i=0, n_elements(nogals)-1 do nogalsind = [nogalsind,where(searchid eq nogals[i])]
	nogalsind = nogalsind[1:n_elements(nogalsind)-1]
	galsind = setdifference(lindgen(n_elements(rname)),nogalsind)

		rname = rname[galsind]
		ra = ra[galsind]
		dec = dec[galsind]
		match_objid = match_objid[galsind]
		searchid = searchid[galsind]

	uname = unique(rname,/sort)
	nb = n_elements(uname)
	unique_searchid2 = unique(searchid,/sort)

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
		endif else print,'No match for ID ',unique_searchid2[i], '#',i
	endfor 

	random = replicate(r,nb)

	random.n_mag_u = ptrarr_mag_u
	random.n_mag_g = ptrarr_mag_g
	random.n_mag_r = ptrarr_mag_r
	random.n_mag_i = ptrarr_mag_i
	random.n_mag_z = ptrarr_mag_z
	random.n_mag_err_u = ptrarr_mag_err_u
	random.n_mag_err_g = ptrarr_mag_err_g
	random.n_mag_err_r = ptrarr_mag_err_r
	random.n_mag_err_i = ptrarr_mag_err_i
	random.n_mag_err_z = ptrarr_mag_err_z
	random.n_absmag_u = ptrarr_absmag_u
	random.n_absmag_g = ptrarr_absmag_g
	random.n_absmag_r = ptrarr_absmag_r
	random.n_absmag_i = ptrarr_absmag_i
	random.n_absmag_z = ptrarr_absmag_z
	random.n_redshift = ptrarr_redshift
	random.n_redshift_err = ptrarr_redshift_err
	random.n_ra = ptrarr_ra
	random.n_dec = ptrarr_dec

	ind_rstructure = uniq(searchid,sort(searchid))
	;bindarr = [uniq(rname),n_elements(rname)]

	random.rname=rname[ind_rstructure]
	random.ra=ra[ind_rstructure]
	random.dec=dec[ind_rstructure]
	random.z=z[ind_rstructure]
	random.searchid=searchid[ind_rstructure]

	ptrarr_matchedid = ptrarr(nb,/allocate_heap)
	for i=0,nb-1 do begin
		*(ptrarr_matchedid[i]) = match_objid[where(searchid eq searchid[ind_rstructure[i]])]
	endfor
	random.match_objid=ptrarr_matchedid

		; For the "redshift" of the point, randomly select from the distribution of neighbor galaxy redshifts

		;nm = n_elements(mag_redshift)
		;randarr = long(randomu(randseed,nb)*nm)
		;random.z = abs(mag_redshift[randarr])

	; For each galaxy, find the 500 kpc angular distance [arcsec]

	ang500 = zang(500., random.z,/silent,/wmap7)
	random.fieldsize = ang500
	Mr_star = -21.18
	sdss_lim = 22.2

	ptrarr_bg_mag = ptrarr(nb,/allocate_heap)

	; Background galaxy array generation

	file_rname_north = [c3_rname,c10_rname]
	searchid_north = [c3_searchid,c10_searchid+5000]
	con_ra_north = [c3_ra, c10_ra]
	con_dec_north = [c3_dec, c10_dec]

	ncon_searchid_north = [cn3_searchid,cn10_searchid+5000]
	ncon_ra_north = [cn3_ra, cn10_ra]
	ncon_dec_north = [cn3_dec, cn10_dec]
	ncon_absmagR_north = [cn3_absmagR, cn10_absmagR]

	file_rname_south = [c3d_rname,c10d_rname]
	searchid_south = [c3d_searchid,c10d_searchid+5000]
	con_ra_south = [c3d_ra, c10d_ra]
	con_dec_south = [c3d_dec, c10d_dec]

	ncon_searchid_south = [cn3d_searchid,cn10d_searchid+5000]
	ncon_ra_south = [cn3d_ra, cn10d_ra]
	ncon_dec_south = [cn3d_dec, cn10d_dec]
	ncon_absmagR_south = [cn3d_absmagR, cn10d_absmagR]

	for i=0,nb-1 do begin

		; Find number of neighbors within projected 500 Mpc

		gcirc, 2, random[i].ra, random[i].dec, *(random[i].n_ra), *(random[i].n_dec), angdist
		ind3 = where(angdist lt ang500[i], ind3count)
		random[i].n500 = ind3count

		; Find the number of neighbors within 500 Mpc and APPARENT magnitude brighter than counting mag

		distance = lumdist(random[i].z,/silent,/wmap7)
		app_mstar = (Mr_star+2) + 5*alog10(distance*1e6) - 5
		countingmag_app = sdss_lim < app_mstar

		;countingmag_app = sdss_lim
		;countingmag_app = app_mstar

		ind3_cmagapp = where((angdist lt ang500[i]) and (*(random[i].n_mag_r) lt countingmag_app), cmagcount_app)
		random[i].cmag_app = countingmag_app

			; Set the counting magnitude to ABSOLUTE, rather than apparent magnitude

			sdss_lim_abs = sdss_lim - 5*alog10(distance*1e6) + 5
			abs_mstar = Mr_star + 2
			countingmag_abs = sdss_lim_abs < abs_mstar

			;countingmag_abs = sdss_lim_abs
			;countingmag_abs = abs_mstar

			abs_rmag = *(random[i].n_absmag_r)
			ind3_cmagabs = where((angdist lt ang500[i]) and $
				(abs_rmag lt countingmag_abs), cmagcount_abs)

		n500_cmag = cmagcount_abs
		cmag = countingmag_abs

		random[i].n500_cmag = n500_cmag
		random[i].cmag = cmag

		; Northern control field

		ind_north = where(random[i].rname eq file_rname_north,icount_north)
		if icount_north gt 0 then begin
			icn_searchid = searchid_north[ind_north[0]]
			icn_ra = con_ra_north[ind_north[0]]
			icn_dec = con_dec_north[ind_north[0]]
			ind_ncon_north = where(ncon_searchid_north eq icn_searchid, magcount_north)

			if magcount_north gt 0 then begin
				ra_temp = ncon_ra_north[ind_ncon_north]
				dec_temp = ncon_dec_north[ind_ncon_north]
				absmagr_temp = ncon_absmagR_north[ind_ncon_north]

				gcirc, 2, icn_ra, icn_dec, ra_temp, dec_temp, con_angdist

				nind = where((absmagR_temp lt random[i].cmag) and (absmagR_temp gt -99) and (con_angdist lt ang500[i]), nc_north)
				;if nc_north gt 0 then *(ptrarr_bg_mag[i]) = absmagr_temp[nind]
				bgtemp_north = nc_north
			endif
		endif else begin
			nc_north = 0
			bgtemp_north = 0
		endelse

		; Southern control field

		ind_south = where(random[i].rname eq file_rname_south,icount_south)
		if icount_south gt 0 then begin
			icn_searchid = searchid_south[ind_south[0]]
			icn_ra = con_ra_south[ind_south[0]]
			icn_dec = con_dec_south[ind_south[0]]
			ind_ncon_south = where(ncon_searchid_south eq icn_searchid, magcount_south)

			if magcount_south gt 0 then begin
				ra_temp = ncon_ra_south[ind_ncon_south]
				dec_temp = ncon_dec_south[ind_ncon_south]
				absmagr_temp = ncon_absmagR_south[ind_ncon_south]

				gcirc, 2, icn_ra, icn_dec, ra_temp, dec_temp, con_angdist

				nind = where((absmagR_temp lt random[i].cmag) and (absmagR_temp gt -99) and (con_angdist lt ang500[i]), nc_south)
				;if nc_south gt 0 then *(ptrarr_bg_mag[i]) = absmagr_temp[nind]
				bgtemp_south = nc_south
			endif
		endif else begin
			nc_south = 0
			bgtemp_south = 0
		endelse

		if nc_north gt 0 then begin
			if nc_south gt 0 then random[i].bg = mean([bgtemp_north,bgtemp_south]) else random[i].bg = bgtemp_north
		endif else random[i].bg = bgtemp_south

	endfor

	;random.bg_mag = ptrarr_bg_mag
		
	save, filename=savdir+'random_structure.sav', random

	if ~keyword_set(quiet) then begin
		print,''
		print,'SDSS random selection data structure written to '+savdir+'random_structure.sav'
		print,''
	endif

	if keyword_set(stop) then stop

end
