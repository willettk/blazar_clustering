;+
; NAME:
;       
;	BZ_STRMAKE
;
; PURPOSE:
;
;	Create structure with SDSS, BZCAT information on known blazar properties
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
;	
;
; REVISION HISTORY
;       Written by K. Willett                Nov 11
;-

; Begin program

pro bz_strmake, stop=stop, readbfiles=readbfiles

	heap_gc		; Clean up heap variables

	csvdir = '~/Astronomy/Research/blazars/csv/'
	savdir = '~/Astronomy/Research/blazars/sav/'
	fitsdir = '~/Astronomy/Research/blazars/fits/'

	spindexfile = csvdir+'bzcat_aro.csv'

	bgbfile_3 = csvdir+'ab_neighbors3_willettk.csv'
	sdssfile_3 = csvdir+'ab_3arcmin_mags_willettk.csv'
	bgbfile_10 = csvdir+'ab_neighbors10_willettk.csv'
	sdssfile_10 = csvdir+'ab_10arcmin_mags_willettk.csv'
	bgbfile_all = csvdir+'ab_neighbors_all.csv'
	sdssfile_all = csvdir+'ab_mags_all.csv'

	spawn,'cat '+bgbfile_3+' '+bgbfile_10+' > '+bgbfile_all
	spawn,'cat '+sdssfile_3+' '+sdssfile_10+' > '+sdssfile_all

	if keyword_set(readbfiles) then begin

		readcol, bgbfile_all, $
			bname,ra,dec,z,btype,searchid,match_objid, $
			format='a,f,f,f,a,i,a', $
			skipline=1
		bbadind=where(strmid(match_objid,19,4) eq 'name')
		match_objid[bbadind] = strmid(match_objid[bbadind],0,19)

		readcol, sdssfile_all, $
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

		readcol, spindexfile, $
			sp_no, sp_name, sp_ra, sp_dec, $
			sp_z, sp_mag_r, sp_btype, sp_radioflux, sp_xrayflux, sp_fermiflux, $
			sp_aro, sp_arx, sp_aox, $
			delimiter=',' , $
			format='i,a,a,a,f,f,a,f,f,f,f,f,f', $
			skipline=1
			
		; Kluge to ensure unique search IDs

		nlines_bgb3 = file_lines(bgbfile_3)
		nlines_sdss3 = file_lines(sdssfile_3)

		searchid[0:nlines_bgb3-2] += 5000
		n_searchid[0:nlines_sdss3-2] += 5000

		; Load information for the control fields

		; North one degree from blazar
		file3arcmin_neighbors = csvdir+'ab_con_neighbors3_willettk.csv'
		file3arcmin_mags = csvdir+'ab_con_3arcmin_mags_willettk.csv'
		file10arcmin_neighbors = csvdir+'ab_con_neighbors10_willettk.csv'
		file10arcmin_mags = csvdir+'ab_con_10arcmin_mags_willettk.csv'

		; South one degree from blazar
		file3darcmin_neighbors = csvdir+'ab_con_d_neighbors3_willettk.csv'
		file3darcmin_mags = csvdir+'ab_con_d_3arcmin_mags_willettk.csv'
		file10darcmin_neighbors = csvdir+'ab_con_d_neighbors10_willettk.csv'
		file10darcmin_mags = csvdir+'ab_con_d_10arcmin_mags_willettk.csv'

		readcol, file3arcmin_neighbors, $
			c3_file_bname,c3_ra,c3_dec,c3_z,c3_btype,c3_search_id,c3_matched_id, $
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file3arcmin_mags, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_search_id,cn3_objid, $
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file10arcmin_neighbors, $
			c10_file_bname,c10_ra,c10_dec,c10_z,c10_btype,c10_search_id,c10_matched_id,$
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file10arcmin_mags, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_search_id,cn10_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file3darcmin_neighbors, $
			c3d_file_bname,c3d_ra,c3d_dec,c3d_z,c3d_btype,c3d_search_id,c3d_matched_id, $
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file3darcmin_mags, $
			cn3d_u,cn3d_g,cn3d_r,cn3d_i,cn3d_z,cn3d_err_u,cn3d_err_g,cn3d_err_r,cn3d_err_i,cn3d_err_z,$
			cn3d_ra,cn3d_dec,cn3d_type,cn3d_redshift,cn3d_redshift_err,$
			cn3d_absmagU,cn3d_absmagG,cn3d_absmagR,cn3d_absmagI,cn3d_absmagZ,$
			cn3d_kcorrU,cn3d_kcorrG,cn3d_kcorrR,cn3d_kcorrI,cn3d_kcorrZ,$
			cn3d_search_id,cn3d_objid, $
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file10darcmin_neighbors, $
			c10d_file_bname,c10d_ra,c10d_dec,c10d_z,c10d_btype,c10d_search_id,c10d_matched_id,$
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file10darcmin_mags, $
			cn10d_u,cn10d_g,cn10d_r,cn10d_i,cn10d_z,cn10d_err_u,cn10d_err_g,cn10d_err_r,cn10d_err_i,cn10d_err_z,$
			cn10d_ra,cn10d_dec,cn10d_type,cn10d_redshift,cn10d_redshift_err,$
			cn10d_absmagU,cn10d_absmagG,cn10d_absmagR,cn10d_absmagI,cn10d_absmagZ,$
			cn10d_kcorrU,cn10d_kcorrG,cn10d_kcorrR,cn10d_kcorrI,cn10d_kcorrZ,$
			cn10d_search_id,cn10d_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		save, filename=savdir+'blazar_controlfields.sav', $
			c3_file_bname,c3_ra,c3_dec,c3_z,c3_btype,c3_search_id,c3_matched_id, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_search_id,cn3_objid, $
			c10_file_bname,c10_ra,c10_dec,c10_z,c10_btype,c10_search_id,c10_matched_id, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_search_id,cn10_objid, $
			c3d_file_bname,c3d_ra,c3d_dec,c3d_z,c3d_btype,c3d_search_id,c3d_matched_id, $
			cn3d_u,cn3d_g,cn3d_r,cn3d_i,cn3d_z,cn3d_err_u,cn3d_err_g,cn3d_err_r,cn3d_err_i,cn3d_err_z,$
			cn3d_ra,cn3d_dec,cn3d_type,cn3d_redshift,cn3d_redshift_err,$
			cn3d_absmagU,cn3d_absmagG,cn3d_absmagR,cn3d_absmagI,cn3d_absmagZ,$
			cn3d_kcorrU,cn3d_kcorrG,cn3d_kcorrR,cn3d_kcorrI,cn3d_kcorrZ,$
			cn3d_search_id,cn3d_objid, $
			c10d_file_bname,c10d_ra,c10d_dec,c10d_z,c10d_btype,c10d_search_id,c10d_matched_id, $
			cn10d_u,cn10d_g,cn10d_r,cn10d_i,cn10d_z,cn10d_err_u,cn10d_err_g,cn10d_err_r,cn10d_err_i,cn10d_err_z,$
			cn10d_ra,cn10d_dec,cn10d_type,cn10d_redshift,cn10d_redshift_err,$
			cn10d_absmagU,cn10d_absmagG,cn10d_absmagR,cn10d_absmagI,cn10d_absmagZ,$
			cn10d_kcorrU,cn10d_kcorrG,cn10d_kcorrR,cn10d_kcorrI,cn10d_kcorrZ,$
			cn10d_search_id,cn10d_objid
	
		save, filename=savdir+'ab_neighbors_all.sav', $
			bname, ra, dec, z, btype, searchid, match_objid, $
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
			n_objid, $
			sp_no, sp_name , sp_ra, sp_dec, $
			sp_z, sp_mag_r, sp_btype, sp_radioflux, sp_xrayflux, sp_fermiflux, $
			sp_aro, sp_arx, sp_aox

	endif else begin
		restore, savdir+'blazar_controlfields.sav'
		restore, savdir+'ab_neighbors_all.sav'
	endelse

	b = {bname:'', btype:'', $
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
		fieldsize:0.,$
		cmag:0.,$
		cmag_app:0.,$
		n500:0, $
		n500_cmag:0, $
		flux_radio:0., $
		flux_xray:0., $
		sp_mag_r:0., $
		spindex_ro:0., $
		spindex_rx:0., $
		spindex_ox:0.}

	; Match the SDSS data to its blazar parent ID

	unique_searchid1 = unique(searchid,/sort)
	unique_nsearchid1 = unique(n_searchid,/sort)

	; Remove blazars that had no galaxies within the search radii (kluge for now; should retain info)

	nogals = setdifference(unique_searchid1, unique_nsearchid1)
	nogalsind = [0]
	for i=0, n_elements(nogals)-1 do nogalsind = [nogalsind,where(searchid eq nogals[i])]
	nogalsind = nogalsind[1:n_elements(nogalsind)-1]
	galsind = setdifference(lindgen(n_elements(bname)),nogalsind)

		bname = bname[galsind]
		ra = ra[galsind]
		dec = dec[galsind]
		z = z[galsind]
		btype = btype[galsind]
		match_objid = match_objid[galsind]
		searchid = searchid[galsind]


	; Find unique parent blazars

	name_uniq = unique(bname,/sort)
	nb = n_elements(name_uniq)
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
		endif else print,'No match for search ID ',unique_searchid2[i], '#',i
	endfor 

	bz = replicate(b,nb)

	bz.n_mag_u = ptrarr_mag_u
	bz.n_mag_g = ptrarr_mag_g
	bz.n_mag_r = ptrarr_mag_r
	bz.n_mag_i = ptrarr_mag_i
	bz.n_mag_z = ptrarr_mag_z
	bz.n_mag_err_u = ptrarr_mag_err_u
	bz.n_mag_err_g = ptrarr_mag_err_g
	bz.n_mag_err_r = ptrarr_mag_err_r
	bz.n_mag_err_i = ptrarr_mag_err_i
	bz.n_mag_err_z = ptrarr_mag_err_z
	bz.n_absmag_u = ptrarr_absmag_u
	bz.n_absmag_g = ptrarr_absmag_g
	bz.n_absmag_r = ptrarr_absmag_r
	bz.n_absmag_i = ptrarr_absmag_i
	bz.n_absmag_z = ptrarr_absmag_z
	bz.n_redshift = ptrarr_redshift
	bz.n_redshift_err = ptrarr_redshift_err
	bz.n_ra = ptrarr_ra
	bz.n_dec = ptrarr_dec

	blazarind = uniq(searchid,sort(searchid))

	bz.bname=bname[blazarind]
	bz.btype=btype[blazarind]
	bz.ra=ra[blazarind]
	bz.dec=dec[blazarind]
	bz.z=z[blazarind]
	bz.searchid=searchid[blazarind]

	ptrarr_matchedid = ptrarr(nb,/allocate_heap)
	for i=0,nb-1 do begin
		*(ptrarr_matchedid[i]) = match_objid[where(searchid eq searchid[blazarind[i]])]
	endfor
	bz.match_objid=ptrarr_matchedid

	; Insert the multi-wavelength data from BzCat

	match, strtrim(sp_name,2), strtrim(bz.bname,2), spind, bspind, count=spcount
	if spcount gt 0 then begin
		bz[bspind].flux_radio = sp_radioflux[spind]
		bz[bspind].flux_xray = sp_xrayflux[spind]
		bz[bspind].spindex_ox = sp_aox[spind]
		bz[bspind].spindex_rx = sp_arx[spind]
		bz[bspind].spindex_ro = sp_aro[spind]
		bz[bspind].sp_mag_r = sp_mag_r[spind]
	endif

	; For each galaxy, find the 500 kpc angular distance [arcsec]

	ang500 = zang(500., bz.z,/silent,/wmap7)
	bz.fieldsize = ang500
	Mr_star = -21.1837		; Blanton et al. (2003)
	sdss_lim = 22.2			; SDSS apparent limiting magnitude

	;ptrarr_bg_mag = ptrarr(nb,/allocate_heap)

	; Background galaxy array generation

	file_bname_north = [c3_file_bname,c10_file_bname]
	search_id_north = [c3_search_id,c10_search_id+5000]
	con_ra_north = [c3_ra, c10_ra]
	con_dec_north = [c3_dec, c10_dec]

	ncon_search_id_north = [cn3_search_id,cn10_search_id+5000]
	ncon_ra_north = [cn3_ra, cn10_ra]
	ncon_dec_north = [cn3_dec, cn10_dec]
	ncon_absmagR_north = [cn3_absmagR, cn10_absmagR]

	file_bname_south = [c3d_file_bname,c10d_file_bname]
	search_id_south = [c3d_search_id,c10d_search_id+5000]
	con_ra_south = [c3d_ra, c10d_ra]
	con_dec_south = [c3d_dec, c10d_dec]

	ncon_search_id_south = [cn3d_search_id,cn10d_search_id+5000]
	ncon_ra_south = [cn3d_ra, cn10d_ra]
	ncon_dec_south = [cn3d_dec, cn10d_dec]
	ncon_absmagR_south = [cn3d_absmagR, cn10d_absmagR]

	for i=0,nb-1 do begin

		; Find number of neighbors within projected 500 Mpc

		gcirc, 2, bz[i].ra, bz[i].dec, *(bz[i].n_ra), *(bz[i].n_dec), angdist
		ind3 = where(angdist lt ang500[i], ind3count)
		bz[i].n500 = ind3count

		; Find the number of neighbors within 500 Mpc and APPARENT magnitude brighter than counting mag

		distance = lumdist(bz[i].z,/silent,/wmap7)
		app_mstar = (Mr_star+2) + 5*alog10(distance*1e6) - 5
		countingmag_app = sdss_lim < app_mstar

		;countingmag_app = sdss_lim
		;countingmag_app = app_mstar

		ind3_cmagapp = where((angdist lt ang500[i]) and (*(bz[i].n_mag_r) lt countingmag_app), cmagcount_app)
		bz[i].cmag_app = countingmag_app

			; Set the counting magnitude to ABSOLUTE, rather than apparent magnitude

			sdss_lim_abs = sdss_lim - 5*alog10(distance*1e6) + 5
			abs_mstar = Mr_star + 2
			countingmag_abs = sdss_lim_abs < abs_mstar

			;countingmag_abs = sdss_lim_abs
			;countingmag_abs = abs_mstar

			abs_rmag = *(bz[i].n_absmag_r)
			ind3_cmagabs = where((angdist lt ang500[i]) and $
				(abs_rmag lt countingmag_abs), cmagcount_abs)

		n500_cmag = cmagcount_abs
		cmag = countingmag_abs

		bz[i].n500_cmag = n500_cmag
		bz[i].cmag = cmag

		; Northern control field

		ind_north = where(bz[i].bname eq file_bname_north,icount_north)
		if icount_north gt 0 then begin
			icn_search_id = search_id_north[ind_north[0]]
			icn_ra = con_ra_north[ind_north[0]]
			icn_dec = con_dec_north[ind_north[0]]
			ind_ncon_north = where(ncon_search_id_north eq icn_search_id, magcount_north)

			if magcount_north gt 0 then begin
				ra_temp = ncon_ra_north[ind_ncon_north]
				dec_temp = ncon_dec_north[ind_ncon_north]
				absmagr_temp = ncon_absmagR_north[ind_ncon_north]

				gcirc, 2, icn_ra, icn_dec, ra_temp, dec_temp, con_angdist

				nind = where((absmagR_temp lt bz[i].cmag) and (absmagR_temp gt -99) and (con_angdist lt ang500[i]), nc_north)
				;if nc_north gt 0 then *(ptrarr_bg_mag[i]) = absmagr_temp[nind]
				bgtemp_north = nc_north
			endif
		endif else begin
			nc_north = 0
			bgtemp_north = 0
		endelse

		; Southern control field

		ind_south = where(bz[i].bname eq file_bname_south,icount_south)
		if icount_south gt 0 then begin
			icn_search_id = search_id_south[ind_south[0]]
			icn_ra = con_ra_south[ind_south[0]]
			icn_dec = con_dec_south[ind_south[0]]
			ind_ncon_south = where(ncon_search_id_south eq icn_search_id, magcount_south)

			if magcount_south gt 0 then begin
				ra_temp = ncon_ra_south[ind_ncon_south]
				dec_temp = ncon_dec_south[ind_ncon_south]
				absmagr_temp = ncon_absmagR_south[ind_ncon_south]

				gcirc, 2, icn_ra, icn_dec, ra_temp, dec_temp, con_angdist

				nind = where((absmagR_temp lt bz[i].cmag) and (absmagR_temp gt -99) and (con_angdist lt ang500[i]), nc_south)
				;if nc_south gt 0 then *(ptrarr_bg_mag[i]) = absmagr_temp[nind]
				bgtemp_south = nc_south
			endif
		endif else begin
			nc_south = 0
			bgtemp_south = 0
		endelse

		if nc_north gt 0 then begin
			if nc_south gt 0 then bz[i].bg = mean([bgtemp_north,bgtemp_south]) else bz[i].bg = bgtemp_north
		endif else bz[i].bg = bgtemp_south

	endfor

	;bz.bg_mag = ptrarr_bg_mag
		
	; Blazar type indices

	bllac = where(bz.btype eq 'BLLac')
	fsrq = where(bz.btype eq 'FSRQ')
	candidate = where(bz.btype eq 'BLLac_candidate')
	uncertain = where(bz.btype eq 'blazar_uncertain')

	save, filename=savdir+'ab_structure.sav', bz, $
		bllac, fsrq, candidate, uncertain

	mwrfits, fitsdir+'bz_strmake.fits', bz

	print,''
	print,'Blazar structure written to '+savdir+'ab_structure.sav'
	print,''

	if keyword_set(stop) then stop

end
