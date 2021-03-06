
;+
; NAME:
;       
;	GZ_STRMAKE
;
; PURPOSE:
;
;	Create structure with SDSS, BZCAT information on SDSS-selected radio galaxies
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
;	IDL> gz_strmake
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Apr 12
;-

; Begin program

pro gz_strmake, stop=stop, readbfiles=readbfiles

	; ELLIPTICAL GALAXIES

	heap_gc		; Clean up heap variables

	csvdir = '~/Astronomy/Research/blazars/csv/'
	savdir = '~/Astronomy/Research/blazars/sav/'
	gzdir = '~/Astronomy/Research/blazars/gz/'

	el_file_3 = gzdir+'el_zns_modz_3arcmin_n_willettk.csv'
	el_sdssfile_3 = gzdir+'el_zns_modz_3arcmin_sdss_willettk.csv'
	el_file_10 = gzdir+'el_zns_modz_10arcmin_n_willettk.csv'
	el_sdssfile_10 = gzdir+'el_zns_modz_10arcmin_sdss_willettk.csv'

	el_neighbors_all = gzdir+'el_neighbors_all.csv'
	el_sdss_all = gzdir+'el_sdss_all.csv'

	spawn,'cat '+el_file_3+' '+el_file_10+' > '+el_neighbors_all
	spawn,'cat '+el_sdssfile_3+' '+el_sdssfile_10+' > '+el_sdss_all

	if keyword_set(readbfiles) then begin

		readcol, el_neighbors_all, $
			gname,ra,dec,z,p_el,searchid,match_objid, $
			format='a,f,f,f,f,i,a', $
			skipline=1
		bbadind=where(strmid(match_objid,19,4) eq 'name')
		match_objid[bbadind] = strmid(match_objid[bbadind],0,19)

		readcol, el_sdss_all, $
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

		nlines_bgb3 = file_lines(el_file_3)
		nlines_sdss3 = file_lines(el_sdssfile_3)

		searchid[0:nlines_bgb3-2] += 5000
		n_searchid[0:nlines_sdss3-2] += 5000

	;	searchid[where(searchid eq 5638)] = 5635

		; Load information for the control fields

		; North one degree from blazar
		file_connorth_3arcmin_neighbors = gzdir+'el_zns_modz_3arcmin_connorth_n_willettk.csv'
		file_connorth_3arcmin_sdss = gzdir+'el_zns_modz_3arcmin_connorth_sdss_willettk.csv'
		file_connorth_10arcmin_neighbors = gzdir+'el_zns_modz_10arcmin_connorth_n_willettk.csv'
		file_connorth_10arcmin_sdss = gzdir+'el_zns_modz_10arcmin_connorth_sdss_willettk.csv'

		; South one degree from blazar
		file_consouth_3arcmin_neighbors = gzdir+'el_zns_modz_3arcmin_consouth_n_willettk.csv'
		file_consouth_3arcmin_sdss = gzdir+'el_zns_modz_3arcmin_consouth_sdss_willettk.csv'
		file_consouth_10arcmin_neighbors = gzdir+'el_zns_modz_10arcmin_consouth_n_willettk.csv'
		file_consouth_10arcmin_sdss = gzdir+'el_zns_modz_10arcmin_consouth_sdss_willettk.csv'

		readcol, file_connorth_3arcmin_neighbors, $
			c3_file_gname,c3_ra,c3_dec,c3_z,c3_p_el,c3_search_id,c3_matched_id, $
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_connorth_3arcmin_sdss, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_search_id,cn3_objid, $
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_connorth_10arcmin_neighbors, $
			c10_file_gname,c10_ra,c10_dec,c10_z,c10_p_el,c10_search_id,c10_matched_id,$
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_connorth_10arcmin_sdss, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_search_id,cn10_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_consouth_3arcmin_neighbors, $
			c3d_file_gname,c3d_ra,c3d_dec,c3d_z,c3d_p_el,c3d_search_id,c3d_matched_id, $
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_consouth_3arcmin_sdss, $
			cn3d_u,cn3d_g,cn3d_r,cn3d_i,cn3d_z,cn3d_err_u,cn3d_err_g,cn3d_err_r,cn3d_err_i,cn3d_err_z,$
			cn3d_ra,cn3d_dec,cn3d_type,cn3d_redshift,cn3d_redshift_err,$
			cn3d_absmagU,cn3d_absmagG,cn3d_absmagR,cn3d_absmagI,cn3d_absmagZ,$
			cn3d_kcorrU,cn3d_kcorrG,cn3d_kcorrR,cn3d_kcorrI,cn3d_kcorrZ,$
			cn3d_search_id,cn3d_objid, $
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_consouth_10arcmin_neighbors, $
			c10d_file_gname,c10d_ra,c10d_dec,c10d_z,c10d_p_el,c10d_search_id,c10d_matched_id,$
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_consouth_10arcmin_sdss, $
			cn10d_u,cn10d_g,cn10d_r,cn10d_i,cn10d_z,cn10d_err_u,cn10d_err_g,cn10d_err_r,cn10d_err_i,cn10d_err_z,$
			cn10d_ra,cn10d_dec,cn10d_type,cn10d_redshift,cn10d_redshift_err,$
			cn10d_absmagU,cn10d_absmagG,cn10d_absmagR,cn10d_absmagI,cn10d_absmagZ,$
			cn10d_kcorrU,cn10d_kcorrG,cn10d_kcorrR,cn10d_kcorrI,cn10d_kcorrZ,$
			cn10d_search_id,cn10d_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		save, filename=savdir+'el_controlfields.sav', $
			c3_file_gname,c3_ra,c3_dec,c3_z,c3_p_el,c3_search_id,c3_matched_id, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_search_id,cn3_objid, $
			c10_file_gname,c10_ra,c10_dec,c10_z,c10_p_el,c10_search_id,c10_matched_id, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_search_id,cn10_objid, $
			c3d_file_gname,c3d_ra,c3d_dec,c3d_z,c3d_p_el,c3d_search_id,c3d_matched_id, $
			cn3d_u,cn3d_g,cn3d_r,cn3d_i,cn3d_z,cn3d_err_u,cn3d_err_g,cn3d_err_r,cn3d_err_i,cn3d_err_z,$
			cn3d_ra,cn3d_dec,cn3d_type,cn3d_redshift,cn3d_redshift_err,$
			cn3d_absmagU,cn3d_absmagG,cn3d_absmagR,cn3d_absmagI,cn3d_absmagZ,$
			cn3d_kcorrU,cn3d_kcorrG,cn3d_kcorrR,cn3d_kcorrI,cn3d_kcorrZ,$
			cn3d_search_id,cn3d_objid, $
			c10d_file_gname,c10d_ra,c10d_dec,c10d_z,c10d_p_el,c10d_search_id,c10d_matched_id, $
			cn10d_u,cn10d_g,cn10d_r,cn10d_i,cn10d_z,cn10d_err_u,cn10d_err_g,cn10d_err_r,cn10d_err_i,cn10d_err_z,$
			cn10d_ra,cn10d_dec,cn10d_type,cn10d_redshift,cn10d_redshift_err,$
			cn10d_absmagU,cn10d_absmagG,cn10d_absmagR,cn10d_absmagI,cn10d_absmagZ,$
			cn10d_kcorrU,cn10d_kcorrG,cn10d_kcorrR,cn10d_kcorrI,cn10d_kcorrZ,$
			cn10d_search_id,cn10d_objid
	
		save, filename=savdir+'el_neighbors_all.sav', $
			gname, ra, dec, z, p_el, searchid, match_objid, $
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

	endif else begin
		restore, savdir+'el_controlfields.sav'
		restore, savdir+'el_neighbors_all.sav'
	endelse

	b = {gname:'', p_el:'', $
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
		n500_cmag:0}

	; Match the SDSS data to its galaxy parent ID

	unique_searchid1 = unique(searchid,/sort)
	unique_nsearchid1 = unique(n_searchid,/sort)

	; Remove galaxies that had no galaxies within the search radii (kluge for now; should retain info)

	nogals = setdifference(unique_searchid1, unique_nsearchid1)
	nogalsind = [0]
	for i=0, n_elements(nogals)-1 do nogalsind = [nogalsind,where(searchid eq nogals[i])]
	nogalsind = nogalsind[1:n_elements(nogalsind)-1]
	galsind = setdifference(lindgen(n_elements(gname)),nogalsind)

		gname = gname[galsind]
		ra = ra[galsind]
		dec = dec[galsind]
		z = z[galsind]
		p_el = p_el[galsind]
		match_objid = match_objid[galsind]
		searchid = searchid[galsind]


	; Find unique parent galaxies

	name_uniq = unique(gname,/sort)
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

	el = replicate(b,nb)

	el.n_mag_u = ptrarr_mag_u
	el.n_mag_g = ptrarr_mag_g
	el.n_mag_r = ptrarr_mag_r
	el.n_mag_i = ptrarr_mag_i
	el.n_mag_z = ptrarr_mag_z
	el.n_mag_err_u = ptrarr_mag_err_u
	el.n_mag_err_g = ptrarr_mag_err_g
	el.n_mag_err_r = ptrarr_mag_err_r
	el.n_mag_err_i = ptrarr_mag_err_i
	el.n_mag_err_z = ptrarr_mag_err_z
	el.n_absmag_u = ptrarr_absmag_u
	el.n_absmag_g = ptrarr_absmag_g
	el.n_absmag_r = ptrarr_absmag_r
	el.n_absmag_i = ptrarr_absmag_i
	el.n_absmag_z = ptrarr_absmag_z
	el.n_redshift = ptrarr_redshift
	el.n_redshift_err = ptrarr_redshift_err
	el.n_ra = ptrarr_ra
	el.n_dec = ptrarr_dec

	ugalind = uniq(searchid,sort(searchid))

	el.gname=gname[ugalind]
	el.p_el=p_el[ugalind]
	el.ra=ra[ugalind]
	el.dec=dec[ugalind]
	el.z=z[ugalind]
	el.searchid=searchid[ugalind]

	ptrarr_matchedid = ptrarr(nb,/allocate_heap)
	for i=0,nb-1 do begin
		*(ptrarr_matchedid[i]) = match_objid[where(searchid eq searchid[ugalind[i]])]
	endfor
	el.match_objid=ptrarr_matchedid

	; Insert the multi-wavelength data from BzCat

	;match, strtrim(sp_name,2), strtrim(bz.gname,2), spind, bspind, count=spcount
	;if spcount gt 0 then begin
	;	bz[bspind].flux_radio = sp_radioflux[spind]
	;	bz[bspind].flux_xray = sp_xrayflux[spind]
	;	bz[bspind].spindex_ox = sp_aox[spind]
	;	bz[bspind].spindex_rx = sp_arx[spind]
	;	bz[bspind].spindex_ro = sp_aro[spind]
	;	bz[bspind].sp_mag_r = sp_mag_r[spind]
	;endif

	; For each galaxy, find the 500 kpc angular distance [arcsec]

	ang500 = zang(500., el.z,/silent,/wmap7)
	el.fieldsize = ang500
	Mr_star = -21.1837		; Blanton et al. (2003)
	sdss_lim = 22.2			; SDSS apparent limiting magnitude

	;ptrarr_bg_mag = ptrarr(nb,/allocate_heap)

	; Background galaxy array generation

	file_gname_north = [c3_file_gname,c10_file_gname]
	search_id_north = [c3_search_id,c10_search_id+5000]
	con_ra_north = [c3_ra, c10_ra]
	con_dec_north = [c3_dec, c10_dec]

	ncon_search_id_north = [cn3_search_id,cn10_search_id+5000]
	ncon_ra_north = [cn3_ra, cn10_ra]
	ncon_dec_north = [cn3_dec, cn10_dec]
	ncon_absmagR_north = [cn3_absmagR, cn10_absmagR]

	file_gname_south = [c3d_file_gname,c10d_file_gname]
	search_id_south = [c3d_search_id,c10d_search_id+5000]
	con_ra_south = [c3d_ra, c10d_ra]
	con_dec_south = [c3d_dec, c10d_dec]

	ncon_search_id_south = [cn3d_search_id,cn10d_search_id+5000]
	ncon_ra_south = [cn3d_ra, cn10d_ra]
	ncon_dec_south = [cn3d_dec, cn10d_dec]
	ncon_absmagR_south = [cn3d_absmagR, cn10d_absmagR]

	for i=0,nb-1 do begin

		; Find number of neighbors within projected 500 Mpc

		gcirc, 2, el[i].ra, el[i].dec, *(el[i].n_ra), *(el[i].n_dec), angdist
		ind3 = where(angdist lt ang500[i], ind3count)
		el[i].n500 = ind3count

		; Find the number of neighbors within 500 Mpc and APPARENT magnitude brighter than counting mag

		distance = lumdist(el[i].z,/silent,/wmap7)
		app_mstar = (Mr_star+2) + 5*alog10(distance*1e6) - 5
		countingmag_app = sdss_lim < app_mstar

		;countingmag_app = sdss_lim
		;countingmag_app = app_mstar

		ind3_cmagapp = where((angdist lt ang500[i]) and (*(el[i].n_mag_r) lt countingmag_app), cmagcount_app)
		el[i].cmag_app = countingmag_app

			; Set the counting magnitude to ABSOLUTE, rather than apparent magnitude

			sdss_lim_abs = sdss_lim - 5*alog10(distance*1e6) + 5
			abs_mstar = Mr_star + 2
			countingmag_abs = sdss_lim_abs < abs_mstar

			;countingmag_abs = sdss_lim_abs
			;countingmag_abs = abs_mstar

			abs_rmag = *(el[i].n_absmag_r)
			ind3_cmagabs = where((angdist lt ang500[i]) and $
				(abs_rmag lt countingmag_abs), cmagcount_abs)

		n500_cmag = cmagcount_abs
		cmag = countingmag_abs

		el[i].n500_cmag = n500_cmag
		el[i].cmag = cmag

		; Northern control field

		ind_north = where(el[i].gname eq file_gname_north,icount_north)
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

				nind = where((absmagR_temp lt el[i].cmag) and (absmagR_temp gt -99) and (con_angdist lt ang500[i]), nc_north)
				;if nc_north gt 0 then *(ptrarr_bg_mag[i]) = absmagr_temp[nind]
				bgtemp_north = nc_north
			endif
		endif else begin
			nc_north = 0
			bgtemp_north = 0
		endelse

		; Southern control field

		ind_south = where(el[i].gname eq file_gname_south,icount_south)
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

				nind = where((absmagR_temp lt el[i].cmag) and (absmagR_temp gt -99) and (con_angdist lt ang500[i]), nc_south)
				;if nc_south gt 0 then *(ptrarr_bg_mag[i]) = absmagr_temp[nind]
				bgtemp_south = nc_south
			endif
		endif else begin
			nc_south = 0
			bgtemp_south = 0
		endelse

		if nc_north gt 0 then begin
			if nc_south gt 0 then el[i].bg = mean([bgtemp_north,bgtemp_south]) else el[i].bg = bgtemp_north
		endif else el[i].bg = bgtemp_south

	endfor

	save, filename=savdir+'el_modz_structure.sav', el

	print,''
	print,'Elliptical (modz) galaxy structure written to '+savdir+'el_modz_structure.sav'
	print,''

	; SPIRAL GALAXIES

	heap_gc		; Clean up heap variables

	sp_file_3 = gzdir+'sp_zns_modz_3arcmin_n_willettk.csv'
	sp_sdssfile_3 = gzdir+'sp_zns_modz_3arcmin_sdss_willettk.csv'
	sp_file_10 = gzdir+'sp_zns_modz_10arcmin_n_willettk.csv'
	sp_sdssfile_10 = gzdir+'sp_zns_modz_10arcmin_sdss_willettk.csv'

	sp_neighbors_all = gzdir+'sp_neighbors_all.csv'
	sp_sdss_all = gzdir+'sp_sdss_all.csv'

	spawn,'cat '+sp_file_3+' '+sp_file_10+' > '+sp_neighbors_all
	spawn,'cat '+sp_sdssfile_3+' '+sp_sdssfile_10+' > '+sp_sdss_all

	if keyword_set(readbfiles) then begin

		readcol, sp_neighbors_all, $
			gname,ra,dec,z,p_sp,searchid,match_objid, $
			format='a,f,f,f,f,i,a', $
			skipline=1
		bbadind=where(strmid(match_objid,19,4) eq 'name')
		match_objid[bbadind] = strmid(match_objid[bbadind],0,19)

		readcol, sp_sdss_all, $
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

		nlines_bgb3 = file_lines(sp_file_3)
		nlines_sdss3 = file_lines(sp_sdssfile_3)

		searchid[0:nlines_bgb3-2] += 5000
		n_searchid[0:nlines_sdss3-2] += 5000

	;	searchid[where(searchid eq 5638)] = 5635

		; Load information for the control fields

		; North one degree from blazar
		file_connorth_3arcmin_neighbors = gzdir+'sp_zns_modz_3arcmin_connorth_n_willettk.csv'
		file_connorth_3arcmin_sdss = gzdir+'sp_zns_modz_3arcmin_connorth_sdss_willettk.csv'
		file_connorth_10arcmin_neighbors = gzdir+'sp_zns_modz_10arcmin_connorth_n_willettk.csv'
		file_connorth_10arcmin_sdss = gzdir+'sp_zns_modz_10arcmin_connorth_sdss_willettk.csv'

		; South one degree from blazar
		file_consouth_3arcmin_neighbors = gzdir+'sp_zns_modz_3arcmin_consouth_n_willettk.csv'
		file_consouth_3arcmin_sdss = gzdir+'sp_zns_modz_3arcmin_consouth_sdss_willettk.csv'
		file_consouth_10arcmin_neighbors = gzdir+'sp_zns_modz_10arcmin_consouth_n_willettk.csv'
		file_consouth_10arcmin_sdss = gzdir+'sp_zns_modz_10arcmin_consouth_sdss_willettk.csv'

		readcol, file_connorth_3arcmin_neighbors, $
			c3_file_gname,c3_ra,c3_dec,c3_z,c3_p_sp,c3_search_id,c3_matched_id, $
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_connorth_3arcmin_sdss, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_search_id,cn3_objid, $
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_connorth_10arcmin_neighbors, $
			c10_file_gname,c10_ra,c10_dec,c10_z,c10_p_sp,c10_search_id,c10_matched_id,$
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_connorth_10arcmin_sdss, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_search_id,cn10_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_consouth_3arcmin_neighbors, $
			c3d_file_gname,c3d_ra,c3d_dec,c3d_z,c3d_p_sp,c3d_search_id,c3d_matched_id, $
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_consouth_3arcmin_sdss, $
			cn3d_u,cn3d_g,cn3d_r,cn3d_i,cn3d_z,cn3d_err_u,cn3d_err_g,cn3d_err_r,cn3d_err_i,cn3d_err_z,$
			cn3d_ra,cn3d_dec,cn3d_type,cn3d_redshift,cn3d_redshift_err,$
			cn3d_absmagU,cn3d_absmagG,cn3d_absmagR,cn3d_absmagI,cn3d_absmagZ,$
			cn3d_kcorrU,cn3d_kcorrG,cn3d_kcorrR,cn3d_kcorrI,cn3d_kcorrZ,$
			cn3d_search_id,cn3d_objid, $
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		readcol, file_consouth_10arcmin_neighbors, $
			c10d_file_gname,c10d_ra,c10d_dec,c10d_z,c10d_p_sp,c10d_search_id,c10d_matched_id,$
			skip=1, $
			delimiter=',', $
			format='a,f,f,f,a,a,a'

		readcol, file_consouth_10arcmin_sdss, $
			cn10d_u,cn10d_g,cn10d_r,cn10d_i,cn10d_z,cn10d_err_u,cn10d_err_g,cn10d_err_r,cn10d_err_i,cn10d_err_z,$
			cn10d_ra,cn10d_dec,cn10d_type,cn10d_redshift,cn10d_redshift_err,$
			cn10d_absmagU,cn10d_absmagG,cn10d_absmagR,cn10d_absmagI,cn10d_absmagZ,$
			cn10d_kcorrU,cn10d_kcorrG,cn10d_kcorrR,cn10d_kcorrI,cn10d_kcorrZ,$
			cn10d_search_id,cn10d_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		save, filename=savdir+'sp_controlfields.sav', $
			c3_file_gname,c3_ra,c3_dec,c3_z,c3_p_sp,c3_search_id,c3_matched_id, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_search_id,cn3_objid, $
			c10_file_gname,c10_ra,c10_dec,c10_z,c10_p_sp,c10_search_id,c10_matched_id, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_search_id,cn10_objid, $
			c3d_file_gname,c3d_ra,c3d_dec,c3d_z,c3d_p_sp,c3d_search_id,c3d_matched_id, $
			cn3d_u,cn3d_g,cn3d_r,cn3d_i,cn3d_z,cn3d_err_u,cn3d_err_g,cn3d_err_r,cn3d_err_i,cn3d_err_z,$
			cn3d_ra,cn3d_dec,cn3d_type,cn3d_redshift,cn3d_redshift_err,$
			cn3d_absmagU,cn3d_absmagG,cn3d_absmagR,cn3d_absmagI,cn3d_absmagZ,$
			cn3d_kcorrU,cn3d_kcorrG,cn3d_kcorrR,cn3d_kcorrI,cn3d_kcorrZ,$
			cn3d_search_id,cn3d_objid, $
			c10d_file_gname,c10d_ra,c10d_dec,c10d_z,c10d_p_sp,c10d_search_id,c10d_matched_id, $
			cn10d_u,cn10d_g,cn10d_r,cn10d_i,cn10d_z,cn10d_err_u,cn10d_err_g,cn10d_err_r,cn10d_err_i,cn10d_err_z,$
			cn10d_ra,cn10d_dec,cn10d_type,cn10d_redshift,cn10d_redshift_err,$
			cn10d_absmagU,cn10d_absmagG,cn10d_absmagR,cn10d_absmagI,cn10d_absmagZ,$
			cn10d_kcorrU,cn10d_kcorrG,cn10d_kcorrR,cn10d_kcorrI,cn10d_kcorrZ,$
			cn10d_search_id,cn10d_objid
	
		save, filename=savdir+'sp_neighbors_all.sav', $
			gname, ra, dec, z, p_sp, searchid, match_objid, $
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

	endif else begin
		restore, savdir+'sp_controlfields.sav'
		restore, savdir+'sp_neighbors_all.sav'
	endelse

	b = {gname:'', p_sp:'', $
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
		n500_cmag:0}

	; Match the SDSS data to its galaxy parent ID

	unique_searchid1 = unique(searchid,/sort)
	unique_nsearchid1 = unique(n_searchid,/sort)

	; Remove galaxies that had no galaxies within the search radii (kluge for now; should retain info)

	nogals = setdifference(unique_searchid1, unique_nsearchid1)
	nogalsind = [0]
	for i=0, n_elements(nogals)-1 do nogalsind = [nogalsind,where(searchid eq nogals[i])]
	nogalsind = nogalsind[1:n_elements(nogalsind)-1]
	galsind = setdifference(lindgen(n_elements(gname)),nogalsind)

		gname = gname[galsind]
		ra = ra[galsind]
		dec = dec[galsind]
		z = z[galsind]
		p_sp = p_sp[galsind]
		match_objid = match_objid[galsind]
		searchid = searchid[galsind]


	; Find unique parent galaxies

	name_uniq = unique(gname,/sort)
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

	sp = replicate(b,nb)

	sp.n_mag_u = ptrarr_mag_u
	sp.n_mag_g = ptrarr_mag_g
	sp.n_mag_r = ptrarr_mag_r
	sp.n_mag_i = ptrarr_mag_i
	sp.n_mag_z = ptrarr_mag_z
	sp.n_mag_err_u = ptrarr_mag_err_u
	sp.n_mag_err_g = ptrarr_mag_err_g
	sp.n_mag_err_r = ptrarr_mag_err_r
	sp.n_mag_err_i = ptrarr_mag_err_i
	sp.n_mag_err_z = ptrarr_mag_err_z
	sp.n_absmag_u = ptrarr_absmag_u
	sp.n_absmag_g = ptrarr_absmag_g
	sp.n_absmag_r = ptrarr_absmag_r
	sp.n_absmag_i = ptrarr_absmag_i
	sp.n_absmag_z = ptrarr_absmag_z
	sp.n_redshift = ptrarr_redshift
	sp.n_redshift_err = ptrarr_redshift_err
	sp.n_ra = ptrarr_ra
	sp.n_dec = ptrarr_dec

	ugalind = uniq(searchid,sort(searchid))

	sp.gname=gname[ugalind]
	sp.p_sp=p_sp[ugalind]
	sp.ra=ra[ugalind]
	sp.dec=dec[ugalind]
	sp.z=z[ugalind]
	sp.searchid=searchid[ugalind]

	ptrarr_matchedid = ptrarr(nb,/allocate_heap)
	for i=0,nb-1 do begin
		*(ptrarr_matchedid[i]) = match_objid[where(searchid eq searchid[ugalind[i]])]
	endfor
	sp.match_objid=ptrarr_matchedid

	; Insert the multi-wavelength data from BzCat

	;match, strtrim(sp_name,2), strtrim(bz.gname,2), spind, bspind, count=spcount
	;if spcount gt 0 then begin
	;	bz[bspind].flux_radio = sp_radioflux[spind]
	;	bz[bspind].flux_xray = sp_xrayflux[spind]
	;	bz[bspind].spindex_ox = sp_aox[spind]
	;	bz[bspind].spindex_rx = sp_arx[spind]
	;	bz[bspind].spindex_ro = sp_aro[spind]
	;	bz[bspind].sp_mag_r = sp_mag_r[spind]
	;endif

	; For each galaxy, find the 500 kpc angular distance [arcsec]

	ang500 = zang(500., sp.z,/silent,/wmap7)
	sp.fieldsize = ang500
	Mr_star = -21.1837		; Blanton et al. (2003)
	sdss_lim = 22.2			; SDSS apparent limiting magnitude

	;ptrarr_bg_mag = ptrarr(nb,/allocate_heap)

	; Background galaxy array generation

	file_gname_north = [c3_file_gname,c10_file_gname]
	search_id_north = [c3_search_id,c10_search_id+5000]
	con_ra_north = [c3_ra, c10_ra]
	con_dec_north = [c3_dec, c10_dec]

	ncon_search_id_north = [cn3_search_id,cn10_search_id+5000]
	ncon_ra_north = [cn3_ra, cn10_ra]
	ncon_dec_north = [cn3_dec, cn10_dec]
	ncon_absmagR_north = [cn3_absmagR, cn10_absmagR]

	file_gname_south = [c3d_file_gname,c10d_file_gname]
	search_id_south = [c3d_search_id,c10d_search_id+5000]
	con_ra_south = [c3d_ra, c10d_ra]
	con_dec_south = [c3d_dec, c10d_dec]

	ncon_search_id_south = [cn3d_search_id,cn10d_search_id+5000]
	ncon_ra_south = [cn3d_ra, cn10d_ra]
	ncon_dec_south = [cn3d_dec, cn10d_dec]
	ncon_absmagR_south = [cn3d_absmagR, cn10d_absmagR]

	for i=0,nb-1 do begin

		; Find number of neighbors within projected 500 Mpc

		gcirc, 2, sp[i].ra, sp[i].dec, *(sp[i].n_ra), *(sp[i].n_dec), angdist
		ind3 = where(angdist lt ang500[i], ind3count)
		sp[i].n500 = ind3count

		; Find the number of neighbors within 500 Mpc and APPARENT magnitude brighter than counting mag

		distance = lumdist(sp[i].z,/silent,/wmap7)
		app_mstar = (Mr_star+2) + 5*alog10(distance*1e6) - 5
		countingmag_app = sdss_lim < app_mstar

		;countingmag_app = sdss_lim
		;countingmag_app = app_mstar

		ind3_cmagapp = where((angdist lt ang500[i]) and (*(sp[i].n_mag_r) lt countingmag_app), cmagcount_app)
		sp[i].cmag_app = countingmag_app

			; Set the counting magnitude to ABSOLUTE, rather than apparent magnitude

			sdss_lim_abs = sdss_lim - 5*alog10(distance*1e6) + 5
			abs_mstar = Mr_star + 2
			countingmag_abs = sdss_lim_abs < abs_mstar

			;countingmag_abs = sdss_lim_abs
			;countingmag_abs = abs_mstar

			abs_rmag = *(sp[i].n_absmag_r)
			ind3_cmagabs = where((angdist lt ang500[i]) and $
				(abs_rmag lt countingmag_abs), cmagcount_abs)

		n500_cmag = cmagcount_abs
		cmag = countingmag_abs

		sp[i].n500_cmag = n500_cmag
		sp[i].cmag = cmag

		; Northern control field

		ind_north = where(sp[i].gname eq file_gname_north,icount_north)
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

				nind = where((absmagR_temp lt sp[i].cmag) and (absmagR_temp gt -99) and (con_angdist lt ang500[i]), nc_north)
				;if nc_north gt 0 then *(ptrarr_bg_mag[i]) = absmagr_temp[nind]
				bgtemp_north = nc_north
			endif
		endif else begin
			nc_north = 0
			bgtemp_north = 0
		endelse

		; Southern control field

		ind_south = where(sp[i].gname eq file_gname_south,icount_south)
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

				nind = where((absmagR_temp lt sp[i].cmag) and (absmagR_temp gt -99) and (con_angdist lt ang500[i]), nc_south)
				;if nc_south gt 0 then *(ptrarr_bg_mag[i]) = absmagr_temp[nind]
				bgtemp_south = nc_south
			endif
		endif else begin
			nc_south = 0
			bgtemp_south = 0
		endelse

		if nc_north gt 0 then begin
			if nc_south gt 0 then sp[i].bg = mean([bgtemp_north,bgtemp_south]) else sp[i].bg = bgtemp_north
		endif else sp[i].bg = bgtemp_south

	endfor

	save, filename=savdir+'sp_modz_structure.sav', sp

	print,''
	print,'Spiral (modz) galaxy structure written to '+savdir+'sp_modz_structure.sav'
	print,''

	if keyword_set(stop) then stop
end
