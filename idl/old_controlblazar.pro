pro old_controlblazar

; Old code for directly matching control fields and computing background counts for blazars in BGB_BLAZARS.pro

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;;;;;;;;;;;; Try the control field method of computing number of background galaxies ;;;;;;;;;;;;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	if keyword_set(readfiles) then begin

		file3arcmin_neighbors = blazardir+'ab_con_d_neighbors3_willettk.csv'
		file3arcmin_mags = blazardir+'ab_con_d_3arcmin_mags_willettk.csv'
		file10arcmin_neighbors = blazardir+'ab_con_d_neighbors10_willettk.csv'
		file10arcmin_mags = blazardir+'ab_con_d_10arcmin_mags_willettk.csv'

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
			cn3_mag_search_id,cn3_mag_objid, $
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
			cn10_mag_search_id,cn10_mag_objid,$
			skip=1, $
			delimiter=',', $
			format='f,f,f,f,f,f,f,f,f,f,f,f,a,f,f,f,f,f,f,f,f,f,f,f,f,a,a'

		save, filename=blazardir+'blazar_controlfield.sav', $
			c3_file_bname,c3_ra,c3_dec,c3_z,c3_btype,c3_search_id,c3_matched_id, $
			cn3_u,cn3_g,cn3_r,cn3_i,cn3_z,cn3_err_u,cn3_err_g,cn3_err_r,cn3_err_i,cn3_err_z,$
			cn3_ra,cn3_dec,cn3_type,cn3_redshift,cn3_redshift_err,$
			cn3_absmagU,cn3_absmagG,cn3_absmagR,cn3_absmagI,cn3_absmagZ,$
			cn3_kcorrU,cn3_kcorrG,cn3_kcorrR,cn3_kcorrI,cn3_kcorrZ,$
			cn3_mag_search_id,cn3_mag_objid, $
			c10_file_bname,c10_ra,c10_dec,c10_z,c10_btype,c10_search_id,c10_matched_id, $
			cn10_u,cn10_g,cn10_r,cn10_i,cn10_z,cn10_err_u,cn10_err_g,cn10_err_r,cn10_err_i,cn10_err_z,$
			cn10_ra,cn10_dec,cn10_type,cn10_redshift,cn10_redshift_err,$
			cn10_absmagU,cn10_absmagG,cn10_absmagR,cn10_absmagI,cn10_absmagZ,$
			cn10_kcorrU,cn10_kcorrG,cn10_kcorrR,cn10_kcorrI,cn10_kcorrZ,$
			cn10_mag_search_id,cn10_mag_objid
	
	endif else restore, blazardir+'blazar_controlfield.sav'

	; ra, dec, file_bname, search_id, absmagR

	ra = [c3_ra, c10_ra]
	dec = [c3_dec, c10_dec]
	file_bname = [c3_file_bname,c10_file_bname]
	search_id = [c3_search_id,c10_search_id+5000]

	mag_search_id = [cn3_mag_search_id,cn10_mag_search_id+5000]
	mag_ra = [cn3_ra, cn10_ra]
	mag_dec = [cn3_dec, cn10_dec]
	absmagR = [cn3_absmagR, cn10_absmagR]

	match, bz.bname, file_bname[rem_dup(file_bname)], ia, ib, count=matchcount
	unique_bname = (file_bname[rem_dup(file_bname)])[ib]
	unique_search_id = (search_id[rem_dup(file_bname)])[ib]
	unique_ra = (ra[rem_dup(file_bname)])[ib]
	unique_dec = (dec[rem_dup(file_bname)])[ib]

	ncarr = intarr(matchcount)
	for i=0,matchcount-1 do begin
		matchind = where(mag_search_id eq unique_search_id[i],n)
		if n gt 0 then begin
			ang500 = zang(500., bz[ia[i]].z,/silent,/wmap7)
			gcirc, 2, unique_ra[i], unique_dec[i], mag_ra[matchind], mag_dec[matchind], angdist

			nind = where((absmagR[matchind] lt bz[ia[i]].cmag) and (absmagR[matchind] gt -99) and (angdist lt ang500), nc)
			ncarr[i] = nc

		endif
	endfor

	bz = bz[ia]
	nb = ncarr
		
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;


			; Different method (which didn't work - big underprediction) using method of Smith, O'Dea, & Baum (1995)
			
			bg_inds = where((angdist gt ang1000[i]) and $
				(*(bz[i].n_absmag_r) lt cmag), bgcount)
			bg_mag = (*(bz[i].n_absmag_r))[bg_inds]
			if bz[i].z gt 0.162 then outer_radius = 3. else outer_radius = 10.
			area = !dpi * (outer_radius^2 - (ang1000[i]/60.)^2)	; arcmin
			bz[i].bg = bgcount
			bz[i].area = area
			*(ptrarr_bg_mag[i]) = bg_mag

end
