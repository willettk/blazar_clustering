;+
; NAME:
;       
;	GZ_CONTROLFIELDS
;
; PURPOSE:
;
;	Generate control fields, split into 3 and 10 arcminute samples for Galaxy Zoo-selected objects
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
;	There are limits to the GZ catalog in redshift space; I cannot go out to the same distance for user-classified as I can for blazars. 
;
;	ZooNoSpec: 40 elliptical galaxies with photometric redshifts z > 0.5, with the maximum at z=0.605
;
;	ZooSpec: 28 elliptical galaxies with spectroscopic redshifts z > 0.35; maximum is 4.0, but likely unreliable. 
;	         370 elliptical galaxies with spectroscopic redshifts z > 0.30
;
; REVISION HISTORY
;       Written by K. Willett                Apr 12
;-

file_sp='~/Astronomy/Research/blazars/zoonospec_spiral_willettk.csv'

readcol, file_sp, $
	objID_zns_sp, ra_zns_sp, dec_zns_sp, p_sp_zns,z_zns_sp,zerr_zns_sp, $
	;objID,ra,dec,redshift,redshiftErr,p_el_debiased,p_cs_debiased,search_id, $
	/skip, $
	/silent, $
	format='a,f,f,f,f,f'

name_array = 'sp_zns_'+strtrim(lindgen(n_elements(objid_zns_sp)),2)

ra=ra_zns_sp
dec=dec_zns_sp
redshift = z_zns_sp
redshifterr = zerr_zns_sp

; Trim spiral galaxies based on the blazar redshift distribution

	blazardir = '~/Astronomy/Research/blazars/'
	restore,blazardir+'ab_structure.sav'

	lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	neighbors = 10
	ind = where(bz.z lt highz and bz.z gt lowz and bz.n500_cmag ge neighbors and bz.bg gt 0)
	bz=bz[ind]

result = histogram(bz.z,bins=delz,omin=omin,omax=omax)
bins = fillarr(delz,omin,omax)

; For each bin in the MODIFIED blazar redshift distribution, randomly select the same number of spiral galaxies to fill it.

breakz = 0.162d

binarr10 = fillarr(delz,lowz,breakz)
totlow = where(bz.z le breakz,nlow)
n10 = n_elements(binarr10)
sp_10arcmin_ind = lonarr(nlow)
j=0

for i=0, n10-2 do begin
	bjunk = where(bz.z gt binarr10[i] and bz.z le binarr10[i+1],nb10)
	kjunk = where(redshift gt binarr10[i] and redshift le binarr10[i+1],nk10)
	if (nb10 gt 0) and (nk10 gt 0) and (nk10 gt nb10) then sp_10arcmin_ind[j:j+nb10-1] = kjunk[cgRandomIndices(nk10,nb10,seed=seed10)] else print,binarr10[i],nb10,nk10
	j+=nb10
endfor
bjunk_10last = where(bz.z gt binarr10[n10-1] and bz.z le breakz,nb10_last)
kjunk_10last = where(redshift gt binarr10[n10-1] and redshift le breakz,nk10_last)
if (nb10_last gt 0) and (nk10_last gt 0) then sp_10arcmin_ind[j:j+nb10_last-1] = kjunk_10last[cgRandomIndices(nk10_last,nb10_last,seed=seed2)]

binarr3 = fillarr(delz,breakz,max(bz.z))
tothigh = where(bz.z gt breakz,nhigh)
n3 = n_elements(binarr3)
sp_3arcmin_ind = lonarr(nhigh)
j=0

for i=0, n3-2 do begin
	bjunk = where(bz.z gt binarr3[i] and bz.z le binarr3[i+1],nb3)
	kjunk = where(redshift gt binarr3[i] and redshift le binarr3[i+1],nk3)
	if (nb3 gt 0) and (nk3 gt 0) and (nk3 gt nb3) then sp_3arcmin_ind[j:j+nb3-1] = kjunk[cgRandomIndices(nk3,nb3,seed=seed3)] else print,binarr3[i],nb3,nk3
	j+=nb3
endfor
bjunk_3last = where(bz.z gt binarr3[n3-1] and bz.z le max(bz.z),nb3_last)
kjunk_3last = where(redshift gt binarr3[n3-1] and redshift lt max(bz.z),nk3_last)
if (nb3_last gt 0) and (nk3_last gt 0) then sp_3arcmin_ind[j:j+nb3_last-1] = kjunk[cgRandomIndices(nk3_last,nb3_last,seed=seed4)]

sp_3arcmin_ind = sp_3arcmin_ind[where(sp_3arcmin_ind ne 0)]	; Weird kluge; unhappy about this. 

; Write a CSV file on which SDSS can search for photometry and nearby neighbors

spiral_arr3 = [transpose(string(name_array[sp_3arcmin_ind],format='(a15)')), transpose(string(ra_zns_sp[sp_3arcmin_ind])), transpose(string(dec_zns_sp[sp_3arcmin_ind])), transpose(string(redshift[sp_3arcmin_ind])),transpose(string(p_sp_zns[sp_3arcmin_ind]))]
spiral_arr10 = [transpose(string(name_array[sp_10arcmin_ind],format='(a15)')), transpose(string(ra_zns_sp[sp_10arcmin_ind])), transpose(string(dec_zns_sp[sp_10arcmin_ind])), transpose(string(redshift[sp_10arcmin_ind])),transpose(string(p_sp_zns[sp_10arcmin_ind]))]

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/sp_zns_modz_3arcmin.cat', /get_lun
printf, lun1, 'name              ra              dec            z           p_sp'
printf, lun1, spiral_arr3
close, lun1
free_lun, lun1

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/sp_zns_modz_10arcmin.cat', /get_lun
printf, lun1, 'name              ra              dec            z           p_sp'
printf, lun1, spiral_arr10
close, lun1
free_lun, lun1

; Make control field samples for both the 3 and 10 arcminute spiral galaxies

ra3 = ra[sp_3arcmin_ind]
dec3 = dec[sp_3arcmin_ind]

ranew3_north = ra3
ranew3_south = ra3
decnew3_north = dec3
decnew3_south = dec3
highdec = where(dec3 ge 89., hdcount)
lowdec = where(dec3 lt 89., ldcount)
if hdcount gt 0 then begin
	ranew3_north[highdec] = raflip(ra3[highdec],/deg)
	decnew3_south[highdec] = dec3 - 1.
endif
if ldcount gt 0 then begin
	decnew3_north[lowdec] = dec3 + 1.
	decnew3_south[lowdec] = dec3 - 1.
endif

sparr3_connorth = [transpose(string(name_array[sp_3arcmin_ind],format='(a15)')), transpose(string(ranew3_north)), transpose(string(decnew3_north)), transpose(string(redshift[sp_3arcmin_ind])),transpose(string(p_sp_zns[sp_3arcmin_ind]))]
sparr3_consouth = [transpose(string(name_array[sp_3arcmin_ind],format='(a15)')), transpose(string(ranew3_south)), transpose(string(decnew3_south)), transpose(string(redshift[sp_3arcmin_ind])),transpose(string(p_sp_zns[sp_3arcmin_ind]))]

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/sp_zns_modz_3arcmin_connorth.cat', /get_lun
printf, lun1, 'sp_id              ra              dec            z         p_sp'
printf, lun1, sparr3_connorth
close, lun1
free_lun, lun1

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/sp_zns_modz_3arcmin_consouth.cat', /get_lun
printf, lun1, 'sp_id              ra              dec            z         p_sp'
printf, lun1, sparr3_consouth
close, lun1
free_lun, lun1

ra10 = ra_array[sp_10arcmin_ind]
dec10 = dec_array[sp_10arcmin_ind]

ranew10_north = ra10
ranew10_south = ra10
decnew10_north = dec10
decnew10_south = dec10
highdec = where(dec10 ge 89., hdcount)
lowdec = where(dec10 lt 89., ldcount)
if hdcount gt 0 then begin
	ranew10_north[highdec] = raflip(ra10[highdec],/deg)
	decnew10_south[highdec] = dec10 - 1.
endif
if ldcount gt 0 then begin
	decnew10_north[lowdec] = dec10 + 1.
	decnew10_south[lowdec] = dec10 - 1.
endif

sparr10_connorth = [transpose(string(name_array[sp_10arcmin_ind],format='(a15)')), transpose(string(ranew10_north)), transpose(string(decnew10_north)), transpose(string(redshift[sp_10arcmin_ind])),transpose(string(p_sp_zns[sp_10arcmin_ind]))]
sparr10_consouth = [transpose(string(name_array[sp_10arcmin_ind],format='(a15)')), transpose(string(ranew10_south)), transpose(string(decnew10_south)), transpose(string(redshift[sp_10arcmin_ind])),transpose(string(p_sp_zns[sp_10arcmin_ind]))]

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/sp_zns_modz_10arcmin_connorth.cat', /get_lun
printf, lun1, 'sp_id              ra              dec            z        p_sp'
printf, lun1, sparr10_connorth
close, lun1
free_lun, lun1

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/sp_zns_modz_10arcmin_consouth.cat', /get_lun
printf, lun1, 'sp_id              ra              dec            z        p_sp'
printf, lun1, sparr10_consouth
close, lun1
free_lun, lun1

; Plot positions in an Aitoff projection

!p.multi=[0,2,2]

cghistoplot, redshift, /trad, charsize=1.5, bin=0.02, xr=[0,1], /ylog, /min_val,title='Spiral galaxies',xtitle='Redshift',datacolor='blue'
cghistoplot, /oplot, bz.z, datacolor='green', bin=0.02

aitoff_grid, label=2, /new
aitoff, [ra3,ra10], [dec3,dec10], x, y
aitoff, ranew3_north, decnew3_north, x3n, y3n
aitoff, ranew3_south, decnew3_south, x3s, y3s
aitoff, ranew10_north, decnew10_north, x10n, y10n
aitoff, ranew10_south, decnew10_south, x10s, y10s
cgplots, x, y, psym=2, color="Yellow"
cgplots, x10n, y10n, psym=2, color="Red"
cgplots, x10s, y10s, psym=2, color="Green"
cgplots, x3n, y3n, psym=2, color="Red"
cgplots, x3s, y3s, psym=2, color="Green"


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Elliptical galaxies

file_el='~/Astronomy/Research/blazars/zoonospec_elliptical_willettk.csv'

readcol, file_el, $
	objID_zns_el, ra_zns_el, dec_zns_el, p_el_zns,z_zns_el,zerr_zns_el, $
	;objID,ra,dec,redshift,redshiftErr,p_el_debiased,p_cs_debiased,search_id, $
	/skip, $
	/silent, $
	format='a,f,f,f,f,f'

name_array = 'el_zns_'+strtrim(lindgen(n_elements(objid_zns_el)),2)

ra=ra_zns_el
dec=dec_zns_el
redshift = z_zns_el
redshifterr = zerr_zns_el

; Trim elliptical galaxies based on the blazar redshift distribution

	blazardir = '~/Astronomy/Research/blazars/'
	restore,blazardir+'ab_structure.sav'

	lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	neighbors = 10
	ind = where(bz.z lt highz and bz.z gt lowz and bz.n500_cmag ge neighbors and bz.bg gt 0)
	bz=bz[ind]

result = histogram(bz.z,bins=delz,omin=omin,omax=omax)
bins = fillarr(delz,omin,omax)

; For each bin in the MODIFIED blazar redshift distribution, randomly select the same number of elliptical galaxies to fill it.

breakz = 0.162d

binarr10 = fillarr(delz,lowz,breakz)
totlow = where(bz.z le breakz and bz.z gt lowz,nlow)
n10 = n_elements(binarr10)
el_10arcmin_ind = lonarr(nlow)
j=0

for i=0, n10-2L do begin
	bjunk = where(bz.z gt binarr10[i] and bz.z le binarr10[i+1],nb10)
	kjunk = where(redshift gt binarr10[i] and redshift le binarr10[i+1],nk10)
	if (nb10 gt 0) and (nk10 gt 0) and (nk10 gt nb10) then el_10arcmin_ind[j:j+nb10-1] = kjunk[cgRandomIndices(nk10,nb10,seed=seed10)] else print,binarr10[i],nb10,nk10
	j+=nb10
endfor
bjunk_10last = where(bz.z gt binarr10[n10-1] and bz.z le breakz,nb10_last)
kjunk_10last = where(redshift gt binarr10[n10-1] and redshift le breakz,nk10_last)
if (nb10_last gt 0) and (nk10_last gt 0) then el_10arcmin_ind[j:j+nb10_last-1] = kjunk_10last[cgRandomIndices(nk10_last,nb10_last,seed=seed2)]

el_10arcmin_ind = el_10arcmin_ind[where(el_10arcmin_ind ne 0)]	; Weird kluge; unhappy about this. 

binarr3 = fillarr(delz,breakz,max(bz.z))
tothigh = where(bz.z gt breakz,nhigh)
n3 = n_elements(binarr3)
el_3arcmin_ind = lonarr(nhigh)
j=0

for i=0, n3-2L do begin
	bjunk = where(bz.z gt binarr3[i] and bz.z le binarr3[i+1],nb3)
	kjunk = where(redshift gt binarr3[i] and redshift le binarr3[i+1],nk3)
	if (nb3 gt 0) and (nk3 gt 0) and (nk3 gt nb3) then el_3arcmin_ind[j:j+nb3-1] = kjunk[cgRandomIndices(nk3,nb3,seed=seed3)] else print,binarr3[i],nb3,nk3
	j+=nb3
endfor
bjunk_3last = where(bz.z gt binarr3[n3-1] and bz.z le max(bz.z),nb3_last)
kjunk_3last = where(redshift gt binarr3[n3-1] and redshift lt max(bz.z),nk3_last)
if (nb3_last gt 0) and (nk3_last gt 0) then el_3arcmin_ind[j:j+nb3_last-1] = kjunk[cgRandomIndices(nk3_last,nb3_last,seed=seed4)]

el_3arcmin_ind = el_3arcmin_ind[where(el_3arcmin_ind ne 0)]	; Weird kluge; unhappy about this. 

; Write a CSV file on which SDSS can search for photometry and nearby neighbors

elliptical_arr3 = [transpose(string(name_array[el_3arcmin_ind],format='(a15)')), transpose(string(ra_zns_el[el_3arcmin_ind])), transpose(string(dec_zns_el[el_3arcmin_ind])), transpose(string(redshift[el_3arcmin_ind])),transpose(string(p_el_zns[el_3arcmin_ind]))]
elliptical_arr10 = [transpose(string(name_array[el_10arcmin_ind],format='(a15)')), transpose(string(ra_zns_el[el_10arcmin_ind])), transpose(string(dec_zns_el[el_10arcmin_ind])), transpose(string(redshift[el_10arcmin_ind])),transpose(string(p_el_zns[el_10arcmin_ind]))]

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/el_zns_modz_3arcmin.cat', /get_lun
printf, lun1, 'name              ra              dec            z           p_el'
printf, lun1, elliptical_arr3
close, lun1
free_lun, lun1

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/el_zns_modz_10arcmin.cat', /get_lun
printf, lun1, 'name              ra              dec            z           p_el'
printf, lun1, elliptical_arr10
close, lun1
free_lun, lun1

; Make control field samples for both the 3 and 10 arcminute elliptical galaxies

ra3 = ra[el_3arcmin_ind]
dec3 = dec[el_3arcmin_ind]

ranew3_north = ra3
ranew3_south = ra3
decnew3_north = dec3
decnew3_south = dec3
highdec = where(dec3 ge 89., hdcount)
lowdec = where(dec3 lt 89., ldcount)
if hdcount gt 0 then begin
	ranew3_north[highdec] = raflip(ra3[highdec],/deg)
	decnew3_south[highdec] = dec3 - 1.
endif
if ldcount gt 0 then begin
	decnew3_north[lowdec] = dec3 + 1.
	decnew3_south[lowdec] = dec3 - 1.
endif

elarr3_connorth = [transpose(string(name_array[el_3arcmin_ind],format='(a15)')), transpose(string(ranew3_north)), transpose(string(decnew3_north)), transpose(string(redshift[el_3arcmin_ind])),transpose(string(p_el_zns[el_3arcmin_ind]))]
elarr3_consouth = [transpose(string(name_array[el_3arcmin_ind],format='(a15)')), transpose(string(ranew3_south)), transpose(string(decnew3_south)), transpose(string(redshift[el_3arcmin_ind])),transpose(string(p_el_zns[el_3arcmin_ind]))]

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/el_zns_modz_3arcmin_connorth.cat', /get_lun
printf, lun1, 'el_id              ra              dec            z        p_el'
printf, lun1, elarr3_connorth
close, lun1
free_lun, lun1

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/el_zns_modz_3arcmin_consouth.cat', /get_lun
printf, lun1, 'el_id              ra              dec            z        p_el'
printf, lun1, elarr3_consouth
close, lun1
free_lun, lun1

ra10 = ra_array[el_10arcmin_ind]
dec10 = dec_array[el_10arcmin_ind]

ranew10_north = ra10
ranew10_south = ra10
decnew10_north = dec10
decnew10_south = dec10
highdec = where(dec10 ge 89., hdcount)
lowdec = where(dec10 lt 89., ldcount)
if hdcount gt 0 then begin
	ranew10_north[highdec] = raflip(ra10[highdec],/deg)
	decnew10_south[highdec] = dec10 - 1.
endif
if ldcount gt 0 then begin
	decnew10_north[lowdec] = dec10 + 1.
	decnew10_south[lowdec] = dec10 - 1.
endif

elarr10_connorth = [transpose(string(name_array[el_10arcmin_ind],format='(a15)')), transpose(string(ranew10_north)), transpose(string(decnew10_north)), transpose(string(redshift[el_10arcmin_ind])),transpose(string(p_el_zns[el_10arcmin_ind]))]
elarr10_consouth = [transpose(string(name_array[el_10arcmin_ind],format='(a15)')), transpose(string(ranew10_south)), transpose(string(decnew10_south)), transpose(string(redshift[el_10arcmin_ind])),transpose(string(p_el_zns[el_10arcmin_ind]))]

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/el_zns_modz_10arcmin_connorth.cat', /get_lun
printf, lun1, 'el_id              ra              dec            z       p_el'
printf, lun1, elarr10_connorth
close, lun1
free_lun, lun1

openw, width=200, lun1, '~/Astronomy/Research/blazars/gz/el_zns_modz_10arcmin_consouth.cat', /get_lun
printf, lun1, 'el_id              ra              dec            z       p_el'
printf, lun1, elarr10_consouth
close, lun1
free_lun, lun1

; Plot positions in an Aitoff projection

cghistoplot, redshift, /trad, charsize=1.5, bin=0.02, xr=[0,1], /ylog, /min_val,title='Elliptical galaxies',xtitle='Redshift',datacolor='red'
cghistoplot, /oplot, bz.z, datacolor='green', bin=0.02

aitoff_grid, label=2, /new
aitoff, [ra3,ra10], [dec3,dec10], x, y
aitoff, ranew3_north, decnew3_north, x3n, y3n
aitoff, ranew3_south, decnew3_south, x3s, y3s
aitoff, ranew10_north, decnew10_north, x10n, y10n
aitoff, ranew10_south, decnew10_south, x10s, y10s
cgplots, x, y, psym=2, color="Yellow"
cgplots, x10n, y10n, psym=2, color="Red"
cgplots, x10s, y10s, psym=2, color="Green"
cgplots, x3n, y3n, psym=2, color="Red"
cgplots, x3s, y3s, psym=2, color="Green"

end
