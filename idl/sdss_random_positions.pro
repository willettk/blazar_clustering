
;+
; NAME:
;       
;	SDSS_RANDOM_POSITIONS
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
;       Written by K. Willett                
;-

;pro sdss_random_positions, n

file='~/Astronomy/Research/blazars/allblazars_remdup.cat'
n = numlines(file) - 1

; http://mathworld.wolfram.com/SpherePointPicking.html
; Create uniformly distributed arrays on the northern hemisphere in decimal degrees

ra_array = 360. * randomu(seed_ra, n)
dec_array = asin(randomu(seed_dec, n)) * 180./!dpi
name_array = 'RND'+strtrim(indgen(n)+1,2)

; Assign a redshift to each RANDOM point selected from the blazar distribution

readcol, file, $
	bname, bra, bdec, bz, btype, $
	format='a,f,f,f,a' , $
	/silent, $
	/skip

rand_3arcmin_ind = where(bz gt 0.043 and bz lt 0.162)
rand_10arcmin_ind = where(bz ge 0.162)

; Make control field samples for both the 3 and 10 arcminute random positions

ra3 = ra_array[rand_3arcmin_ind]
dec3 = dec_array[rand_3arcmin_ind]

ranew3_north = ra3
ranew3_south = ra3
decnew3_north = dec3
decnew3_south = dec3
highdec = where(dec3 ge 89., hdcount)
lowdec = where(dec3 lt 89., ldcount)
if hdcount gt 0 then begin
	ranew3_north[highdec] = raflip(ra3[highdec],/deg)
	decnew3_south[highdec] = dec3[highdec] - 1.
endif
if ldcount gt 0 then begin
	decnew3_north[lowdec] = dec3[lowdec] + 1.
	decnew3_south[lowdec] = dec3[lowdec] - 1.
endif

randarr3_connorth = [transpose(string(name_array[rand_3arcmin_ind],format='(a14)')), transpose(string(ranew3_north)), transpose(string(decnew3_north)), transpose(string(spec_redshift[rand_3arcmin_ind]))]
randarr3_consouth = [transpose(string(name_array[rand_3arcmin_ind],format='(a14)')), transpose(string(ranew3_north)), transpose(string(decnew3_south)), transpose(string(spec_redshift[rand_3arcmin_ind]))]

openw, lun1, '~/Astronomy/Research/blazars/rand_3arcmin_connorth.cat', /get_lun
printf, lun1, 'rand_id              ra              dec            z'
printf, lun1, randarr3_connorth
close, lun1
free_lun, lun1

openw, lun1, '~/Astronomy/Research/blazars/rand_3arcmin_consouth.cat', /get_lun
printf, lun1, 'rand_id              ra              dec            z'
printf, lun1, randarr3_consouth
close, lun1
free_lun, lun1

ra10 = ra_array[rand_10arcmin_ind]
dec10 = dec_array[rand_10arcmin_ind]

ranew10_north = ra10
ranew10_south = ra10
decnew10_north = dec10
decnew10_south = dec10
highdec = where(dec10 ge 89., hdcount)
lowdec = where(dec10 lt 89., ldcount)
if hdcount gt 0 then begin
	ranew10_north[highdec] = raflip(ra10[highdec],/deg)
	decnew10_south[highdec] = dec10[highdec] - 1.
endif
if ldcount gt 0 then begin
	decnew10_north[lowdec] = dec10[lowdec] + 1.
	decnew10_south[lowdec] = dec10[lowdec] - 1.
endif

randarr10_connorth = [transpose(string(name_array[rand_10arcmin_ind],format='(a14)')), transpose(string(ranew10_north)), transpose(string(decnew10_north)), transpose(string(spec_redshift[rand_10arcmin_ind]))]
randarr10_consouth = [transpose(string(name_array[rand_10arcmin_ind],format='(a14)')), transpose(string(ranew10_south)), transpose(string(decnew10_south)), transpose(string(spec_redshift[rand_10arcmin_ind]))]

openw, lun1, '~/Astronomy/Research/blazars/rand_10arcmin_connorth.cat', /get_lun
printf, lun1, 'rand_id              ra              dec            z'
printf, lun1, randarr10_connorth
close, lun1
free_lun, lun1

openw, lun1, '~/Astronomy/Research/blazars/rand_10arcmin_consouth.cat', /get_lun
printf, lun1, 'rand_id              ra              dec            z'
printf, lun1, randarr10_consouth
close, lun1
free_lun, lun1

; Plot random positions in an Aitoff projection

!p.multi=[0,1,1]
aitoff_grid, label=2, /new			; Can't change the center coordinates of map 
       					; maybe switch to matplotlib/kapteyn packages?
aitoff, ra_array, dec_array, x, y
aitoff, ranew10_north, decnew10_north, x10n, y10n
aitoff, ranew10_south, decnew10_south, x10s, y10s
cgplots, x, y, psym=3, color="Yellow"
cgplots, x10n, y10n, psym=3, color="Red"
cgplots, x10s, y10s, psym=3, color="Green"

end
