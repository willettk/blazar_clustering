
;+
; NAME:
;       
;	Q0332
;
; PURPOSE:
;
;	
;
; INPUTS:
;
;       phi - integrated luminosity function
;
;       z - redshift of central quasar/galaxy
;
;       D - angular size distance [Mpc]
;
;       scale - converts 170 kpc to arcsec at given redshift
;
;       mi - Faber magnitude?
;
;       mi0 - zero=point magnitude?
;
;       ar0 - magnitude of some kind for main field?
;
;       ar1 - magnitude of some kind for dedicated offset field?
;
;       rg0 - index of radio galaxy within data file?
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
;       Written by Cristina Ramos Almeida, IAC
;       See Ramos Almeida et al. (2013) for details on implementation
;-

PRO q0332, phi, z, D, scale, mi, mi0, ar0, ar1, rg0

; Default values

; phi = 0.005
; z = 0.310
; D = 948.1 (from CosmoCalc given the redshift)
; scale = 37.292351585219656 (from zang.pro given redshift)
; mi = doesn't matter (not used)
; mi0 = 
; ar0 = 
; ar1 = 
; rg0 = 


readcol, 'q0332.ASC',    id, x, y, aper, best,ebest, bk, $
isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent



r170kpc=170/scale			; in arcsec

rg=rg0-1
x_gal=x[rg]			; position of the radio galaxy in the image. 
y_gal=y[rg]

r=sqrt((x-x_gal)^2+(y-y_gal)^2)
r=r*0.146				;in arcsec


mr=mi0+ar0

stop

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r170kpc)


nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

print, ''
print, 'Number of galaxies in the target field within a 170 kpc  radius =',nele
print, ''

print, ''
print, "mr' Faber=",mi
print, "[m_auto,flag,star]=",best[rg],flag[rg],star[rg]
print, 'Decide whether or not to subtract the radio galaxy from nele_170'
print, ''

stop

nele_170=nele




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   DEDICATED OFFSET FIELD   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, 'q0332_offset.ASC',    id, x, y, aper, best,ebest, bk, $
isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent

x_cen=1555.069   	; center of the image - (CRPIX1,CRPIX2). 
y_cen=1152.044 

r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

r_image=168.192


mr=mi0+ar1


ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0

for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field within a 170 kpc radius =',nele
print, ''

nele_170_offset=nele


print, 'SUMMARY:'
print, '--------------------------------------'
print, 'Galaxies within 170 kpc target field =', nele_170
print, 'Galaxies within 170 kpc offset field =', nele_170_offset


stop


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   AVERAGED OFFSET FIELD    ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0023_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


x_cen=1555.069   	; center of the image - (CRPIX1,CRPIX2). 
y_cen=1152.044 

r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

r_image=168.192

ar=0.0416
mr=mi0+ar


ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 1 within a 170 kpc radius =',nele
print, ''

nele1=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0034_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent

r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.0647
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0

for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 2 within a 170 kpc radius =',nele
print, ''

nele2=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0035_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent

r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.08068
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 3 within a 170 kpc radius =',nele
print, ''

nele3=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0038_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.2693
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 4 within a 170 kpc radius =',nele
print, ''

nele4=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0039_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.0256
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 5 within a 170 kpc radius =',nele
print, ''

nele5=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0043_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.03538
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 6 within a 170 kpc radius =',nele
print, ''

nele6=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0213_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.05718
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 7 within a 170 kpc radius =',nele
print, ''

nele7=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0347_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.728
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 8 within a 170 kpc radius =',nele
print, ''

nele8=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0349_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.0253
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 9 within a 170 kpc radius =',nele
print, ''

nele9=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0404_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.8408
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 10 within a 170 kpc radius =',nele
print, ''

nele10=nele



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0442_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.0851
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 11 within a 170 kpc radius =',nele
print, ''

nele11=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0620_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.217179
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 12 within a 170 kpc radius =',nele
print, ''

nele12=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0625_35_offset.ASC', id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.172479
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 13 within a 170 kpc radius =',nele
print, ''

nele13=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0625_53_offset.ASC', id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.266
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 14 within a 170 kpc radius =',nele
print, ''

nele14=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0806_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.1788
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 15 within a 170 kpc radius =',nele
print, ''

nele15=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0859_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.3656
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 16 within a 170 kpc radius =',nele
print, ''

nele16=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0915_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.1323
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 17 within a 170 kpc radius =',nele
print, ''

nele17=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks0945_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.057
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 18 within a 170 kpc radius =',nele
print, ''

nele18=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks1151_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent

r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.2068
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 19 within a 170 kpc radius =',nele
print, ''

nele19=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks1355_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.227
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 20 within a 170 kpc radius =',nele
print, ''

nele20=nele



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks1559_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.240
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 21 within a 170 kpc radius =',nele
print, ''

nele21=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks1648_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.2506
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 22 within a 170 kpc radius =',nele
print, ''

nele22=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks1733_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.24017
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 23 within a 170 kpc radius =',nele
print, ''

nele23=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks1814_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.2255
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 24 within a 170 kpc radius =',nele
print, ''

nele24=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks1839_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.155
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 25 within a 170 kpc radius =',nele
print, ''

nele25=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks1932_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.1542
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 26 within a 170 kpc radius =',nele
print, ''

nele26=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks1934_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.238
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 27 within a 170 kpc radius =',nele
print, ''

nele27=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks1949_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.518
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 28 within a 170 kpc radius =',nele
print, ''

nele28=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/pks1954_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.1638
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 29 within a 170 kpc radius =',nele
print, ''

nele29=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks2135_14_offset.ASC', id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.147
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 30 within a 170 kpc radius =',nele
print, ''

nele30=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks2211_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.0661
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 31 within a 170 kpc radius =',nele
print, ''

nele31=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks2221_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.2453
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 32 within a 170 kpc radius =',nele
print, ''

nele32=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks2314_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.1667
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 33 within a 170 kpc radius =',nele
print, ''

nele33=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
readcol, '../catalogs/offset/pks2356_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.0419658
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 34 within a 170 kpc radius =',nele
print, ''

nele34=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



readcol, '../catalogs/offset/q0114_offset.ASC', id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong,ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.077265
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 35 within a 170 kpc radius =',nele
print, ''

nele35=nele

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/q0123_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.07658
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 36 within a 170 kpc radius =',nele
print, ''

nele36=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/q0217_00_offset.ASC',   id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.094786
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 37 within a 170 kpc radius =',nele
print, ''

nele37=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/q0217_01_offset.ASC',   id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.07504
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 38 within a 170 kpc radius =',nele
print, ''

nele38=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/q0218_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.07709
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 39 within a 170 kpc radius =',nele
print, ''

nele39=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/q0227_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.071538
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 40 within a 170 kpc radius =',nele
print, ''

nele40=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/q0234_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.0851
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 41 within a 170 kpc radius =',nele
print, ''

nele41=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/q0249_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.1452
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 42 within a 170 kpc radius =',nele
print, ''

nele42=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



readcol, '../catalogs/offset/q0320_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.296496
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 43 within a 170 kpc radius =',nele
print, ''

nele43=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



readcol, '../catalogs/offset/q0332_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.31897
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 44 within a 170 kpc radius =',nele
print, ''

nele44=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



readcol, '../catalogs/offset/q0334_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.2582
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 45 within a 170 kpc radius =',nele
print, ''

nele45=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



readcol, '../catalogs/offset/q0848_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.0916
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 46 within a 170 kpc radius =',nele
print, ''

nele46=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



readcol, '../catalogs/offset/q0904_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.07786
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 47 within a 170 kpc radius =',nele
print, ''

nele47=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/q0923_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.1068
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 48 within a 170 kpc radius =',nele
print, ''

nele48=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/q0924_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.1246
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 49 within a 170 kpc radius =',nele
print, ''

nele49=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/q0948_offset.ASC', id, x, y, aper, best,ebest, bk, isoarea, ra, dec, $
theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.2256
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 50 within a 170 kpc radius =',nele
print, ''

nele50=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol, '../catalogs/offset/q2358_offset.ASC',      id, x, y, aper, best,ebest, bk, isoarea, ra, dec,$
 theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.11128
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind) 
if ind[0] lt 0 then nele=0


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 51 within a 170 kpc radius =',nele
print, ''

nele51=nele



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

readcol, '../catalogs/offset/pks0105_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent

r=sqrt((x-x_cen)^2+(y-y_cen)^2)
r=r*0.146

ar=0.060
mr=mi0+ar

ind=where(star le 0.85 and flag lt 4 and best le (mr+2) and best ge (mr-1) and r le r_image)

nele=n_elements(ind)


for i=0,nele-1 do print, id(ind(i)),ra(ind(i)),dec(ind(i)),best(ind(i)),flag(ind(i)),fwhm(ind(i)),elong(ind(i)),star(ind(i))

nele=nele/(!pi*r_image*r_image)
nele=nele*(!pi*r170kpc*r170kpc)


print, ''
print, 'Number of galaxies in the offset field 52 within a 170 kpc radius =',nele
print, ''

nele52=nele


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;readcol, '../catalogs/offset/pks2250_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks1938_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks1602_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks1306_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks1547_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks1136_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks0117_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks0252_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks0235_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks2135_20_offset.ASC', id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent
;readcol, '../catalogs/offset/pks0409_offset.ASC',    id, x, y, best,ebest, bk, isoarea, ra, dec, theta, flag, fwhm, elong, ellip, star, format="i,f,f,f,f,f,f,f,f,f,i,f,f,f,f", /silent


a=[nele1,nele2,nele3,nele4,nele5,nele6,nele7,nele8,nele9,nele10,nele11,nele12,nele13,nele14,nele15,nele16,nele17,nele18,nele19,nele20,$
nele21,nele22,nele23,nele24,nele25,nele26,nele27,nele28,nele29,nele30,nele31,nele32,nele33,nele34,nele35,nele36,nele37,nele38,nele39,nele40,$
nele41,nele42,nele43,nele44,nele45,nele46,nele47,nele48,nele49,nele50,nele51,nele52]




limit=mean(a)+sqrt(mean(a))
ind=where(a lt limit)
a=a(ind)

hlp, a

stop


print, 'SUMMARY:'
print, '--------------------------------------'
print, 'Galaxies target field =', nele_170
print, 'Galaxies dedicated offset field =', nele_170_offset
print, 'Galaxies averaged offset field =', mean(a)
print, 'Galaxies median offset field =', median(a),'+/-',stddev(a)
print, '--------------------------------------'



theta= r170kpc*!pi/648000.
Nt=nele_170


Nb=nele_170_offset
Agq=(Nt/Nb-1)*((3-1.77)/2)*theta^(0.77) 
Ng=Nb/(!pi*theta*theta) 
Bgq=(Agq*Ng/(3.78*phi))*(D/(1+z))^(1.77-3)
print, 'dedicated offset -> Bgq =', Bgq

Nb=mean(a)
Agq=(Nt/Nb-1)*((3-1.77)/2)*theta^(0.77) 
Ng=Nb/(!pi*theta*theta) 
Bgq1=(Agq*Ng/(3.78*phi))*(D/(1+z))^(1.77-3)
print, 'averaged offset -> Bgq =', Bgq1

Nb=median(a)
Agq=(Nt/Nb-1)*((3-1.77)/2)*theta^(0.77) 
Ng=Nb/(!pi*theta*theta) 
Bgq2=(Agq*Ng/(3.78*phi))*(D/(1+z))^(1.77-3)
print, 'median offset -> Bgq =', Bgq2

fname='../RESULTS/qsos_rev.dat'
openw,1,fname, /append

printf, 1 , String(z, FORMAT='(F0.3)'), String(nele_170, nele_170_offset,Bgq, mean(a), Bgq1,median(a),stddev(a), Bgq2,FORMAT='(F0.2)')

close, 1

END





