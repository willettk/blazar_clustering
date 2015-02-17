
;+
; NAME:
;       
;	CONFIG_READ
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
;       Written by K. Willett                Apr 12
;-

file='~/Astronomy/Research/blazars/config.fit'

f = mrdfits(file,1,/silent)

r10 = where(f.z gt 0.043 and f.z lt 0.162)
r3 = where(f.z ge 0.162)

name='config_s'+strtrim(f.s,2)+'c'+strtrim(f.config,2)

array3 = [string(transpose(name[r3])),string(transpose(f[r3]._raj2000)),string(transpose(f[r3]._dej2000)),string(transpose(f[r3].mt)),string(transpose(f[r3].z)),string(transpose(f[r3].e_z))]
array10 = [string(transpose(name[r10])),string(transpose(f[r10]._raj2000)),string(transpose(f[r10]._dej2000)),string(transpose(f[r10].mt)),string(transpose(f[r10].z)),string(transpose(f[r10].e_z))]

config3 = '~/Astronomy/Research/blazars/config3.cat'
openw, lun1, config3, /get_lun
printf,lun1,'cname             ra            dec           mt          z       z_err'
printf,lun1,array3
close, lun1
free_lun, lun1

config10 = '~/Astronomy/Research/blazars/config10.cat'
openw, lun1, config10, /get_lun
printf,lun1,'cname             ra            dec           mt          z       z_err'
printf,lun1,array10
close, lun1
free_lun, lun1

; Make control field samples for both the 3 and 10 arcminute radio galaxies

ra3 = f[r3]._raj2000
dec3 = f[r3]._dej2000
ra10 = f[r10]._raj2000
dec10 = f[r10]._dej2000

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

config3_connorth = [string(transpose(name[r3])),string(transpose(ranew3_north)),string(transpose(decnew3_north)),string(transpose(f[r3].mt)),string(transpose(f[r3].z)),string(transpose(f[r3].e_z))]
config3_consouth = [string(transpose(name[r3])),string(transpose(ranew3_south)),string(transpose(decnew3_south)),string(transpose(f[r3].mt)),string(transpose(f[r3].z)),string(transpose(f[r3].e_z))]

openw, lun1, '~/Astronomy/Research/blazars/config3_connorth.cat', /get_lun
printf,lun1,'cname             ra            dec           mt          z       z_err'
printf, lun1, config3_connorth
close, lun1
free_lun, lun1

openw, lun1, '~/Astronomy/Research/blazars/config3_consouth.cat', /get_lun
printf,lun1,'cname             ra            dec           mt          z       z_err'
printf, lun1, config3_consouth
close, lun1
free_lun, lun1

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

config10_connorth = [string(transpose(name[r10])),string(transpose(ranew10_north)),string(transpose(decnew10_north)),string(transpose(f[r10].mt)),string(transpose(f[r10].z)),string(transpose(f[r10].e_z))]
config10_consouth = [string(transpose(name[r10])),string(transpose(ranew10_south)),string(transpose(decnew10_south)),string(transpose(f[r10].mt)),string(transpose(f[r10].z)),string(transpose(f[r10].e_z))]

openw, lun1, '~/Astronomy/Research/blazars/config10_connorth.cat', /get_lun
printf,lun1,'cname             ra            dec           mt          z       z_err'
printf, lun1, config10_connorth
close, lun1
free_lun, lun1

openw, lun1, '~/Astronomy/Research/blazars/config10_consouth.cat', /get_lun
printf,lun1,'cname             ra            dec           mt          z       z_err'
printf, lun1, config10_consouth
close, lun1
free_lun, lun1

end
