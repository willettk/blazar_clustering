
;+
; NAME:
;       
;	TEVCAT_CONVERT
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
;       Written by K. Willett                 Aug 11
;-

file = '~/Astronomy/Research/blazars/tevcat_blazars.txt'

readcol, file, name, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, z, btype, $
	format='a,i,i,f,i,i,f,f,a'

n = n_elements(name)

ra_deg = fltarr(n)
dec_deg = fltarr(n)

for i=0, n-1 do begin
	ra_deg[i] = ten(ra_h[i], ra_m[i], ra_s[i])*15.
	dec_deg[i] = ten(dec_d[i], dec_m[i], dec_s[i])
	print, ra_deg[i], dec_deg[i]
endfor

writefile='~/Astronomy/Research/blazars/tevcat_blazars_dec.txt'
openw, lun, writefile, /get_lun
for i=0,n-1 do printf, lun, name[i], string(ra_deg[i]), string(dec_deg[i]), string(z[i]), '   ', btype[i]
close, lun & free_lun, lun

end
