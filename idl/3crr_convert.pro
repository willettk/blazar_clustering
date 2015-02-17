
;+
; NAME:
;       
;	3CRR_CONVERT
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

file='~/Astronomy/Research/blazars/3crr.txt'

readcol, file, $
	name_3crr, name_iau, z, fr, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, $
	format='a,a,f,i,i,i,f,i,i,f', /silent

n = n_elements(z)

ra_deg = fltarr(n)
dec_deg = fltarr(n)

for i=0,n-1 do begin
	ra_deg[i]=ten(ra_h[i],ra_m[i],ra_s[i])*15.
	dec_deg[i]=ten(dec_d[i],dec_m[i],dec_s[i])
endfor

writearray = [transpose(name_3crr), transpose(name_iau), transpose(string(z)), transpose(string(fr)), transpose(string(ra_deg)), transpose(string(dec_deg))]

writefile='~/Astronomy/Research/blazars/3crr_dec.txt'
openw, lun, writefile, /get_lun
printf, lun, writearray
close, lun
free_lun, lun

end
