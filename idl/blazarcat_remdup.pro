
;+
; NAME:
;       
;	BLAZARCAT_REMDUP
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

file = '~/Astronomy/Research/blazars/allblazars.cat'

readcol, file, bname, ra, dec, z, btype, $
	format='a,f,f,f,a', /silent

n = n_elements(bname)

count=0
i = 0

while i lt n do begin

	gcirc, 2, ra[i], dec[i], ra, dec, angdist
	deltaz = abs(z[i] - z)

	angind = where(angdist lt 2.,aa)
	if aa gt 1 then begin
		angind_remdup = setdifference(angind,[i])
		acount = n_elements(angind_remdup)

		if acount gt 0 then begin

			for j=0, acount-1 do begin

				if deltaz[angind_remdup[j]] lt 0.1 then begin
					;print, bname[i],         ra[i],         dec[i],         angdist[angind_remdup[j]], deltaz[angind_remdup[j]]
					;print, bname[angind_remdup[j]], ra[angind_remdup[j]], dec[angind_remdup[j]], angdist[angind_remdup[j]], deltaz[angind_remdup[j]]
					;print,''

					; Remove the second element from the list

					bname = bname[setdifference(indgen(n_elements(bname)),[angind_remdup[j]])]
					ra = ra[setdifference(indgen(n_elements(ra)),[angind_remdup[j]])]
					dec = dec[setdifference(indgen(n_elements(dec)),[angind_remdup[j]])]
					z = z[setdifference(indgen(n_elements(z)),[angind_remdup[j]])]
					btype = btype[setdifference(indgen(n_elements(btype)),[angind_remdup[j]])]

					n -= 1
					count+=1
				endif

			endfor

		endif
	endif

	i += 1

endwhile

print,''
print, 'Total number of duplicates found: ',count,'/',ntot
print,''

writefile = '~/Astronomy/Research/blazars/allblazars_remdup.cat'
openw, lun1, writefile, /get_lun
for i=0,n-1 do printf, lun, bname[i], string(ra[i]), string(dec[i]), string(z[i]), '   ', btype[i]
close, lun1 & free_lun, lun1

end
