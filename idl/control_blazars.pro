
;+
; NAME:
;       
;	CONTROL_BLAZARS
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
;	Added 3d and 10d samples, doubling the number of control fields 
;		by adding one at 1 degree south of each blazar. 		Mar 12
;-

readfile3='~/Astronomy/Research/blazars/allblazars_3arcmin.cat'
readfile10='~/Astronomy/Research/blazars/allblazars_10arcmin.cat'
writefile3='~/Astronomy/Research/blazars/allblazars_con_3arcmin.cat'
writefile10='~/Astronomy/Research/blazars/allblazars_con_10arcmin.cat'
writefile3d='~/Astronomy/Research/blazars/allblazars_con_d_3arcmin.cat'
writefile10d='~/Astronomy/Research/blazars/allblazars_con_d_10arcmin.cat'

; 3 arcmin sample

readcol, readfile3, $
	name3,ra3,dec3,z3,btype3, $
	format='a,f,f,f,a', $
	skipline=1, $
	/silent

decnew3 = dec3
decnew3d = dec3
highdec = where(dec3 ge 89., hdcount)
lowdec = where(dec3 lt 89., ldcount)
if hdcount gt 0 then decnew3[highdec] = dec3 - 1.
if ldcount gt 0 then begin
	decnew3[lowdec] = dec3 + 1.
	decnew3d[lowdec] = dec3 - 1.
endif

writearray3 = [transpose(name3), transpose(string(ra3)), transpose(string(decnew3)), transpose(string(z3)), transpose(string(btype3))]
writearray3d = [transpose(name3), transpose(string(ra3)), transpose(string(decnew3d)), transpose(string(z3)), transpose(string(btype3))]

openw, lun, writefile3, /get_lun
printf, lun, 'name		ra		dec		z	btype'
printf, lun, writearray3
close, lun
free_lun, lun

openw, lun, writefile3d, /get_lun
printf, lun, 'name		ra		dec		z	btype'
printf, lun, writearray3d
close, lun
free_lun, lun

; 10 arcmin sample

readcol, readfile10, $
	name10,ra10,dec10,z10,btype10, $
	format='a,f,f,f,a', $
	skipline=1, $
	/silent

decnew10 = dec10
decnew10d = dec10
highdec = where(dec10 ge 89., hdcount)
lowdec = where(dec10 lt 89., ldcount)
if hdcount gt 0 then decnew10[highdec] = dec10 - 1.
if ldcount gt 0 then begin
	decnew10[lowdec] = dec10 + 1.
	decnew10d[lowdec] = dec10 - 1.
endif


writearray10 = [transpose(name10), transpose(string(ra10)), transpose(string(decnew10)), transpose(string(z10)), transpose(string(btype10))]
writearray10d = [transpose(name10), transpose(string(ra10)), transpose(string(decnew10d)), transpose(string(z10)), transpose(string(btype10))]

openw, lun, writefile10, /get_lun
printf, lun, 'name		ra		dec		z	btype'
printf, lun, writearray10
close, lun
free_lun, lun

openw, lun, writefile10d, /get_lun
printf, lun, 'name		ra		dec		z	btype'
printf, lun, writearray10d
close, lun
free_lun, lun

stop

end
