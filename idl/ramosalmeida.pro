;+
; NAME:
;       
;	RAMOSALMEIDA
;
; PURPOSE:
;
;	Compare Bgb values given in Ramos Almeida+2013 to those I calculate based on same Nt, Nb. Way of
;   checking the code of BGB.pro
;
; INPUTS:
;
;
; OUTPUTS:
;
;
; KEYWORDS:
;
;
; EXAMPLE:
;
;
; NOTES:
;
;	
; REVISION HISTORY
;       Written by K. Willett                Jul 14
;-


readcol, '/Users/willettk/Astronomy/Research/blazars/ramosalmeida2013.csv', PKSID, z, Optical, Radio, Nt, Nb, Bgq, Navb, sigma, Bavgq,deltaBavgq, Nmedb, Bmedgq, format='a,f,a,a,f,f,f,f,f,f,f,f',/skipline,/silent

Mr_star_arr = [-21.43,-22.08,-23.27,-23.59,-23.84]
zarr = (indgen(5)+1)/5.
cmag_arr = fltarr(n_elements(z))
for i = 0,n_elements(z) - 1 do begin
    zind = where(z[i] lt zarr,zc)
    if zc gt 0 then begin
        cmag_arr[i] = Mr_star_arr[zind[0]]+2
    endif
endfor

bgb,nt,nb,zang(170.,z),z,cmag_arr,bgb_ra,bgb_ra_err,lf='ramosalmeida'

cgplot,bgq,bgq-bgb_ra,psym=9,xtitle='B [CRA]',ytitle='B [CRA] - B [KWW]',charsize=1.5 & cgplot,/over,[-10000,10000],[0,0],linestyle=2

end


