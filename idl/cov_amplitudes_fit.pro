
;+
; NAME:
;       
;	COV_AMPLITUDES_FIT
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

phistararr = fillarr(0.01,0.,1.)
n = n_elements(phistararr)
slopearr = fltarr(n)
linfun = 'P[1] * X + P[0]'
for i=0,n-1 do begin

	cov_amplitudes,phistararr[i], barr_papers, barr_calc, /noplot
	funcreturn = mpfitexpr(linfun, barr_papers, barr_calc, sqrt(abs(barr_calc)), p, /quiet)
	slopearr[i] = funcreturn[1]

endfor

bestphistar = phistararr[closeto(slopearr,1)]

cov_amplitudes, bestphistar
cgplot, indgen(5000)-1000, indgen(5000)-1000, /overplot, linestyle=1

end
