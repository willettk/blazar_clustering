
;+
; NAME:
;       
;	YG87_FIG1
;
; PURPOSE:
;
;	Reproduce Yee & Green (1987), Figure 1
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

!p.multi=[0,1,1]

; Read in digitized data points from Fig. 1 (WebPlotDigitizer)

datafile='~/Astronomy/Research/blazars/yg87_fig1.csv'
readcol, datafile, $
	/quiet, $
	mr, npermag, $
	format='f,f'

errfile='~/Astronomy/Research/blazars/yg87_fig1_err.csv'
readcol, errfile, $
	/quiet, $
	mr_err, npermag_col, $
	format='f,f'

npermag_err = fltarr(5)
for i=0,4 do npermag_err[i] = mean([npermag_col[i*2]-npermag[i],npermag[i]-npermag_col[i*2+1]])

ploterror, mr, npermag, npermag_err, $
	color="blue", $
	psym=16, $
	xtitle='M!IR!N (obs)', $
	ytitle='Relative number/mag', $
	xrange=[-19,-25], /xstyle, $
	yrange=[1d-1,2d2], /ystyle, $
	type=1
	
; Fit a Schechter function

parinfo = replicate({value:0.,fixed:0},3)
parinfo[*].value = [1.,-21.73,-1.0]
parinfo[2].fixed = 1
result = mpfitfun('schechter_diff', mr, npermag, npermag_err, parinfo=parinfo,/quiet)

print,result

xarr = fillarr(0.1,-24,-20)
cgplot, xarr, schechter_diff(xarr, result), $
	color='black', $
	/over, $
	thick=2

end
