
;+
; NAME:
;       
;	MEYER_BE
;
; PURPOSE:
;
;	Plot the blazar envelope from Meyer et al. (2011)
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
;       Written by K. Willett                Oct 11
;-

bdir = '~/Astronomy/Research/blazars/meyer/'
mfile2 = bdir+'meyer2011_table2.txt'
mfile3 = bdir+'meyer2011_table3.txt'

;QSOB0003-06   00 06 13.89 -06 23 35.33 Q 0.347  13.16 46.22 41.68 44.21 VLA     f,h       

readcol,mfile2, $
	name2, rah, ram, ras, decd, dem, decs, $
	opt, z, lognu, loglp, logl300, loglkin, $
	/silent, $
	format='a,i,i,f,i,i,f,a,f,f,f,f,f,a,a'

readcol,mfile3, $
	name3, rah, ram, ras, decd, dem, decs, $
	morph, z, lognu_rg, loglp_rg, logl300_rg, loglkin_rg, $
	skipline=6, $
	/silent, $
	format='a,i,i,f,i,i,f,a,f,f,f,f,f,a,a'

loglkin    = 0.64 * (float(logl300)   - 40) + 43.54
loglkin_rg = 0.64 * (float(logl300_rg)- 40) + 43.54

b = where(opt eq 'B')
q = where(opt eq 'Q')
frI = where(morph eq '1')
frII = where(morph eq '2')

b_k = where(opt eq 'B' and loglkin lt 43.5)
b_b = where(opt eq 'B' and loglkin gt 43.5 and loglkin lt 44.5)
b_g = where(opt eq 'B' and loglkin gt 44.5 and loglkin lt 45.5)
b_o = where(opt eq 'B' and loglkin gt 45.5)

g_k = where(opt eq 'Q' and loglkin lt 43.5)
g_b = where(opt eq 'Q' and loglkin gt 43.5 and loglkin lt 44.5)
g_g = where(opt eq 'Q' and loglkin gt 44.5 and loglkin lt 45.5)
g_o = where(opt eq 'Q' and loglkin gt 45.5)

frI_k = where(morph eq '1' and loglkin_rg lt 43.5)
frI_b = where(morph eq '1' and loglkin_rg gt 43.5 and loglkin_rg lt 44.5)
frI_g = where(morph eq '1' and loglkin_rg gt 44.5 and loglkin_rg lt 45.5)
frI_o = where(morph eq '1' and loglkin_rg gt 45.5)

frII_k = where(morph eq '2' and loglkin_rg lt 43.5)
frII_b = where(morph eq '2' and loglkin_rg gt 43.5 and loglkin_rg lt 44.5)
frII_g = where(morph eq '2' and loglkin_rg gt 44.5 and loglkin_rg lt 45.5)
frII_o = where(morph eq '2' and loglkin_rg gt 45.5)

cgplot, lognu[b_k],loglp[b_k],$
	/xs,/ys,$
	psym=16,$
	color="Black", $
	xtitle='log '+greek('nu')+'!Ipeak!N [Hz]',$
	ytitle='log L!Ipeak!N [erg/s]',$
	xr=[12,18],$
	yr=[41,48] 

cgplot, lognu[b_b],loglp[b_b],$
	/over, $
	color="Blue", $
	psym=16

cgplot, lognu[b_g],loglp[b_g],$
	/over, $
	color="Forest Green", $
	psym=16

cgplot, lognu[b_o],loglp[b_o],$
	/over, $
	color="Orange", $
	psym=16

cgplot, lognu[g_b],loglp[g_b],$
	/over, $
	color="Black", $
	psym=17

cgplot, lognu[g_b],loglp[g_b],$
	/over, $
	color="Blue", $
	psym=17

cgplot, lognu[g_g],loglp[g_g],$
	/over, $
	color="Forest Green", $
	psym=17

cgplot, lognu[g_o],loglp[g_o],$
	/over, $
	color="Orange", $
	psym=17

; Radio galaxies

cgplot, lognu_rg[frI_k],loglp_rg[frI_k],$
	/over, $
	color="Black", $
	psym=15

cgplot, lognu_rg[frI_b],loglp_rg[frI_b],$
	/over, $
	color="Blue", $
	psym=15

cgplot, lognu_rg[frI_g],loglp_rg[frI_g],$
	/over, $
	color="Forest Green", $
	psym=15

;cgplot, lognu_rg[frI_o],loglp_rg[frI_o],$
;	/over, $
;	color="Orange", $
;	psym=15

cgplot, lognu_rg[frII_b],loglp_rg[frII_b],$
	/over, $
	color="Black", $
	psym=18

cgplot, lognu_rg[frII_b],loglp_rg[frII_b],$
	/over, $
	color="Blue", $
	psym=18

cgplot, lognu_rg[frII_g],loglp_rg[frII_g],$
	/over, $
	color="Forest Green", $
	psym=18

cgplot, lognu_rg[frII_o],loglp_rg[frII_o],$
	/over, $
	color="Orange", $
	psym=18

stop

end
