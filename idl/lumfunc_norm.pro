
;+
; NAME:
;       
;	LUMFUNC_NORM
;
; PURPOSE:
;
;	Plot the observed luminosity function for a field around blazar; find relative normalization of phi*
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
;       Written by K. Willett                Nov 11
;-

file = '~/Desktop/BZB3_willettk.csv'

readcol, file, objID,distance,u,r,g,i,z,redshift,redshifterr, $
	format='a,f,f,f,f,f,f,f,f', skipline=1, /quiet

z_blazar = 0.244

dl = lumdist(redshift,/silent) * 1d6

!p.multi=[0,2,3]
!p.font=-1
cs=3

cghistoplot, u, xtitle='m!Iu!N', charsize=cs & vline, mean(u), linestyle=2, /data, color=cgcolor("Blue")
cghistoplot, g, xtitle='m!Ig!N', charsize=cs & vline, mean(g), linestyle=2, /data, color=cgcolor("Blue") , /noerase
cghistoplot, r, xtitle='m!Ir!N', charsize=cs & vline, mean(r), linestyle=2, /data, color=cgcolor("Blue") , /noerase
cghistoplot, i, xtitle='m!Ii!N', charsize=cs & vline, mean(i), linestyle=2, /data, color=cgcolor("Blue") , /noerase
cghistoplot, z, xtitle='m!Iz!N', charsize=cs & vline, mean(z), linestyle=2, /data, color=cgcolor("Blue") , /noerase
cghistoplot, redshift, xtitle='Redshift', charsize=cs 
vline, mean(redshift), linestyle=2, /data, color=cgcolor("Blue"), /noerase

absmag_u = u - 5 * alog10(dl) + 5
absmag_g = g - 5 * alog10(dl) + 5
absmag_r = r - 5 * alog10(dl) + 5
absmag_i = i - 5 * alog10(dl) + 5
absmag_z = z - 5 * alog10(dl) + 5

cghistoplot, absmag_u, xtitle='M!Iu!N', charsize=cs, /oprob & vline, mean(absmag_u), linestyle=2, /data, color=cgcolor("Blue")
cghistoplot, absmag_g, xtitle='M!Ig!N', charsize=cs, /oprob & vline, mean(absmag_g), linestyle=2, /data, color=cgcolor("Blue") , /noerase
cghistoplot, absmag_r, xtitle='M!Ir!N', charsize=cs, /oprob & vline, mean(absmag_r), linestyle=2, /data, color=cgcolor("Blue") , /noerase
cghistoplot, absmag_i, xtitle='M!Ii!N', charsize=cs, /oprob & vline, mean(absmag_i), linestyle=2, /data, color=cgcolor("Blue") , /noerase
cghistoplot, absmag_z, xtitle='M!Iz!N', charsize=cs, /oprob & vline, mean(absmag_z), linestyle=2, /data, color=cgcolor("Blue") , /noerase
cghistoplot, redshift, xtitle='Redshift', charsize=cs & vline, mean(redshift), linestyle=2, /data, color=cgcolor("Blue"), /noerase

!p.multi=[0,1,1]

end
