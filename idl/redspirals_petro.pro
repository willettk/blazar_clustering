;+
; NAME:
;       
;	REDSPIRALS_PETRO
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

file='~/Astronomy/proposals/redspirals/redspirals_petro_willettk.csv'

readcol, file, $
	dr7objid,dr8objid,petroR50_r,petroR90_r,z, $
	format='a,a,f,f,f', $
	comment='%', $
	/silent, $
	skip=1
	
cghistoplot, petroR90_r, $
	xtitle='R!Ipetro!N', $
	datacolor='red', $
	charsize=1.5

cghistoplot, petroR50_r, $
	/oplot, $
	datacolor='blue'

sort50 = reverse(petroR50_r[sort(petroR50_r)])
match,sort50[0:49],petror50_r,a,b
;print,[transpose(dr7objid[b]),transpose(string(petror50_r[b],format='(f4.1)'))]

d=[9.4	, 11.3	, 20.4	, 8.7	, 25.3	, 7.8	, 9.7	, 12.3	, 16.8	, 20.5	, 17.9	, 5.2	, 5.4	, 6.1	, 5.5]	
cts=[13,   50     , 73     , 133    , 571    , 1575   , 27     , 100    , 12     , 12     , 10     , 379    , 20     , 13     , 766  ]
lx=(10d)^([38.3,	38.4,	38.6,	37.9, 39.3, 38.9, 38.5, 39.2, 38.4, 38.3, 38.4, 37.9, 37.3, 37.1, 38.2])
t = [2.2	, 10.0	, 26.6	, 42.1	, 67.1	, 39.6	, 3.0	, 3.5	, 4.9	, 9.0	, 4.9	, 49.1	, 9.6	, 13.2	, 58.3	]

zarr = fillarr(0.001,0.001,0.09)
fac = mean(t*lx/(cts * d^2))
dl = lumdist(zarr,/wmap7,/silent)
d0 = lumdist(min(z),/wmap7,/silent)
d1 = lumdist(0.05,/wmap7,/silent)
d2 = lumdist(0.085,/wmap7,/silent)

!p.multi=[0,2,1]

c0 = 25.
l0 = 1d40

cgplot,zarr,fac * c0 / l0 * dl^2, xtitle='z',ytitle='Time [ksec]', charsize=1.5
cgplots,[min(z),min(z)],!y.crange,linestyle=2
cgplots,[0.05,0.05],!y.crange,linestyle=2
cgplots,[0.085,0.085],!y.crange,linestyle=2

fac2 = fac * c0/l0
print,fac2

totaltime = 0.
sortzind = sort(z)
zind = 0

while totaltime lt 300 do begin
	newz = z[sortzind[zind]]
	newtime = fac * c0 / l0 * (lumdist(newz,/wmap7,/silent))^2
	totaltime += newtime
	zind += 1
endwhile

newinds = sortzind[0:zind-1]

ctime = totaltime - newtime
ngal = zind
prad = petroR50_r[newinds]

cgplot,z[newinds],prad, $
	psym=7,charsize = 1.5, $
	xtitle='z',$
	ytitle='Petrosian radius [arcsec]'

print,''
print,'Total exposure time [ksec]: ',ctime
print,'Total number of galaxies: ',ngal
range,prad
print, transpose(dr7objid[newinds])

print,''

end
