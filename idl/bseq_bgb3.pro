;+
; NAME:
;       
;	BSEQ_BGB3
;
; PURPOSE:
;
;	Look at the "blazar sequence" as function of B_gB for our sample using new Eileen Meyer data (Sep 2012)
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

pro bseq_bgb3, ps=ps, stop=stop


paperdir='~/Astronomy/Research/blazars/paper/'
b = mrdfits('~/Astronomy/Research/blazars/fits/meyer_bz_allseds.fits',1,/silent)

cs = 2

n = n_elements(b)

	if n_elements(lowz) eq 0 then lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	if n_elements(highz) eq 0 then highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	if n_elements(neighbors) eq 0 then neighbors = 10

	ztemp = fltarr(n)
	for i = 0,n-1 do begin
		znum = strnumber(b[i].used_redshift,val)
		if znum eq 1 then ztemp[i] = val else ztemp[i] = b[i].z
	endfor

;	ind = where(ztemp lt highz and ztemp gt lowz and b.n500_cmag ge neighbors and b.bg gt 0 and (b.sed_code eq 'Uex' or b.sed_code eq 'Tex'))
	ind = where(ztemp gt lowz and ztemp lt highz and (b.sed_code eq 'Uex' or b.sed_code eq 'Tex'),nb)
	b = b[ind]

	b.btype = strtrim(b.btype,2)
	bllac = where(b.btype eq 'BLLac' or b.btype eq 'Plotkin_blazar' or b.btype eq 'HBL' or b.btype eq 'lBLLac')
	fsrq = where(b.btype eq 'FSRQ')
	uncertain = where(b.btype eq 'BLLac_candidate' or b.btype eq 'blazar_uncertain' or b.btype eq 'Blazar')

	z = ztemp[ind]
	fieldsize = b.fieldsize
	nt = float(b.n500_cmag)
	nb = float(b.bg)
	countingmag = b.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bmeyer, bmeyer_err

bmin = min(bmeyer)
bmax = max(bmeyer)
	bmin = -500.
	bmax = 1500.

!p.multi=[0,3,1]

if keyword_set(ps) then begin
	ps_start, filename=paperdir+'bgb_blazarsequence_allseds.eps', /quiet, /color, /encap, ysize=10, xsize=12
	cs = 2.5
	ct = 1
	labelsize=1.3
	th = 5
	hsize=600
endif else begin
	cs = 4
	ct = 2
	labelsize = 1.5
	hsize=10
endelse
	
cgloadct, 13
minsize = 0.5
sizescale = 5.

; Overplot tracks from Meyer et al. (2011)

readcol, '/Applications/Dexter/meyer_sequence.gif.tracka', xa, ya, format='f,f', /skipline, /silent
readcol, '/Applications/Dexter/meyer_sequence.gif.trackb', xb, yb, format='f,f', /skipline, /silent

slope = (ya[1] - ya[0]) / (xa[1] - xa[0])
intercept = ya[0] - slope * xa[0]

;cgplots, xa, ya, color='black', linestyle=2
;cgplots, xb, yb, color='black', linestyle=2

cgplot, b.nupeak, b.lpeak, $
	background='white', $
	position=[0.08,0.15,0.30,0.95], $
	charsize=cs, $
	charthick = ct, $
	thick=th, $
	/nodata, $
	yr=[43.95,47], /ystyle, $
	xr=[12,17], /xstyle, $
	title='BL Lacs', $
	xtitle='log ('+greek('nu')+'!Ipeak!N) [Hz]', $
	ytitle='log ('+greek('nu')+'L!I'+greek('nu')+'!N) [erg s!E-1!N]'

for j=0,n_elements(bllac) - 1 do begin
	cgplot, b[bllac[j]].nupeak, b[bllac[j]].lpeak, $
		/over, $
		symsize=(bmeyer[bllac[j]] - bmin) / (bmax-bmin) * sizescale + minsize, $
		color=fix((z[bllac[j]] - min(z)) / (max(z)-min(z)) * 255.), $
		psym=16
endfor

cgarrow, (!Y.crange[1] - intercept) / slope, !y.crange[1],(!Y.crange[0] - intercept) / slope, !y.crange[0], color='black', linestyle=2, thick=th,/data,/solid,hsize=hsize
cgplots, [15.80,16.47,17.00], [44.24,44.50,44.82], color='black', linestyle=2, thick=th
cgarrow, 15.80,44.24, 15.05,44.00,color='black', linestyle=2, thick=th, /data, hsize=hsize,/solid

cgplot, b.nupeak, b.lpeak, $
	background='white', $
	position=[0.37,0.15,0.58,0.95], $
	charsize=cs, $
	charthick = ct, $
	thick=th, $
	/nodata, $
	yr=[43.95,47], /ystyle, $
	xr=[12,17], /xstyle, $
	title='FSRQs', $
	xtitle='log ('+greek('nu')+'!Ipeak!N) [Hz]'

for j=0,n_elements(fsrq) - 1 do begin
	cgplot, b[fsrq[j]].nupeak, b[fsrq[j]].lpeak, $
		/over, $
		symsize=(bmeyer[fsrq[j]] - bmin) / (bmax-bmin) * sizescale + minsize, $
		color=fix((z[fsrq[j]] - min(z)) / (max(z)-min(z)) * 255.), $
		psym=15
endfor

cgarrow, (!Y.crange[1] - intercept) / slope, !y.crange[1],(!Y.crange[0] - intercept) / slope, !y.crange[0], color='black', linestyle=2, thick=th,/data,/solid,hsize=hsize
cgplots, [15.80,16.47,17.00], [44.24,44.50,44.82], color='black', linestyle=2, thick=th
cgarrow, 15.80,44.24, 15.05,44.00,color='black', linestyle=2, thick=th, /data, hsize=hsize,/solid

cgplot, b.nupeak, b.lpeak, $
	background='white', $
	position=[0.62,0.15,0.80,0.95], $
	charsize=cs, $
	charthick = ct, $
	thick=th, $
	/nodata, $
	yr=[43.95,47], /ystyle, $
	xr=[12,17], /xstyle, $
	title='Uncertain blazars', $
	xtitle='log ('+greek('nu')+'!Ipeak!N) [Hz]'

for j=0,n_elements(uncertain) - 1 do begin
	cgplot, b[uncertain[j]].nupeak, b[uncertain[j]].lpeak, $
		/over, $
		symsize=(bmeyer[uncertain[j]] - bmin) / (bmax-bmin) * sizescale + minsize, $
		color=fix((z[uncertain[j]] - min(z)) / (max(z)-min(z)) * 255.), $
		psym=14
endfor



;al_legend, /top,/left, psym=[16,15,34], ['BL Lac','FSRQ','Uncertain'], charsize=labelsize*1.4, symsize=2, outline_color='black', textcolor='black', colors='black'
cgcolorbar, position=[0.9,0.5,0.95,0.95], range=[min(z),max(z)], /vert, title='blazar redshift', color='black'

s1 = (-500 - bmin) / (bmax-bmin) * sizescale + minsize
s2 = (0 - bmin) / (bmax-bmin) * sizescale + minsize
s3 = (500 - bmin) / (bmax-bmin) * sizescale + minsize
s4 = (1000 - bmin) / (bmax-bmin) * sizescale + minsize
al_legend, position=[0.84,0.35], psym=16, ['-500','0','500','1000'], charsize=cs, symsize=[s1,s2,s3,s4],/normal,spacing = 3, outline_color='black', textcolor='black', colors='black'
cgtext, /normal, 0.895, 0.37, 'B!IgB!N',charsize=cs

; Kluges to make sure lines don't run off the edge of my truncated plot

cgarrow, (!Y.crange[1] - intercept) / slope, !y.crange[1],(!Y.crange[0] - intercept) / slope, !y.crange[0], color='black', linestyle=2, thick=th,/data,/solid,hsize=hsize
cgplots, [15.80,16.47,17.00], [44.24,44.50,44.82], color='black', linestyle=2, thick=th
cgarrow, 15.80,44.24, 15.05,44.00,color='black', linestyle=2, thick=th, /data, hsize=hsize,/solid

;cgtext, 13.7, 46.7, /data, 'Single-component jet', charsize=labelsize
;cgtext, 15.85, 44.2, /data, 'Decelerating jets', charsize=labelsize


if keyword_set(ps) then ps_end

; Split between low and high-peaked environments?

lowpeak = where(b.nupeak lt 14.5)
highpeak = where(b.nupeak gt 14.5)
highpower = where(b.lpeak gt 45.0)
lowpower = where(b.lpeak le 45.0)

print,''
print,'B_gB lowpeak: ',mean(bmeyer[lowpeak]),stddev(bmeyer[lowpeak]), median(bmeyer[lowpeak])
print,'B_gB highpeak: ',mean(bmeyer[highpeak]),stddev(bmeyer[highpeak]), median(bmeyer[highpeak])
print,''
print,'redshift lowpeak: ',mean(z[lowpeak]),stddev(z[lowpeak])
print,'redshift highpeak: ',mean(z[highpeak]),stddev(z[highpeak])
print,''
kstwo, bmeyer[lowpeak],bmeyer[highpeak],d_bgb,prob_bgb
print,'K-S BgB, split by peak:', prob_bgb, probgauss(prob_bgb)		; Probability of B_gB being drawn from the same distribution is 0.3% (3 sigma)
kstwo, z[lowpeak],z[highpeak],d_z,prob_z			; Probability of redshift being drawn from the same distribution is 2.4% (2.3 sigma)
print,'K-S redshift, split by peak:', prob_z, probgauss(prob_z)
kstwo, b[lowpeak].lpeak,b[highpeak].lpeak,d_lpeak,prob_lpeak			; Probability of redshift being drawn from the same distribution is 2.4% (2.3 sigma)
print,'K-S L_peak, split by peak:', prob_lpeak, probgauss(prob_lpeak)
kstwo, bmeyer[lowpower],bmeyer[highpower],d_lpower,prob_lpower			; Probability of redshift being drawn from the same distribution is 2.4% (2.3 sigma)
print,'K-S L_peak, split by power:', prob_lpower, probgauss(prob_lpower)
; Blazar type in each
;print,'High peaked'
print,''

;stop
;
;; Correlations
;
;mbh_flt = float(b.mass_bh)
;mbhind = where(mbh_flt ne 0.)
;lext_flt = float(b.lext)
;lextind = where(lext_flt ne 0.)
;
;!p.multi=[0,2,2]
;cgplot,bmeyer,b.lpeak,psym=16,xr=[-1000,1500],yr=[42,47],title='L!Ipeak!N', charsize=2
;cgplot,bmeyer,b.nupeak,psym=16,xr=[-1000,1500],yr=[12,17],title=greek('nu')+'!Ipeak!N', charsize=2
;cgplot,bmeyer,b.mass_bh,psym=16,xr=[-1000,1500],yr=[6,10],title='M!IBH!N', charsize=2
;cgplot,bmeyer,b.lext,psym=16,xr=[-1000,1500],yr=[37,44],title='L!Iext!N', charsize=2
;
;print,'Bgb vs lpeak:   ',string(correlate(bmeyer,float(b.lpeak)),format='(f7.4)')
;print,'Bgb vs nupeak:  ',string(correlate(bmeyer,float(b.nupeak)),format='(f7.4)')
;print,'Bgb vs BH mass: ',string(correlate(bmeyer[mbhind],mbh_flt[mbhind]),format='(f7.4)')
;print,'Bgb vs L_ext:   ',string(correlate(bmeyer[lextind],lext_flt[lextind]),format='(f7.4)')


if keyword_set(stop) then stop

end


