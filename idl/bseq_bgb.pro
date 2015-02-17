
;+
; NAME:
;       
;	BSEQ_BGB
;
; PURPOSE:
;
;	Look at the "blazar sequence" as function of B_gB for our sample
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

pro bseq_bgb, ps=ps, stop=stop

figdir='~/Astronomy/meetings/heidelberg2012/'
paperdir='~/Astronomy/Research/blazars/paper/'
b = mrdfits('~/Astronomy/Research/blazars/fits/bz_meyer.fits',1,/silent)

cs = 2

n = n_elements(b)

	if n_elements(lowz) eq 0 then lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	if n_elements(highz) eq 0 then highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	if n_elements(neighbors) eq 0 then neighbors = 10
	ind = where(b.z_2 lt highz and b.z_2 gt lowz and b.n500_cmag ge neighbors and b.bg gt 0)
	b = b[ind]

	bllac = where(b.btype eq 'BLLac' or b.btype eq 'Plotkin_blazar' or b.btype eq 'HBL' or b.btype eq 'lBLLac')
	fsrq = where(b.btype eq 'FSRQ')
	uncertain = where(b.btype eq 'BLLac_candidate' or b.btype eq 'blazar_uncertain' or b.btype eq 'Blazar')

	z = b.z_2
	fieldsize = b.fieldsize
	nt = float(b.n500_cmag)
	nb = float(b.bg)
	countingmag = b.cmag

	bgb, nt, nb, fieldsize, z, countingmag, bmeyer, bmeyer_err

bmax = max(bmeyer)
bmin = min(bmeyer)

if keyword_set(ps) then begin
	ps_start, filename=paperdir+'bgb_blazarsequence.eps', /quiet, /color, /encap, ysize=10, xsize=12
	cs = 2.5
	labelsize=1.3
	th = 5
	hsize=600
endif else begin
	cs = 2
	labelsize = 1.5
	hsize=10
endelse
	
cgplot, b.nupeak, b.lpeak, $
	position=[0.15,0.15,0.8,0.95], $
	charsize=cs, $
	thick=th, $
	/nodata, $
	yr=[43.95,47], /ystyle, $
	xr=[12,17], /xstyle, $
	title='Original Heidelberg plot', $
	xtitle='log ('+greek('nu')+'!Ipeak!N) [Hz]', $
	ytitle='log ('+greek('nu')+'L!I'+greek('nu')+'!N) [erg s!E-1!N]'

cgloadct, 13
minsize = 1.
for i=0,n_elements(bmeyer) - 1 do begin
	case strtrim(b[i].btype,2) of
		'BLLac': sym = 16
		'Plotkin_blazar': sym = 16
		'FSRQ': sym = 15
		else: sym=34
	endcase
	cgplot, b[i].nupeak, b[i].lpeak, $
		/over, $
		symsize=(bmeyer[i] - bmin + minsize) / (bmax-bmin) * 5. + 0, $
		color=fix((z[i] - min(z)) / (max(z)-min(z)) * 255.), $
		psym=sym
endfor

al_legend, /top,/left, psym=[16,15,34], ['BL Lac','FSRQ','Uncertain'], charsize=labelsize*1.4, symsize=2
cgcolorbar, position=[0.9,0.5,0.95,0.95], range=[min(z),max(z)], /vert, title='blazar redshift'

s1 = (-500 - bmin + minsize) / (bmax-bmin) * 5. + 0
s2 = (0 - bmin + minsize) / (bmax-bmin) * 5. + 0
s3 = (500 - bmin + minsize) / (bmax-bmin) * 5. + 0
s4 = (1000 - bmin + minsize) / (bmax-bmin) * 5. + 0
al_legend, position=[0.84,0.35], psym=16, ['-500','0','500','1000'], charsize=cs, symsize=[s1,s2,s3,s4],/normal,spacing = 3
cgtext, /normal, 0.895, 0.37, 'B!IgB!N',charsize=cs

; Overplot tracks from Meyer et al. (2011)

readcol, '/Applications/Dexter/meyer_sequence.gif.tracka', xa, ya, format='f,f', /skipline, /silent
readcol, '/Applications/Dexter/meyer_sequence.gif.trackb', xb, yb, format='f,f', /skipline, /silent

slope = (ya[1] - ya[0]) / (xa[1] - xa[0])
intercept = ya[0] - slope * xa[0]

;cgplots, xa, ya, color='black', linestyle=2
;cgplots, xb, yb, color='black', linestyle=2

; Kluges to make sure lines don't run off the edge of my truncated plot

cgarrow, (!Y.crange[1] - intercept) / slope, !y.crange[1],(!Y.crange[0] - intercept) / slope, !y.crange[0], color='black', linestyle=2, thick=th,/data,/solid,hsize=hsize
cgplots, [15.80,16.47,17.00], [44.24,44.50,44.82], color='black', linestyle=2, thick=th
cgarrow, 15.80,44.24, 15.05,44.00,color='black', linestyle=2, thick=th, /data, hsize=hsize,/solid

;cgtext, 13.7, 46.7, /data, 'Single-component jet', charsize=labelsize
;cgtext, 15.85, 44.2, /data, 'Decelerating jets', charsize=labelsize


if keyword_set(ps) then ps_end

; Split between LSP and HSP environments?

lsp = where(b.nupeak lt 15.)
hsp = where(b.nupeak gt 15.)
print,'B_gB LSP: ',mean(bmeyer[lsp]),stddev(bmeyer[lsp])
print,'B_gB HSP: ',mean(bmeyer[hsp]),stddev(bmeyer[hsp])
print,'redshift LSP: ',mean(z[lsp]),stddev(z[lsp])
print,'redshift HSP: ',mean(z[hsp]),stddev(z[hsp])
kstwo, bmeyer[lsp],bmeyer[hsp],d_bgb,prob_bgb
print,'K-S BgB, split by peak:', prob_bgb, probgauss(prob_bgb)		; Probability of B_gB being drawn from the same distribution is 0.3% (3 sigma)
kstwo, z[lsp],z[hsp],d_z,prob_z			; Probability of redshift being drawn from the same distribution is 2.4% (2.3 sigma)
print,'K-S redshift, split by peak:', prob_z, probgauss(prob_z)
kstwo, b[lsp].lpeak,b[hsp].lpeak,d_lpeak,prob_lpeak			; Probability of redshift being drawn from the same distribution is 2.4% (2.3 sigma)
print,'K-S L_peak, split by peak:', prob_lpeak, probgauss(prob_lpeak)
; Blazar type in each
;print,'High peaked'

if keyword_set(stop) then stop

end
