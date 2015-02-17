
;+
; NAME:
;       
;	KIMBALL_SPECGAL
;
; PURPOSE:
;
;	Find SDSS-detected radio galaxies with RA, dec and classified by radio morphology
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

all4radio_spec = '~/Astronomy/Research/blazars/kimball_specgal.csv'

readcol, all4radio_spec, $
	kimball_id, ra, dec, first_peak_flux, first_flux, nvss_flux, near_modelmag_r, spec_type, spec_redshift, spec_redshifterr, spec_redshift_warning, $
	/silent, $
	skipline=1, $
	delimiter=',', $
	format='a,f,f,f,f,f,f,i,f,f,l'

ngals = n_elements(kimball_id)

; Compute the radio flux "AB magnitude"

first_flux = float(first_flux)
first_peak_flux = float(first_peak_flux)
nvss_flux = float(nvss_flux)

tfirst = -2.5 * alog10(first_flux * 1d-3 /3631.)
tnvss = -2.5 * alog10(nvss_flux * 1d-3 /3631.)

deltat = tfirst - tnvss

theta = sqrt(first_flux / first_peak_flux)

!p.multi=[0,2,2]

cs = 1

cghistoplot, deltat, xtitle='t!IFIRST!N - t!INVSS!N', xr=[-1,2], charsize=cs
cgplots, [0.35,0.35], !y.crange, linestyle=2, /data, thick=2
cgtext, 1.0, 400, 'Complex', /data
cgtext, -0.8, 400, 'Simple',charsize=2.0, /data
cgtext, -0.8, 400-50, '(Compact + Resolved)',charsize=1.0, /data

cghistoplot, alog10(theta^2), xtitle='log '+greek('theta')+'!E2!N', xr=[-0.20,0.50], charsize=cs
cgplots, [0.05,0.05], !y.crange, linestyle=2, /data, thick=2
cgtext, 0.2, 400, 'Resolved', /data
cgtext, -0.15, 400, 'Compact', /data

ind_complex = where(deltat gt 0.35, ncomplex)
ind_compact = where(deltat lt 0.35 and alog10(theta^2) lt 0.05, ncompact)
ind_resolved = where(deltat lt 0.35 and alog10(theta^2) gt 0.05, nresolved)

delz = 0.02

cghistoplot, spec_redshift, $
	xr=[0,1], $
	bin=delz, $
	charsize=cs, $
	xtitle='Redshift'

print,''
;print,'Number of complex  sources: ', ncomplex
;print,'Number of compact  sources: ', ncompact
;print,'Number of resolved sources: ', nresolved

; Trim KI galaxies based on the blazar redshift distribution

	blazardir = '~/Astronomy/Research/blazars/'
	restore,blazardir+'ab_structure.sav'

	lowz = 0.043	; Redshift at which angular size of 10 arcmin = 500 kpc
	highz = 0.75	; Redshift at which an M* galaxy reaches the SDSS completeness limit
	neighbors = 10
	ind = where(bz.z lt highz and bz.z gt lowz and bz.n500_cmag ge neighbors and bz.bg gt 0)
	bz=bz[ind]

bins = fillarr(delz,lowz,highz)
result = histogram(bz.z,bins=delz,omin=omin,omax=omax)

;cghistoplot,bz.z,bin=delz, datacolor='blue',/oplot
cgplot,bins+delz/2.,result,color='blue',/over,psym=10
cgtext,0.7,150,'Spec. gals',color='red', charsize=cs
cgtext,0.7,130,'Blazars',color='blue', charsize=cs

specind = where(spec_redshift lt highz and spec_redshift gt lowz)
kstwo, bz.z, spec_redshift[specind], d, prob
;print,'KS result:',d,prob,sqrt(2d) * inverf(1d - prob)


; For each bin in the MODIFIED blazar redshift distribution, randomly select the same number of KI specgals to fill it.
; Necessary since at z>0.5, there are not more spec gal targets than blazars in each bin. Distribution could be matched IF
; I decreased the total number of galaxies in the sample selection. That would mean scaling to a few bins at the outer edge, though.  

breakz = 0.162d

binarr10 = fillarr(delz,lowz,breakz)
totlow = where(bz.z le breakz,nlow)
n10 = n_elements(binarr10)
ki_10arcmin_ind = intarr(nlow)
j=0

for i=0, n10-2 do begin
	bjunk = where(bz.z gt binarr10[i] and bz.z le binarr10[i+1],nb10)
	kjunk = where(spec_redshift gt binarr10[i] and spec_redshift le binarr10[i+1],nk10)
	if (nb10 gt 0) and (nk10 gt 0) then ki_10arcmin_ind[j:j+nb10-1] = kjunk[cgRandomIndices(nk10,nb10,seed=seed1)]
	j+=nb10
endfor
bjunk_10last = where(bz.z gt binarr10[n10-1] and bz.z le breakz,nb10_last)
kjunk_10last = where(spec_redshift gt binarr10[n10-1] and spec_redshift le breakz,nk10_last)
if (nb10_last gt 0) and (nk10_last gt 0) then ki_10arcmin_ind[j:j+nb10_last-1] = kjunk_10last[cgRandomIndices(nk10_last,nb10_last,seed=seed2)]

binarr3 = fillarr(delz,breakz,max(bz.z))
tothigh = where(bz.z gt breakz,nhigh)
n3 = n_elements(binarr3)
ki_3arcmin_ind = intarr(nhigh)
j=0

for i=0, n3-2 do begin
	bjunk = where(bz.z gt binarr3[i] and bz.z le binarr3[i+1],nb3)
	kjunk = where(spec_redshift gt binarr3[i] and spec_redshift le binarr3[i+1],nk3)
	if (nb3 gt 0) and (nk3 gt 0) then ki_3arcmin_ind[j:j+nb3-1] = kjunk[cgRandomIndices(nk3,nb3,seed=seed3)] else print,nb01,nk3
	j+=nb3
endfor
bjunk_3last = where(bz.z gt binarr3[n3-1] and bz.z le max(bz.z),nb3_last)
kjunk_3last = where(spec_redshift gt binarr3[n3-1] and spec_redshift lt max(bz.z),nk3_last)
if (nb3_last gt 0) and (nk3_last gt 0) then ki_3arcmin_ind[j:j+nb3_last-1] = kjunk[cgRandomIndices(nk3_last,nb3_last,seed=seed4)]

ki_3arcmin_ind = ki_3arcmin_ind[where(ki_3arcmin_ind ne 0)]	; Weird kluge; unhappy about this. 

cgplot,spec_redshift[ki_10arcmin_ind],charsize=cs,yr=[0,1],xr=[0,nlow+nhigh]
cgplot,indgen(nhigh)+nlow,spec_redshift[ki_3arcmin_ind],/over,color='green'
!p.multi=[0,1,1]

; Write a CSV file on which SDSS can search for photometry and nearby neighbors

kimball_rmorph = strarr(ngals)
kimball_rmorph[ind_complex] = 'COMPLEX'
kimball_rmorph[ind_compact] = 'COMPACT'
kimball_rmorph[ind_resolved] = 'RESOLVED'

kimballarr3 = [transpose(string(kimball_id[ki_3arcmin_ind],format='(a14)')), transpose(string(ra[ki_3arcmin_ind])), transpose(string(dec[ki_3arcmin_ind])), transpose(string(spec_redshift[ki_3arcmin_ind])),transpose(kimball_rmorph[ki_3arcmin_ind])]
kimballarr10 = [transpose(string(kimball_id[ki_10arcmin_ind],format='(a14)')), transpose(string(ra[ki_10arcmin_ind])), transpose(string(dec[ki_10arcmin_ind])), transpose(string(spec_redshift[ki_10arcmin_ind])),transpose(kimball_rmorph[ki_10arcmin_ind])]

;openw, lun1, '~/Astronomy/Research/blazars/kimball_radiogals.cat', /get_lun

openw, lun1, '~/Astronomy/Research/blazars/ki_modz_3arcmin.cat', /get_lun
printf, lun1, 'kimball_id              ra              dec            z           radio_morph'
printf, lun1, kimballarr3
close, lun1
free_lun, lun1

openw, lun1, '~/Astronomy/Research/blazars/ki_modz_10arcmin.cat', /get_lun
printf, lun1, 'kimball_id              ra              dec            z           radio_morph'
printf, lun1, kimballarr10
close, lun1
free_lun, lun1

; Make control field samples for both the 3 and 10 arcminute Kimball radio galaxies

ra3 = ra[ki_3arcmin_ind]
dec3 = dec[ki_3arcmin_ind]

ranew3_north = ra3
ranew3_south = ra3
decnew3_north = dec3
decnew3_south = dec3
highdec = where(dec3 ge 89., hdcount)
lowdec = where(dec3 lt 89., ldcount)
if hdcount gt 0 then begin
	ranew3_north[highdec] = raflip(ra3[highdec],/deg)
	decnew3_south[highdec] = dec3 - 1.
endif
if ldcount gt 0 then begin
	decnew3_north[lowdec] = dec3 + 1.
	decnew3_south[lowdec] = dec3 - 1.
endif

kimballarr3_connorth = [transpose(string(kimball_id[ki_3arcmin_ind],format='(a14)')), transpose(string(ranew3_north)), transpose(string(decnew3_north)), transpose(string(spec_redshift[ki_3arcmin_ind])),transpose(kimball_rmorph[ki_3arcmin_ind])]
kimballarr3_consouth = [transpose(string(kimball_id[ki_3arcmin_ind],format='(a14)')), transpose(string(ranew3_north)), transpose(string(decnew3_south)), transpose(string(spec_redshift[ki_3arcmin_ind])),transpose(kimball_rmorph[ki_3arcmin_ind])]

openw, lun1, '~/Astronomy/Research/blazars/ki_modz_3arcmin_connorth.cat', /get_lun
printf, lun1, 'kimball_id              ra              dec            z           radio_morph'
printf, lun1, kimballarr3_connorth
close, lun1
free_lun, lun1

openw, lun1, '~/Astronomy/Research/blazars/ki_modz_3arcmin_consouth.cat', /get_lun
printf, lun1, 'kimball_id              ra              dec            z           radio_morph'
printf, lun1, kimballarr3_consouth
close, lun1
free_lun, lun1

ra10 = ra[ki_10arcmin_ind]
dec10 = dec[ki_10arcmin_ind]

ranew10_north = ra10
ranew10_south = ra10
decnew10_north = dec10
decnew10_south = dec10
highdec = where(dec10 ge 89., hdcount)
lowdec = where(dec10 lt 89., ldcount)
if hdcount gt 0 then begin
	ranew10_north[highdec] = raflip(ra10[highdec],/deg)
	decnew10_south[highdec] = dec10 - 1.
endif
if ldcount gt 0 then begin
	decnew10_north[lowdec] = dec10 + 1.
	decnew10_south[lowdec] = dec10 - 1.
endif

kimballarr10_connorth = [transpose(string(kimball_id[ki_10arcmin_ind],format='(a14)')), transpose(string(ranew10_north)), transpose(string(decnew10_north)), transpose(string(spec_redshift[ki_10arcmin_ind])),transpose(kimball_rmorph[ki_10arcmin_ind])]
kimballarr10_consouth = [transpose(string(kimball_id[ki_10arcmin_ind],format='(a14)')), transpose(string(ranew10_north)), transpose(string(decnew10_south)), transpose(string(spec_redshift[ki_10arcmin_ind])),transpose(kimball_rmorph[ki_10arcmin_ind])]

openw, lun1, '~/Astronomy/Research/blazars/ki_modz_10arcmin_connorth.cat', /get_lun
printf, lun1, 'kimball_id              ra              dec            z           radio_morph'
printf, lun1, kimballarr10_connorth
close, lun1
free_lun, lun1

openw, lun1, '~/Astronomy/Research/blazars/ki_modz_10arcmin_consouth.cat', /get_lun
printf, lun1, 'kimball_id              ra              dec            z           radio_morph'
printf, lun1, kimballarr10_consouth
close, lun1
free_lun, lun1


print,''

; Plot positions in an Aitoff projection

!p.multi=[0,1,1]
aitoff_grid, label=2, /new			; Can't change the center coordinates of map 
aitoff, ra, dec, x, y
aitoff, ranew10_north, decnew10_north, x10n, y10n
aitoff, ranew10_south, decnew10_south, x10s, y10s
cgplots, x, y, psym=3, color="Yellow"
cgplots, x10n, y10n, psym=3, color="Red"
cgplots, x10s, y10s, psym=3, color="Green"

end
