bsfile = '~/Astronomy/proposals/redspirals_vla/bluespirals_lines_willettk.fit'

bs = mrdfits(bsfile,1,/silent)

;, sdss_lines_file, $
;	dr7objid,dr8objid,$
;	petroR50_r,petroR90_r,$
;	z, $
;	nii_6584_flux,nii_6584_flux_err,$
;	oiii_5007_flux,oiii_5007_flux_err,$
;	neiii_3869_flux,neiii_3869_flux_err,$
;	oii_3726_flux,oii_3726_flux_err,$
;	oii_3729_flux,oii_3729_flux_err, $
;	h_beta_flux,h_beta_flux_err,$
;	h_alpha_flux,h_alpha_flux_err, $
;	modelmag_g, modelmag_z, $
;;	kcorr_g, kcorr_z, $
;	skip=1, $
;	/silent, $
;	format='a,a,f,f,f, f,f,f,f,f,f,f,f,f,f,f,f,f,f, f,f,f,f'

; Look only at the 17 lowest-redshift galaxies that we can examine with Chandra

	;lowzids =[ '587733603730522222', '587738568705507359', '587741488761274453', '587738410327801941', '587736543096734014', '587741602564014184', '587734949666881571', '587725551199780967', '588017978911817806', '587736783606382788', '587738066727534662', '587737827288875214', '587739630096023592', '587733397573009582', '587722982822379684', '587730772799914495', '587734894361706638']
	;
	;match, lowzids, dr7objid, ind1, ind2, count=mcount
	;
	;dr7objid           =dr7objid[ind2]
	;dr8objid           =dr8objid[ind2]
	;petroR50_r         =petroR50_r[ind2]
	;petroR90_r         =petroR90_r[ind2]
	;z                  =z[ind2]
	;nii_6584_flux      =nii_6584_flux[ind2]
	;nii_6584_flux_err  =nii_6584_flux_err[ind2]
	;oiii_5007_flux     =oiii_5007_flux[ind2]
	;oiii_5007_flux_err =oiii_5007_flux_err[ind2]
	;neiii_3869_flux    =neiii_3869_flux[ind2]
	;neiii_3869_flux_err=neiii_3869_flux_err[ind2]
	;oii_3726_flux      =oii_3726_flux[ind2]
	;oii_3726_flux_err  =oii_3726_flux_err[ind2]
	;oii_3729_flux      =oii_3729_flux[ind2]
	;oii_3729_flux_err  =oii_3729_flux_err[ind2]
	;h_beta_flux        =h_beta_flux[ind2]
	;h_beta_flux_err    =h_beta_flux_err[ind2]
	;h_alpha_flux       =h_alpha_flux[ind2]
	;h_alpha_flux_err   =h_alpha_flux_err[ind2]
	;modelmag_g         =modelmag_g[ind2]
	;modelmag_z         =modelmag_z[ind2]

;	noneiii = where(neiii_3869_flux lt 0.)
;	neiii_3869_flux[noneiii] = neiii_3869_flux_err[noneiii]

;sdssid=strtrim(sdssid,2)
;match, dr6id, sdssid, drind, sdssind, count=mcount
;bptind = drind[uniq(dr6id[drind])]	; Proud of this nesting.
;neduind = sdssind[uniq(sdssid[sdssind])]

nsig = 3.
nsig_neiii = 1.
goodind = where($
;	(neiii_3869_flux gt nsig_neiii * neiii_3869_flux_err) and $
	;((oii_3726_flux + oii_3729_flux) gt nsig * sqrt(oii_3726_flux_err^2 + oii_3729_flux_err^2)) and $
	(bs.nii_6584_flux gt nsig * bs.nii_6584_flux_err) and $
	(bs.oiii_5007_flux gt nsig * bs.oiii_5007_flux_err) and $
	(bs.h_alpha_flux gt nsig * bs.h_alpha_flux_err) and $
	(bs.h_beta_flux gt nsig * bs.h_beta_flux_err), $
	goodcount)

bsgood = bs[goodind]

print,''

nha = alog10(bsgood.nii_6584_flux/bs.h_alpha_flux)
ohb = alog10(bsgood.oiii_5007_flux/bs.h_beta_flux)

dir = '~/Astronomy/proposals/redspirals_chandra/'
if keyword_set(ps) then begin
	ps_start, filename=dir+'masters_bluespirals.eps',/color,/landscape, xs=10, ys=5, /encap
	cs = 1.0
	th = 5
	hs = 150
endif else begin
	th = 1
	cs = 1.5
	hs = 10
endelse

!p.multi=[0,2,1]
cgplot, nha, ohb, $
	thick=th, xthick=th, ythick=th, $
	/nodata, $
	xr=[-1.5,0.5], /xstyle, $
	yr=[-1.0,1.5], /ystyle, $
	title='BPT', $
	charsize=cs, $
	xtitle='log ([NII] '+greek('lambda')+'6584/H'+greek('alpha')+')', $
	ytitle='log ([OIII] '+greek('lambda')+'5007/H'+greek('beta')+')', $
	psym=9

xkauffmann=fillarr(0.01,-10,0)
kauffmann = 0.61 /(xkauffmann - 0.05)+1.3

xkewley=fillarr(0.01,-10,0.4)
kewley = 0.61 /(xkewley - 0.47)+1.19

xtrouille = fillarr(0.01,-10,10)
trouille = -1.2 * xtrouille - 0.4

cgplot, xkauffmann, kauffmann, $
	color='blue', $
	/over, $
	linestyle=2, thick=4

cgplot, xkewley, kewley, $
	color='red', $
	/over, $
	linestyle=0, thick=4

x1=-0.2
x2=0.5
y1=0.61 /(-0.2 - 0.47)+1.19
y2=0.98
cgplots, [x1,x2], [y1,y2], linestyle=0		; LINER

xliner = fillarr(0.01,-1,1)

linerslope = (y2-y1)/(x2-x1)
linerint = y1 - (linerslope*x1)

sf = where(ohb lt (0.61 /(nha - 0.05)+1.3) and (nha lt -0.2), scount)
agn = where(ohb gt (0.61 /(nha - 0.47)+1.19) and (ohb gt (linerslope * nha + linerint)), acount)
comp = where((ohb lt (0.61 /(nha - 0.47)+1.19)) and (ohb gt (0.61 /(nha - 0.05)+1.3)), ccount)
liner = where((ohb gt (0.61 /(nha - 0.47)+1.19)) and (ohb lt (linerslope * nha + linerint)), lcount)

if scount gt 0 then cgplot, /overplot, nha[sf], ohb[sf], color='purple', psym=4
if acount gt 0 then cgplot, /overplot, nha[agn], ohb[agn], color='black', psym=1
if ccount gt 0 then cgplot, /overplot, nha[comp], ohb[comp], color='green', psym=17
if lcount gt 0 then cgplot, /overplot, nha[liner], ohb[liner], color='orange', psym=16

cgtext, /data, -0.7, 1.0, 'AGN', charsize=cs
cgtext, /data, 0.15, -0.2, 'LINERs', charsize=cs
cgtext, /data, -1.3, 0.0, 'Star-forming', charsize=cs
cgtext, /data, -0.2, -0.7, 'Comp.', charsize=cs

neo = alog10(bsgood.neiii_3869_flux / (bsgood.oii_3726_flux + bsgood.oii_3729_flux))
gz =  (bsgood.modelmag_g) - (bsgood.modelmag_z)

nsig_ne3 = 3.
badne = where(bsgood.neiii_3869_flux lt nsig_ne3 * bsgood.neiii_3869_flux_err,bnc)
goodne = where(bsgood.neiii_3869_flux gt nsig_ne3 * bsgood.neiii_3869_flux_err,gnc)
;badne = where(neiii_3869_flux lt nsig_ne3 * neiii_3869_flux_err,bnc)

cgplot, neo[goodne], gz[goodne], $
	/nodata, $
	thick=th, xthick=th, ythick=th, $
	xr=[-3.0,0.5], /xstyle, $
	yr=[1.0,2.2], /ystyle, $
	title='TBT', $
	charsize=cs, $
	ytitle='(g - z)', $
	xtitle='log ([NeIII '+greek('lambda')+'3869]/[OII '+greek('lambda')+greek('lambda')+'3726,3729])', $
	psym=9

;if scount gt 0 then cgplot, /overplot, neo[setintersection(sf,goodne)], gz[setintersection(sf,goodne)], color='purple', psym=4, symsize=1.5
if acount gt 0 then cgplot, /overplot, neo[setintersection(agn,goodne)], gz[setintersection(agn,goodne)], color='black', psym=1, symsize=1.0
if ccount gt 0 then cgplot, /overplot, neo[setintersection(comp,goodne)], gz[setintersection(comp,goodne)], color='green', psym=17, symsize=1.0
if lcount gt 0 then cgplot, /overplot, neo[setintersection(liner,goodne)], gz[setintersection(liner,goodne)], color='orange', psym=16, symsize=1.0

cgplot, /overplot, xtrouille, trouille, thick=4

cgtext, /data, -1.3, 2.0, 'AGN', charsize=cs
cgtext, /data, -2.5, 1.2, 'Star-forming', charsize=cs

; Limits on [NeIII]

	;if bnc gt 0 then cgarrow, neo[badne], gz[badne], neo[badne] - 0.1, gz[badne], /data, color=' grey', hsize=hs

; Print list of ra, dec, spectral type

openw, lun, '~/Astronomy/proposals/redspirals_vla/bs_bytype.txt',/get_lun,width=120
printf, lun, 'ra, dec, redshift, objid, spectype,'
printf, lun, [transpose(string(bsgood[liner].ra)), transpose(string(bsgood[liner].dec)),transpose(string(bsgood[liner].redshift)),transpose(string(bsgood[liner].dr7objid)),transpose(replicate('LINER',n_elements(liner)))]+','
printf, lun, [transpose(string(bsgood[sf].ra)), transpose(string(bsgood[sf].dec)),transpose(string(bsgood[sf].redshift)),transpose(string(bsgood[sf].dr7objid)),transpose(replicate('SF',n_elements(sf)))]+','
printf, lun, [transpose(string(bsgood[agn].ra)), transpose(string(bsgood[agn].dec)),transpose(string(bsgood[agn].redshift)),transpose(string(bsgood[agn].dr7objid)),transpose(replicate('AGN',n_elements(agn)))]+','
printf, lun, [transpose(string(bsgood[comp].ra)), transpose(string(bsgood[comp].dec)),transpose(string(bsgood[comp].redshift)),transpose(string(bsgood[comp].dr7objid)),transpose(replicate('Composite',n_elements(comp)))]+','
close, lun

; Restore matched files from TOPCAT

bsradio = mrdfits('~/Astronomy/proposals/redspirals_vla/bs_radio.fits',1,/silent)

rliner = where(strtrim(bsradio.spectype,2) eq 'LINER',nrliner)
rsf = where(strtrim(bsradio.spectype,2) eq 'SF',nrsf)
ragn = where(strtrim(bsradio.spectype,2) eq 'AGN',nragn)
rcomp = where(strtrim(bsradio.spectype,2) eq 'Composite',nrcomp)

print,'Blue spiral radio LINERs:    ',nrliner
print,'Blue spiral radio SF:        ',nrsf
print,'Blue spiral radio AGN:       ',nragn
print,'Blue spiral radio composite: ',nrcomp

print,''
end
