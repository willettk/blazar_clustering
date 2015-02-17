
;+
; NAME:
;       
;	MASTERS_REDSPIRALS
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
pro masters_redspirals,ps=ps

file='~/Astronomy/proposals/redspirals_chandra/masters2010/mnr_16503_sm_TableA1.txt'

readcol, file, $
	objid, ra, dec, redshift, $
	absmag_r, absmag_r_err, color_gr, color_gr_err, redness, ellipticity, fracdeV, $
	format='a,f,f,f,f,f,f,f,f,f', $
	comment='%', $
	/silent, $
	skip=4
	
cgplot, color_gr, absmag_r, $
	psym=16, $
	xtitle='(g-r)', $
	ytitle='M!Ir!N', $
	charsize=1.5

radiofile='~/Astronomy/proposals/redspirals_chandra/redspirals_radio.txt'

;readcol, radiofile, $
;	rowno, sdssid, $
;	ra_sdss, dec_sdss, $
;	angsep, $
;	nedname, ned_ra, ned_dec, $
;	objtype, $
;	z_ned, z_sdss, $
;	flux_nvss, err_nvss, flux_first, $
;	ew_nii, ew_halpha, ew_oiii, ew_hbeta, $
;	format='i,a,f,f,f,a,a,a,a,f,f,f,f,f,f,f,f,f', $
;	delimiter='|', $
;	/silent, $
;	skip=7, numline=17
;

sdss_lines_file = '~/Astronomy/proposals/redspirals_vla/redspirals_lines_willettk.fit'

rs = mrdfits(sdss_lines_file,1,/silent)

;sdssid=strtrim(sdssid,2)
;match, dr6id, sdssid, drind, sdssind, count=mcount
;bptind = drind[uniq(dr6id[drind])]	; Proud of this nesting.
;neduind = sdssind[uniq(sdssid[sdssind])]

nsig = 3.
nsig_neiii = 1.
goodind = where($
;	(neiii_3869_flux gt nsig_neiii * neiii_3869_flux_err) and $
	;((oii_3726_flux + oii_3729_flux) gt nsig * sqrt(oii_3726_flux_err^2 + oii_3729_flux_err^2)) and $
	(rs.nii_6584_flux gt nsig * rs.nii_6584_flux_err) and $
	(rs.oiii_5007_flux gt nsig * rs.oiii_5007_flux_err) and $
	(rs.h_alpha_flux gt nsig * rs.h_alpha_flux_err) and $
	(rs.h_beta_flux gt nsig * rs.h_beta_flux_err), $
	goodcount)

rsgood = rs[goodind]

nha = alog10(rsgood.nii_6584_flux/rsgood.h_alpha_flux)
ohb = alog10(rsgood.oiii_5007_flux/rsgood.h_beta_flux)

dir = '~/Astronomy/proposals/redspirals_chandra/'
if keyword_set(ps) then begin
	ps_start, filename=dir+'masters_redspirals.eps',/color,/landscape, xs=10, ys=5, /encap
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

neo = alog10(rsgood.neiii_3869_flux / (rsgood.oii_3726_flux + rsgood.oii_3729_flux))
gz =  (rsgood.modelmag_g) - (rsgood.modelmag_z)

nsig_ne3 = 3.
badne = where(rsgood.neiii_3869_flux lt nsig_ne3 * rsgood.neiii_3869_flux_err,bnc)
goodne = where(rsgood.neiii_3869_flux gt nsig_ne3 * rsgood.neiii_3869_flux_err,gnc)
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

if bnc gt 0 then cgarrow, neo[badne], gz[badne], neo[badne] - 0.1, gz[badne], /data, color=' grey', hsize=hs

; Computing the Chandra observation time

;totaltime = 0.
;zg = z[goodind]
;ztime = zg[[liner]]
;sortzind = sort(ztime)
;zind = 0
;fac2 = 0.00090765913
;
;	; Find the closest LINER with NeIII emission 
;
;	match, ztime, z, a, b
;	neon = where(neiii_3869_flux gt nsig_ne3 * neiii_3869_flux_err)
;	neonliner = setintersection(b,neon)
;	minind = where(z eq min(z[neonliner]))
;	print,'Min ind = ',minind
;	print,'zmin = ',z[minind]
;	print,'Objid = ',objid[minind]
;
;print,''
;
;objidarr = strarr(14)
;zarr = fltarr(14)
;
;l0 = 1d40
;mpc2cm=3.086d24
;
;ctrate = 1e-4 * [7.487,6.979,6.964,6.158,5.547,5.669,5.842]
;print,transpose(string(25e-3 / ctrate,format='(f5.1)'))
;
;while totaltime lt 300 do begin
;
;	;flux = l0 / (4d * !dpi * (dl*mpc2cm)^2)
;
;	newz = ztime[sortzind[zind]]
;	newtime = fac2 * (lumdist(newz,/wmap7,/silent))^2
;	;if (where(liner eq sortzind[zind]))[0] ne -1 then tag = 'LINER'
;	;if (where(agn eq sortzind[zind]))[0] ne -1 then tag = 'AGN'
;	;if (where(comp eq sortzind[zind]))[0] ne -1 then tag = 'Comp'
;	i = where(z eq newz)
;;	print,string(objid[i])
;	print,string(newz,format='(f5.3)'), '  ', string(newtime,format='(f4.1)'), '  ', string(petroR50_r[i],format='(f4.1)')
;	objidarr[zind] = objid[i]
;	zarr[zind] = newz
;	totaltime += newtime
;	zind += 1
;endwhile
;
;match, zarr, z, a1, b1
;neon_mysamp = setintersection(b1,neon)
;print,''
;print,'LINERs with NeIII: ',objid[neon_mysamp]
;
;newinds = sortzind[0:zind-1]
;
;ctime = totaltime - newtime
;ngal = zind
;
;print,''
;print,'Total exposure time [ksec]: ',ctime
;print,'Total number of galaxies: ',ngal
;
;print,''

if keyword_set(ps) then ps_end

; Read in Karen's bar data

;klmfile = '~/Astronomy/proposals/redspirals_chandra/redspiralbars_klmasters.fit'
;result = mrdfits(klmfile,1,header,/silent)
;
;match, result.objid, objidarr, aa, bb
;i = result[aa].petromag_mi
;g = result[aa].petromag_mg
;mstar = (4.52-i)/2.5 + (-0.963+1.032*(g-i))
;print,[transpose(string(mstar,format='(f7.2)')),transpose(string(result[aa].pbar,format='(f7.2)'))]


; Print list of ra, dec, spectral type

openw, lun, '~/Astronomy/proposals/redspirals_vla/rs_bytype.txt',/get_lun,width=120
printf, lun, 'ra, dec, redshift, objid, spectype,'
printf, lun, [transpose(string(rsgood[liner].ra)), transpose(string(rsgood[liner].dec)),transpose(string(rsgood[liner].redshift)),transpose(string(rsgood[liner].dr7objid)),transpose(replicate('LINER',n_elements(liner)))]+','
printf, lun, [transpose(string(rsgood[sf].ra)), transpose(string(rsgood[sf].dec)),transpose(string(rsgood[sf].redshift)),transpose(string(rsgood[sf].dr7objid)),transpose(replicate('SF',n_elements(sf)))]+','
printf, lun, [transpose(string(rsgood[agn].ra)), transpose(string(rsgood[agn].dec)),transpose(string(rsgood[agn].redshift)),transpose(string(rsgood[agn].dr7objid)),transpose(replicate('AGN',n_elements(agn)))]+','
printf, lun, [transpose(string(rsgood[comp].ra)), transpose(string(rsgood[comp].dec)),transpose(string(rsgood[comp].redshift)),transpose(string(rsgood[comp].dr7objid)),transpose(replicate('Composite',n_elements(comp)))]+','
close, lun

; Match the LINERs to the RA and dec from Karen's table

;dr7good = dr7objid[goodind]
;match, objid, dr7good[liner], a, b
;
;tr = ra[a]
;td = dec[a]
;sortdec = reverse(sort(td))
;raliner = tr[sortdec]
;decliner = td[sortdec]
;
;print,''
;for i=0,n_elements(a) - 1 do begin
;
;	ratemp = sixty(raliner[i]/15.)
;	rastring = [string(ratemp[0:1],format='(i02)'),string(ratemp[2],format='(f04.1)')]
;	dectemp = sixty(decliner[i])
;	if dectemp[0] eq 0 and dectemp[1] lt 0 then decsign = '-' else decsign = '+'
;	if dectemp[0] gt 0 then decsign = '+'
;	if dectemp[0] lt 0 then decsign = '-'
;	decstring = [decsign+string(abs(dectemp[0]),format='(i02)'),string(abs(dectemp[1]),format='(i02)'),string(dectemp[2],format='(f04.1)')]
;	;print,[rastring,decstring]
;
;endfor

stop

;;;;;;;;;;;;; Red elliptical galaxies ;;;;;;;;;;;;;;;;;;;;;

redelliptical_file='~/Astronomy/proposals/redspirals_chandra/re5_willettk.csv'

readcol, redelliptical_file, $
	specobjid,ra,dec,redshift,p_el_debiased,$
	nii_6584_flux,nii_6584_flux_err,oiii_5007_flux,oiii_5007_flux_err,neiii_3869_flux,neiii_3869_flux_err,oii_3726_flux,oii_3726_flux_err,oii_3729_flux,oii_3729_flux_err,h_beta_flux,h_beta_flux_err,h_alpha_flux,h_alpha_flux_err,g,z,absmagG,absmagZ, $
	format='a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f', $
	comment='%', $
	/silent, $
	skip=1

!p.multi=[0,1,2]
	
cgplot, alog10(nii_6584_flux/h_alpha_flux), alog10(oiii_5007_flux/h_beta_flux), $
	xr=[-1.5,0.5], /xstyle, $
	yr=[-1.0,1.5], /ystyle, $
	title='Red elliptical galaxies', $
	xtitle='log ([NII] '+greek('lambda')+'6584/H'+greek('alpha')+')', $
	ytitle='log ([OII] '+greek('lambda')+'5007/H'+greek('beta')+')', $
	charsize=2, $
	psym=9

xkauffmann=fillarr(0.01,-10,0)
xkewley=fillarr(0.01,-10,0.4)
kauffmann = 0.61 /(xkauffmann - 0.05)+1.3
kewley = 0.61 /(xkewley - 0.47)+1.19

cgplot, xkauffmann, kauffmann, $
	/over, $
	linestyle=2

cgplot, xkewley, kewley, $
	/over, $
	linestyle=1

cgplots, [-0.2,0.5], [0.61 /(-0.2 - 0.47)+1.19,0.98], linestyle=0

cgtext, /data, -0.7, 1.0, 'Seyferts', charsize=2.0
cgtext, /data, 0.1, -0.1, 'LINERs', charsize=2.0
cgtext, /data, -1.3, 0.0, 'Star-forming', charsize=2.0
cgtext, /data, -0.2, -0.7, 'Composite', charsize=1.5

ind = where(neiii_3869_flux gt 0. and (oii_3726_flux + oii_3729_flux) gt 0.)

color = absmagG[ind] - absmagZ[ind]
color = G[ind] - z[ind]
tbt_lr = alog10(neiii_3869_flux[ind]/(oii_3726_flux[ind] + oii_3729_flux[ind]))

cgplot, tbt_lr, color, $
	yr=[-3.0,0.5], /xstyle, $
;	xr=[-0.5,2.5], /ystyle, $
	title='Red elliptical galaxies', $
	xtitle='log ([NeIII]/[OII])', $
	ytitle='!E0!N(g-z)', $
	charsize=2, $
	psym=16

x = fillarr(0.1,-10,10)

cgplot, /over, x, (-1.2)*x -0.4, $
	thick=3, linestyle=2



stop

end
