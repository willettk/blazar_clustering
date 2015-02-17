
;+
; NAME:
;       
;	PLOTKIN
;
; PURPOSE:
;
;	Read in and plot data on BL Lacs from Plotkin et al. (2010)
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
;       Written by K. Willett                Jul 11
;-

plotkinfile5 = '~/Astronomy/Research/blazars/plotkin_bllacs_t5.txt'
plotkinfile7 = '~/Astronomy/Research/blazars/plotkin_bllacs_t7.txt'

readcol, plotkinfile5, $
	Name, RAdeg, DEdeg, Class, InNED, umag, gmag, rmag, imag, zmag, $
	;Zspec, f_Zspec, CaII_break, SpIndex, Comments, $
	skipline=49, format='A, F, F, A, A, F, F, F, F, F, F, A, F, F, A,', /quiet


readcol, plotkinfile7, $
	Name7  , FAGN  , logOL , RFlag , Rflux , logRL , XFlag , Xflux , logXL , Ind_RO, Ind_OX, z, $
	skipline=67, format='A, F, F, A, F, F, A, F, F, F, F, f', /quiet

!p.multi=[0,2,2]
cs = 1.2
ss=0.5

cgplot, umag-gmag, gmag-rmag, psym=5, symsize=ss, charsize=cs, xtitle='u-g', ytitle='g-r', yr=[-1,2.5], color="Red"
hline, 1.4, /data, color=cgcolor("Black"), linestyle=2
cgtext, 3.0, 2.5, 'LRG', /data, charsize=1.5
cgplot, gmag-rmag, rmag-imag, psym=5, symsize=ss, charsize=cs, xtitle='g-r', ytitle='r-i', xr=[-1,3], yr=[-0.5,1], color="Red"
vline, 1.4, /data, color=cgcolor("Black"), linestyle=2
cgplots, [0.35,0.35], [0.13,!y.crange[0]], linestyle=2
cgplots, [!x.crange[0],0.35], [0.13,0.13], linestyle=2
cgtext, -0.9, -0.4, 'WD', /data, charsize=1.5
cgtext, 2.0, -0.1, 'LRG', /data, charsize=1.5
cgplot, rmag-imag, imag-zmag, psym=5, symsize=ss, charsize=cs, xtitle='r-i', ytitle='i-z', color="Red"
cgplot, gmag-rmag, imag,      psym=5, symsize=ss, charsize=cs, xtitle='g-r', ytitle='i', yrange=[21,15], xrange=[-1,2.5], color="Red"

zoospec = '~/Astronomy/Research/blazars/plotkin_dr8_zoospec.csv'
zoonospec = '~/Astronomy/Research/blazars/plotkin_dr8_zoonospec.csv'

readcol, zoospec, $
	zs_col0,zs_objID,zs_ra,zs_dec,zs_run,zs_rerun,zs_camcol,zs_field,zs_type,$
	zs_modelMag_u,zs_modelMag_g,zs_modelMag_r,zs_modelMag_i,zs_modelMag_z,$
	zs_p_el,zs_p_cw,zs_p_acw,zs_p_edge,zs_p_dk,zs_p_mg,zs_p_cs, $
	skipline=1, format='a,a,f,f,i,i,i,i,a,f,f,f,f,f,f,f,f,f,f,f,f', /quiet

readcol, zoonospec, $
	zns_col0,zns_objID,zns_ra,zns_dec,zns_run,zns_rerun,zns_camcol,zns_field,zns_type,$
	zns_modelMag_u,zns_modelMag_g,zns_modelMag_r,zns_modelMag_i,zns_modelMag_z,$
	zns_p_el,zns_p_cw,zns_p_acw,zns_p_edge,zns_p_dk,zns_p_mg,zns_p_cs, $
	skipline=1, format='a,a,f,f,i,i,i,i,a,f,f,f,f,f,f,f,f,f,f,f,f', /quiet

!p.multi=[0,1,1]

ps_start, file='~/Astronomy/Research/blazars/plotkin_zoo.eps', /encap, /color, /quiet

cgplot, zs_modelmag_g - zs_modelmag_r, zs_p_el, psym=17, color="red", xtitle='g-r',ytitle='p!Iel!N', thick=4, xthick=4, ythick=4, charthick=2, charsize=2
cgplot, /over, zns_modelmag_g - zns_modelmag_r, zns_p_el, psym=18, color="blue"

ps_end

;;;;;;;;;;
; Plot the environment for a 3 arcminute radius around each Plotkin BL Lac
;;;;;;;;;;

zs_3arcmin = '~/Astronomy/Research/blazars/plotkin_zs_3arcmin.csv'
zns_3arcmin = '~/Astronomy/Research/blazars/plotkin_zns_3arcmin.csv'

readcol, zs_3arcmin, $
	zs3_col0,zs3_objID,zs3_ra,zs3_dec,zs3_type,zs3_modelMag_u,zs3_modelMag_g,zs3_modelMag_r,zs3_modelMag_i,zs3_modelMag_z,zs3_p_el, $
	skipline=1, format='a,a,f,f,a,f,f,f,f,f,f,f', /quiet

readcol, zns_3arcmin, $
	zns3_col0,zns3_objID,zns3_ra,zns3_dec,zns3_type,zns3_modelMag_u,zns3_modelMag_g,zns3_modelMag_r,zns3_modelMag_i,zns3_modelMag_z,zns3_p_el, $
	skipline=1, format='a,a,f,f,a,f,f,f,f,f,f,f', /quiet

; Find unique neighbors to the Plotkin BL Lacs and histogram it

zs3_nos = strmid(zs3_col0,13,99)
zns3_nos = strmid(zns3_col0,13,99)

both_zs3nos = [zs3_nos, zns3_nos]

; Find the object IDs for both samples

both_zs3objids = [zs3_objid, zns3_objid]
both_zsobjids = [zs_objid, zns_objid]

; Match
match, both_zs3objids, both_zsobjids, aa, bb

; Find the indices for objects that are not the parent BL Lac target
norepeats = setdifference(indgen(n_elements(both_zs3objids)),aa)

norepeat_nos = bothnos[norepeats]

h=histogram(fix(bothnos))
h_norepeats = histogram(fix(bothnos[norepeats]))

zz = where(h eq 0, zc)

cghistoplot, h_norepeats, xtitle="Number of BL Lac neighbors w/in 3'", ytitle='N'
;cghistoplot, h[zz], /oplot, /line_fill, orientation=45, spacing=1d-1, polycolor='black'

end
