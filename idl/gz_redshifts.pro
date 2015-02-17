
;+
; NAME:
;       
;	GZ_REDSHIFTS
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

dir='~/Astronomy/Research/blazars/'
file_zns_el = dir+'zoonospec_elliptical_willettk.csv'
file_zns_sp = dir+'zoonospec_spiral_willettk.csv'
file_zs_el  = dir+'zoospec_elliptical_willettk.csv'
file_zs_sp  = dir+'zoospec_spiral_willettk.csv'

readcol, file_zns_el, $
	objID_zns_el,p_el_zns,z_zns_el,zerr_zns_el, $
	format='a,f,f,f', $
	/skip

readcol, file_zns_sp, $
	objID_zns_sp,p_sp_zns,z_zns_sp,zerr_zns_sp, $
	format='a,f,f,f', $
	/skip

readcol, file_zs_el, $
	objID_zs_el,p_el_zs,z_zs_el,zerr_zs_el, $
	format='a,f,f,f', $
	/skip

readcol, file_zs_sp, $
	objID_zs_sp,p_sp_zs,z_zs_sp,zerr_zs_sp, $
	format='a,f,f,f', $
	/skip

file='~/Astronomy/Research/blazars/allblazars_remdup.cat'
readcol, file, $
	bname, bra, bdec, bz, btype, $
	format='a,f,f,f,a' , $
	/silent, $
	/skip

!p.multi=[0,2,1]
binsize = 0.02
cs = 1.5
cghistoplot,z_zns_el, xr=[0,1], bins=binsize, datacolor='red', title='ZooNoSpec', charsize=cs, /ylog, min_val=1
cghistoplot,z_zns_sp, /oplot, bins=binsize, datacolor='blue'
cghistoplot,bz, /oplot, bins=binsize, datacolor='green'

cghistoplot,z_zs_sp, xr=[0,1], bins=binsize, datacolor='blue', title='ZooSpec', charsize=cs, /ylog, min_val=1
cghistoplot,z_zs_el, /oplot, bins=binsize, datacolor='red'
cghistoplot,bz, /oplot, bins=binsize, datacolor='green'

h_zns_el = histogram(z_zns_el,bin=binsize,omin=omin_zns_el,max=1.0)
h_zns_sp = histogram(z_zns_sp,bin=binsize,omin=omin_zns_sp,max=1.0)
h_zs_el = histogram(z_zs_el,bin=binsize,omin=omin_zs_el,max=1.0)
h_zs_sp = histogram(z_zs_sp,bin=binsize,omin=omin_zs_sp,max=1.0)
h_bz = histogram(bz,bin=binsize,omin=omin_bz,max=1.0)

zbin = fillarr(binsize,omin_zns_el,1.0)

print,'ZNS - el',zbin[closeto(h_zns_el - h_bz,0.)], n_elements(where(bz le zbin[closeto(h_zns_el - h_bz,0.)]))
print,'ZNS - sp',zbin[closeto(h_zns_sp - h_bz,0.)], n_elements(where(bz le zbin[closeto(h_zns_sp - h_bz,0.)]))
print,'ZS - el',zbin[closeto(h_zs_el - h_bz,0.)], n_elements(where(bz le zbin[closeto(h_zs_el - h_bz,0.)]))
print,'ZS - sp',zbin[closeto(h_zs_sp - h_bz,0.)], n_elements(where(bz le zbin[closeto(h_zs_sp - h_bz,0.)]))

stop

end
