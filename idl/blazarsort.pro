
;+
; NAME:
;       
;	
;
; PURPOSE:
;
;	Sort the blazar csv output from the Python program
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
;       Written by K. Willett                Sep 11
;-

file_zs='~/Astronomy/Research/blazars/blazars_all_zs.csv'  
file_zns='~/Astronomy/Research/blazars/blazars_all_zns.csv'  

readcol,file_zs, $
	p_mg_zs,g_zs,p_acw_zs,i_zs,type_zs,p_el_zs,ra_zs,r_zs,p_dk_zs,objid_zs,nchild_zs,p_cw_zs,u_zs,p_edge_zs,z_zs,dec_zs,p_cs_zs,blazar_name_zs, $
	format=strjoin(replicate('a,',17))+'a',delimiter=',',skip=1 

readcol,file_zns, $
	p_mg_zns,g_zns,p_acw_zns,i_zns,type_zns,p_el_zns,ra_zns,r_zns,p_dk_zns,objid_zns,nchild_zns,p_cw_zns,u_zns,p_edge_zns,z_zns,dec_zns,p_cs_zns,blazar_name_zns, $
	format=strjoin(replicate('a,',17))+'a',delimiter=',',skip=1 

; Combine the zns and zs catalogs
; Remove all stars (type = 3)
; Sort by blazar_name, then by ra
; Create csv file

p_mg = [p_mg_zs, p_mg_zns]
p_acw = [p_acw_zs, p_acw_zns]
p_el = [p_el_zs, p_el_zns]
p_dk = [p_dk_zs, p_dk_zns]
p_cw = [p_cw_zs, p_cw_zns]
p_edge = [p_edge_zs, p_edge_zns]
p_cs = [p_cs_zs, p_cs_zns]
u = [u_zs, u_zns]
g = [g_zs, g_zns]
r = [r_zs, r_zns]
i = [i_zs, i_zns]
z = [z_zs, z_zns]
type = [type_zs, type_zns]
ra = [ra_zs, ra_zns]
dec = [dec_zs, dec_zns]
nchild = [nchild_zs, nchild_zns]
objid = [objid_zs, objid_zns]
blazar_name = [blazar_name_zs, blazar_name_zns]

zooind = multisort(blazar_name,ra)

objid = objid[zooind]
u = u[zooind]
g = g[zooind]
r = r[zooind]
i = i[zooind]
z = z[zooind]
type = type[zooind]
ra = ra[zooind]
dec = dec[zooind]
nchild = nchild[zooind]
p_mg = p_mg[zooind]
p_acw = p_acw[zooind]
p_el = p_el[zooind]
p_dk = p_dk[zooind]
p_cw = p_cw[zooind]
p_edge = p_edge[zooind]
p_cs = p_cs[zooind]
blazar_name = blazar_name[zooind]

galaxies_only = where(type eq 3)

objid = objid[galaxies_only]
u = u[galaxies_only]
g = g[galaxies_only]
r = r[galaxies_only]
i = i[galaxies_only]
z = z[galaxies_only]
type = type[galaxies_only]
ra = ra[galaxies_only]
dec = dec[galaxies_only]
nchild = nchild[galaxies_only]
p_mg = p_mg[galaxies_only]
p_acw = p_acw[galaxies_only]
p_el = p_el[galaxies_only]
p_dk = p_dk[galaxies_only]
p_cw = p_cw[galaxies_only]
p_edge = p_edge[galaxies_only]
p_cs = p_cs[galaxies_only]
blazar_name = blazar_name[galaxies_only]

zoo_array = [$
transpose(objid      ),$
transpose(u          ),$
transpose(g          ),$
transpose(r          ),$
transpose(i          ),$
transpose(z          ),$
transpose(type       ),$
transpose(ra         ),$
transpose(dec        ),$
transpose(nchild     ),$
transpose(p_mg       ),$
transpose(p_acw      ),$
transpose(p_el       ),$
transpose(p_dk       ),$
transpose(p_cw       ),$
transpose(p_edge     ),$
transpose(p_cs       ),$
transpose(blazar_name)]


writefile_zs = '~/Astronomy/Research/blazars/blazars_all_sorted_zs.csv'  
openw, lun1, writefile_zs, /get_lun
printf, lun1, zoo_array
close, lun1
free_lun, lun1

; Needed to write short VIM script to combine all lines at the end; worked by searching on beginning string 12376.

end
