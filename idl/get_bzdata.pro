
;+
; NAME:
;       
;	GET_BZDATA
;
; PURPOSE:
;
;	Restore SDSS data on blazar neighbors from structure saved in file
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
;       Written by K. Willett                Dec 11
;-

function get_bzdata, bz

restore,'/Users/willettk/Astronomy/Research/blazars/ab_structure.sav'

return, bz

end
