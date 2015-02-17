
;+
; NAME:
;       
;	SCHECHTER_DIFF
;
; PURPOSE:
;
;	The differential Schechter luminosity function in units of absolute magnitude
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

function schechter_diff, x, p

	phi_star = p[0]
	m_star = p[1]
	alpha = p[2]
	m = x

	return, 0.4 * alog(10)  * phi_star * $
	(10^(0.4 * (m_star - m)))^(alpha + 1) * $
	exp(-10^(0.4*(m_star - m)))

end
