
;+
; NAME:
;       
;	BGB
;
; PURPOSE:
;
;	Compute the spatial covariance amplitude for a collection of objects on the sky
;
; INPUTS:
;
;	NT			- number of total galaxies within measured field
;
;	NB			- number of background galaxies for control field
;
;	FIELDSIZE		- radius of field [arcsec]
;
;	Z			- redshift of central galaxy
;
;	COUNTINGMAG		- absolute magnitude down to which the LF is computed
;
;
; OUTPUTS:
;
;	BGB_RESULT 		- spatial covariance amplitude
;
;	BGB_RESULT_ERROR 	- uncertainty in spatial covariance amplitude
;
; KEYWORDS:
;
;	COSMO 	- 		Set to compute luminosity distances according to older cosmologies:
;					'smith' - H0 = 75, q0 = 0.0	(Smith et al. 1995)
;					'wurtz' - H0 = 50, q0 = 0.02	(Wurtz et al. 1997)
;					'yee' - H0 = 50, q0 = 0.50	(Yee & Green 1987)
;				Default is to use the WMAP7 concordance cosmology. 
;
; EXAMPLE:
;
;	IDL> nt = [100,130,145]
;	IDL> nbg = [50, 49, 40]
;	IDL> z = [0.1, 0.2, 0.1]
;	IDL> fieldsize = replicate(0.1,3)
;	IDL> countingmag = [-19.,-21.,-20.]
;	IDL> bgb, nt, nbg, fieldsize, z, countingmag, bgb_result, bgb_result_err
;
; NOTES:
;
;	
;
; REVISION HISTORY
;       Written by K. Willett                Dec 11
;	Background counts fixed to control fields	Mar 12
;-

; Differential form of the Schechter luminosity function in absolute magnitudes

function schechter_mag, x

	common schechter_block, phi_star, mr_star, alpha

	return, 0.4 * alog(10)  * phi_star* (10^(0.4 * (mr_star - x)))^(alpha + 1) * exp(-10^(0.4*(mr_star - x)))

end

pro bgb, nt, nb, fieldsize, z, countingmag, bgb_result, bgb_result_err, cosmo=cosmo, stop=stop, lf = lf

	; Common block passed to differential LF 
	; (necessary since function is called in QSIMP, which takes no variables)

	common schechter_block, phi_star, mr_star, alpha

	nt = float(nt)
	nb = float(nb)

	; Angular covariance amplitude (Yee & Green 1987)

	gamma = 1.77			; Power-law index of galaxy-galaxy covariance relationship
	theta = fieldsize / 206265.
	A_gB = (nt - nb) / nb * (3. - gamma) / 2. * theta^(gamma-1)

	; Spatial covariance amplitude

	i_gamma = 3.87

	if n_elements(cosmo) eq 0 then begin
		dl = lumdist(z,/silent,/wmap7)
		;h=0.71
		h=0.704     ; Incorrectly rounded before; see Hinshaw+09
	endif else begin
		if cosmo eq 'wurtz' then begin
			dl = lumdist(z, q0=0.02, H0=50, Omega_Lambda=0,/silent)
			h = 0.50
		endif else if cosmo eq 'smith' then begin
			dl = lumdist(z, q0=0.0, H0=75, Omega_Lambda=0,/silent)
			h = 0.75
		endif else if cosmo eq 'yee' then begin
			dl = lumdist(z, q0=0.50, H0=50, Omega_Lambda=0,/silent)
			h = 0.50
		endif
	endelse

	da = dl / (1. + z)^2					; Angular diameter distance

	; Luminosity function parameters

	if n_elements(lf) eq 0 then lf = 'dr6ages'
	nz = n_elements(z)

	case lf of
		'dr1': begin
				phi_star = 1.49e-2 * h^3	; Blanton et al. (2003)
				Mr_star = -20.44+5*alog10(h)	; SDSS, pre-DR1
				alpha = -1.05			; 
			     end
		'dr6': begin
				phi_star = 0.90e-2 * h^3	; Montero-Dorta & Prada (2006)
				Mr_star = -20.73+5*alog10(h)	; SDSS DR6
				alpha = -1.23			; 
			     end
		'6dfgs': begin
				phi_star = 10.^(-2.081)	* h^3	; Jones et al. (2006)
				Mr_star = -20.98+5*alog10(h)	; 6dFGS
				alpha = -1.21			; 
			     end
		'ages': begin
				phi_star = 0.0159 * h^3		; Cool et al. (2012)
				Mr_star = -20.58 + 5*alog10(h)	; AGN and Galaxy Evolution survey
				alpha = -1.05			; 
				
				; Evolving LF from Table 4

				zarr = [0.15,0.25,0.35,0.45,0.55,0.65]
				phi_star_arr = [0.0159,0.0152,0.0124,0.0144,0.0108,0.0105] * h^3
				Mr_star_arr = [-20.58,-20.81,-20.81,-20.99,-21.29,-21.38] + 5*alog10(h)
				ind = where(z lt zarr,zc)
				if zc gt 0 then zind = ind[0] else zind = -1

				phi_star_t = fltarr(nz)
				Mr_star_t = fltarr(nz)
				alpha_t = fltarr(nz)

				for i=0,nz - 1 do begin
					ind = where(z[i] lt zarr,zc)
					if zc gt 0 then zind = ind[0] else zind = -1

					phi_star_t[i] = phi_star_arr[zind]
					Mr_star_t[i] = Mr_star_arr[zind]
					alpha_t[i] = alpha
				endfor
			     end
		'dr6ages': begin        ; Use DR6 for low-z, AGES for increasing redshift
				phi_star = 0.90e-2 * h^3	; DR6
				Mr_star = -20.73+5*alog10(h)
				alpha_dr6 = -1.23			; 

				alpha_ages = -1.05			; 
				
				; Evolving LF from Table 4

				zarr = [0.15,0.25,0.35,0.45,0.55,0.65]
				phi_star_arr = [phi_star/(h^3),0.0152,0.0124,0.0144,0.0108,0.0105] * h^3
				Mr_star_arr = [mr_star - 5*alog10(h),-20.81,-20.81,-20.99,-21.29,-21.38] + 5*alog10(h)

				phi_star_t = fltarr(nz)
				Mr_star_t = fltarr(nz)
				alpha_t = fltarr(nz)

				for i=0,nz - 1 do begin
					ind = where(z[i] lt zarr,zc)
					if zc gt 0 then zind = ind[0] else zind = -1

					phi_star_t[i] = phi_star_arr[zind]
					Mr_star_t[i] = Mr_star_arr[zind]
					if zind eq 0 then alpha_t[i] = alpha_dr6 else alpha_t[i] = alpha_ages
				endfor
			     end
		'gama': begin
				phi_star = 0.90 * h^3/100. 	; Loveday et al. (2012)
				Mr_star = -20.73+5*alog10(h)	; GAMA
				alpha = -1.26			; 

				p = 1.6				; step-wise max. likelihood values from Table 5
				q = 0.2

				z0 = 0.1			; normalization redshift of the LF

				phi_star_t = phi_star * 10.^(0.4 * p * (z - z0))	; Evolution according to parameterization of Lin et al. (1999)
				Mr_star_t = mr_star - q*(z - z0)
				alpha_t = fltarr(nz) + alpha
			     end
    'ramosalmeida': begin

				alpha = -1.3
				
				; Evolving LF from Table 4

				zarr = (indgen(5)+1)/5.
				phi_star_arr = [0.0038,0.0037,0.0035,0.0033,0.0031]
				Mr_star_arr_r = [-21.43,-22.08,-22.77,-22.62,-22.87]
				Mr_star_arr_i = [-21.76,-22.47,-23.27,-23.59,-23.84]

				phi_star_t = fltarr(nz)
				Mr_star_t = fltarr(nz)
				alpha_t = fltarr(nz)

				for i=0,nz - 1 do begin
                    if z[i] > 0.4 then Mr_star_arr = Mr_star_arr_r else Mr_star_arr = Mr_star_arr_i
					ind = where(z[i] lt zarr,zc)
					if zc gt 0 then zind = ind[0] else zind = -1

					phi_star_t[i] = phi_star_arr[zind]
					Mr_star_t[i] = Mr_star_arr[zind]
					alpha_t[i] = alpha
				endfor
			     end
	endcase

;	stop

	lowerlim = -28.			; More negative numbers cause underflow errors.
	if lf eq 'ramosalmeida' then lowerlim = countingmag - 1

	; Array + loop. Necessary because variables passed to QSIMP now change, which have to be done as a common block, but QSIMP only expects a single result from the function.  
	n = n_elements(phi_star_t)
	if n gt 1 then begin
		psi = dblarr(n)
		for i=0,n - 1 do begin
			phi_star = phi_star_t[i]
			Mr_star = Mr_star_t[i]
			alpha = alpha_t[i]
			psi[i] = qsimp('schechter_mag', lowerlim, countingmag[i],/double)
		endfor
	endif else psi = qsimp('schechter_mag', lowerlim, countingmag,/double)

	n_net = (nt - nb)
	a_theta = !pi * theta^2

	; Muzzin et al. (2007)
	B_gB = n_net * (3. - gamma) * da^(gamma - 3.) * theta^(gamma - 1.) / (2. * a_theta * i_gamma * psi)

	; Wold et al, Best et al, Zauderer et al. method
	;ng = 1e8		; Taken from Best paper - surface density of all galaxies per sr
	;B_gB = A_gB * ng / i_gamma * da^(gamma - 3.) / psi	; YG87 have the exponent wrong.

	; Compute uncertainty on bgb (Yee & Lopez-Cruz 1999)

	err_bgb = b_gb * sqrt(abs(nt - nb) + 1.3^2 * nb) / abs(nt - nb)

	bgb_result = B_gB
	bgb_result_err = err_bgb

    ; Print out relevant variables
    
    print,''
    print, format='(A,I21)'   ,'n_net: ', n_net
    print, format='(A,F21.2)','gamma: ', gamma
    print, format='(A,F21.2)','da: ', da
    print, format='(A,E21.2)','theta: ', theta
    print, format='(A,E21.2)','a_theta: ', a_theta
    print, format='(A,F21.2)','i_gamma: ', i_gamma
    print, format='(A,F21.2)','counting mag: ', countingmag
    print,''
    print, format='(A,E21.3)','phi_star: ', phi_star
    print, format='(A,F21.3)','Mr_star: ', Mr_star
    print, format='(A,F21.3)','alpha: ', alpha
    print,''
    print, format='(A,F21.4)','h: ', h
    print,''
    print, format='(A,E21.2)','psi: ', psi
    print,''

	if keyword_set(stop) then stop

end
