
; Differential form of the Schechter luminosity function in absolute magnitudes

function schechter_mag, x

	common schechter_block, phi_star, mr_star, alpha

	return, 0.4 * alog(10)  * phi_star* (10^(0.4 * (mr_star - x)))^(alpha + 1) * exp(-10^(0.4*(mr_star - x)))

end

pro psiz

	common schechter_block, phi_star, mr_star, alpha

	zarr = fillarr(0.01,0.1,0.7)
	nz = n_elements(zarr)
	nt = 100.
	nb = 20.

	fieldradius = 500	; kpc
	fieldsize = zang(float(fieldradius),zarr,/silent,/wmap7)
	Mr_star = -21.1837		; Blanton et al. (2003)
	sdss_lim = 22.2			; SDSS apparent limiting magnitude
	sdss_lim_abs = sdss_lim - 5*alog10(lumdist(zarr,/silent,/wmap7)*1e6) + 5
	abs_mstar = Mr_star + 2
	countingmag = sdss_lim_abs < abs_mstar

	psiarr = fltarr(nz)

	;bgb, nt, nb, fieldsize, z, countingmag, bgb_result, bgb_result_err, oldcosmo=oldcosmo, stop=stop

	; Common block passed to differential LF 
	; (necessary since function is called in QSIMP, which takes no variables)

	for i=0, nz - 1 do begin

		; Angular covariance amplitude (Yee & Green 1987)

		gamma = 1.77			; Power-law index of galaxy-galaxy covariance relationship
		theta = fieldsize[i] / 206265.
		A_gB = (nt - nb) / nb * (3 - gamma) / 2. * theta^(gamma-1)

		; Spatial covariance amplitude

		i_gamma = 3.87
		if keyword_set(oldcosmo) then dl = lumdist(zarr[i], q0=0.50, H0=50, Omega_Lambda=0,/silent) else $
			dl = lumdist(zarr[i],/silent,/wmap7)
		da = dl / (1. + zarr[i])^2					; Angular diameter distance

		; Luminosity function parameters (Blanton et al. 2003)

		h = 0.71
		phi_star = 1.49e-2 * h^3	; Normalization constant 
		Mr_star = -20.44+5*alog10(h)	; Characteristic absolute magnitude of LF
		alpha = -1.05			; Slope of the power-law portion of LF

		lowerlim = -28.			; More negative numbers cause underflow errors.

		psi = qsimp('schechter_mag', lowerlim, countingmag[i],/double)

		psiarr[i] = psi

		ng = 1e8		; Taken from Best paper - surface density of all galaxies per sr

		n_net = (nt - nb)
		a_theta = !pi * theta^2

		; Muzzin et al. (2007)
		B_gB = n_net * (3. - gamma) * da^(gamma - 3.) * theta^(gamma - 1.) / (2. * a_theta * i_gamma * psi)

		; Wold et al, Best et al, Zauderer et al. 
		;B_gB = A_gB * ng / i_gamma * da^(gamma - 3.) / psi	; YG87 have the exponent wrong.

		; Compute uncertainty on bgb (Yee & Lopez-Cruz 1999)

		err_bgb = b_gb * sqrt(abs(nt - nb) + 1.3^2 * nb) / abs(nt - nb)

		bgb_result = B_gB
		bgb_result_err = err_bgb

	endfor

	cgplot, zarr, psiarr, $
		yr=[0,0.01], $
		xtitle='Redshift', $
		ytitle='Psi'

	cgplot, zarr, countingmag, $
		xtitle='Redshift', $
		ytitle='Counting magnitude'

	if keyword_set(stop) then stop

end
