
;+
; NAME:
;       
;	COV_AMPLITUDES
;
; PURPOSE:
;
;	Compute the angular and spatial covariance amplitudes for a galaxy field to determine the cluster richness
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
;	Phi_star (normalization of the LF) must be determined in order to calculate B_gB. 
;
; REVISION HISTORY
;       Written by K. Willett                Sep 11
;-

; Differential form of the Schechter luminosity function, depending on normalization constant, characteristic magnitude, and slope

function schechter_mag, x

	common share1, phi_star_share, mr_star, alpha

	return, 0.4 * alog(10)  * phi_star_share * (10^(0.4 * (mr_star - x)))^(alpha + 1) * exp(-10^(0.4*(mr_star - x)))

end

; Begin program

pro cov_amplitudes, phi_star, barr_papers, barr_calc, stop=stop, noplot = noplot, ps=ps

	if n_elements(phi_star) eq 0 then phi_star = 0.5

	; Common block passed to differential LF (necessary since function is called in QSIMP, which takes no variables)

	common share1, phi_star_share, mr_star, alpha

; Create dummy field based off actual objects from Yee & Green (1987) - Paper III, Wurtz et al. (1997)

	; B_gB, total N_gal, background N_gal, radius of field in arcsec, redshift of blazar

	barr_yg87 = [473, 26, -90, 108, 731, 145, -210, 286, -165, 160, -46, 943, 993, 653, 792, 278, 304, 187, 773]
	nt_yg87 = [23.3, 9.2, 10.4, 15.6, 21.8, 8.1, 7.4, 12.0, 5.6, 13.7, 9.2, 25.1, 13.2, 13.8, 18.1, 16.6, 16.2, 11.4, 24.4]
	ng_yg87 = [10.2, 8.9, 13.6, 10.5, 11.7, 5.8, 10.1, 5.8, 10.4, 10.4, 11.6, 12.7, 6.6, 7.5, 8.9, 9.9, 8.4, 8.5, 11.5]
	fieldsize_yg87 = [51.2,50.7,56.8,54.0,56.5,47.8,52.4,54.2,56.8,50.7,54.5,58.3,46.6,51.9,51.6,54.0,47.7,49.2,53.0]
	z_yg87 = [0.450,0.624,0.457,0.367,0.614,0.462,0.611,0.344,0.422,0.420,0.422,0.634,0.652,0.594,0.599,0.361,0.436,0.555,0.590]

	barr_wsey97 = [114, 170, 392, -239, -53, -12, 559, -284, 120, 1262, -17, 694, -129, -283, 1326, 9, 1195, 166, -54, 383, -97, 31, -145, 845, -117, 252, 106, -52, 126, 302, 682, 22, 388, -120, 180, -159, 159, 501, 503, -73, 199, 190, 270, -21, 101]
	nt_wsey97 = [15.0, 15.0, 25.0, 2.0, 5.0, 3.0, 25.0, 9.0, 1.0, 12.1, 9.0, 34.3, 14.0, 3.0, 49.0, 13.0, 47.0, 3.0, 6.0, 13.3, 8.0, 9.0, 8.0, 38.0, 6.0, 14.4, 5.0, 9.0, 11.0, 19.0, 35.0, 8.0, 3.0, 2.0, 2.4, 7.0, 18.0, 6.0, 26.0, 10.0, 13.0, 8.0, 5.0, 0.0, 9.0]
	ng_wsey97 = [11.9, 10.4, 14.5, 8.3 , 6.4 , 3.2 , 10.0, 16.4, 0.1 , 0.2 , 9.4 , 16.2, 17.2, 10.6, 19.2, 12.8, 19.0, 0.7 , 7.4 , 4.5 , 10.3, 8.2 , 11.9, 15.8, 9.1 , 7.8 , 3.0 , 10.4, 7.7 , 10.8, 18.9, 7.4 , 0.1 , 5.1 , 0.1 , 11.3, 13.7, 0.1 , 12.3, 12.0, 7.8 , 3.2 , 0.5 , 0.1 , 6.4]
	fieldradius_wsey97 = [500, 500, 500, 500, 500, 416, 500, 500, 199, 229, 500, 500, 500, 500, 500, 500, 500, 313, 500, 457, 500, 500, 500, 500, 500, 500, 394, 500, 500, 500, 500, 500, 191, 500, 291, 500, 500, 273, 500, 500, 500, 500, 366, 132, 500]
	z_wsey97 = [0.339, 0.299, 0.444, 0.247, 0.190, 0.147, 0.287, 0.512, 0.061, 0.069, 0.267, 0.506, 0.548, 0.306, 0.638, 0.367, 0.615, 0.102, 0.218, 0.164, 0.297, 0.244, 0.344, 0.495, 0.260, 0.237, 0.152, 0.299, 0.235, 0.312, 0.605, 0.222, 0.034, 0.159, 0.055, 0.322, 0.407, 0.050, 0.352, 0.342, 0.237, 0.117, 0.068, 0.028, 0.190]
	fieldsize_wsey97 = zang(float(fieldradius_wsey97),z_wsey97)
	bgb_err_wsey97 = [178, 174, 220, 142, 126, 114, 210, 203, 128, 371, 150, 259, 215, 157, 392, 172, 331, 132, 135, 176, 155, 145, 165, 266, 147, 168, 139, 155, 153, 190, 293, 138, 239, 114, 125, 161, 192, 211, 216, 166, 162, 126, 137, 86, 141]

	barr_papers = [barr_yg87, barr_wsey97]
	nt_papers = [nt_yg87, nt_wsey97]
	ng_papers = [ng_yg87, ng_wsey97]
	fieldsize_papers = [fieldsize_yg87, fieldsize_wsey97]
	z_papers = [z_yg87, z_wsey97]

	ngal = n_elements(barr_papers)

	barr_calc = fltarr(ngal)
	
	; Use the Sebok LF with q_0 = 0.5 (S2 model)

		; Angular covariance amplitude

		gamma = 1.77				; Power-law index of galaxy-galaxy covariance relationship
		theta_papers = fieldsize_papers / 206265.
		A_gB = (nt_papers - ng_papers) / ng_papers * (3 - gamma) / 2. * theta_papers^(gamma-1)

		; Spatial covariance amplitude

		i_gamma = 3.87
		dl_papers = lumdist(z_papers, q0=0.50, H0=50, Omega_Lambda=0,/silent)	; Cosmological parameters from original paper
		da_papers = dl_papers / (1. + z_papers)^2					; Angular diameter distance (not used?)

		; Luminosity function

		;phi_star = 0.5		; Changing phistar affects the slope of the plot; ~0.5 best fits the data. NOT DETERMINED!!!
		phi_star_share = phi_star

		Mr_star = -21.22	; Changing Mstar does not significantly change the plot
		alpha = -1.2		; Changing alpha does not significantly change the plot

		; Lesson: relative normalization of the luminosity function is off; need to find out proper values to generate B_gb

		lowerlim = -25.
		psi = qsimp('schechter_mag', lowerlim, mr_star,/double)
		;print,''
		;print, 'Psi (M_lim = ',strtrim(lowerlim,2),') = ',psi
		;print,''

		B_gB = A_gB * ng_papers / i_gamma * dl_papers^(3. - gamma) / psi

		barr_calc = b_gb

	if ~keyword_set(noplot) then begin

	!p.multi=[0,1,2]	

	if keyword_set(ps) then begin
		ps_start, filename='~/Astronomy/Research/blazars/bgb_compare.ps', /quiet
		cs = 1
	endif else begin
		cs = 2
	endelse
	
	; Plot the Schechter LF

	magarr = fillarr(0.1,-25,-18)
	cgplot, magarr, schechter_mag(magarr), $
		charsize = cs, $
		/ylog, $
		xrange=[max(magarr),min(magarr)], $
		yrange=[1d-2,1d2], $
		xtitle='M!IR!N', $
		ytitle=greek('phi')+'[M]dM'

	cgplot, [mr_star,mr_star], 10^!y.crange, linestyle=2, color="Red", /overplot

	mr_star = -20 & cgplot, magarr, schechter_mag(magarr), linestyle=2, /overplot
	mr_star = -24 & cgplot, magarr, schechter_mag(magarr), linestyle=2, /overplot

	; Plot my calculated B_gB vs. the B_gB quoted in the literature

	yg87ind = indgen(n_elements(barr_yg87))
	wsey97ind = indgen(n_elements(barr_wsey97)) + n_elements(barr_yg87)

	cgplot, barr_papers[yg87ind], barr_calc[yg87ind], $
		psym=16, $
		charsize = cs, $
		xr = [-500,1500], $
		yr = [-500,3500], $
;		/iso, $
		xtitle='Published B!IgB!N', $
		ytitle='Calcuated B!IgB!N'

	cgplot, barr_papers[wsey97ind], barr_calc[wsey97ind], $
		psym=9, $
		/over

	al_legend, ['YG87','WSEY97'], $
		psym = [16,9], $
		/top, /left, $
		charsize=cs


	xr = !x.crange
	yr = !y.crange

	; Plot the one-to-one ratio line

	cgplot, indgen(5000)-1000, indgen(5000)-1000, $
		/overplot, $
		linestyle=2

	; Find the best fit for the phi star?

	if keyword_set(ps) then ps_end

	!p.multi=[0,1,1]	

	endif

	if keyword_set(stop) then stop

end
