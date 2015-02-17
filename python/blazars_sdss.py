"""
Module to retrive and plot SDSS DR8 data for galaxies near blazars in the blazar sample

	search_sdss()			- search the SDSS DR8 database using an SQL query for a list of objects
	plot_results()			- plot color-magnitude diagrams and histograms of galaxy type
	fried_plot()			- Reproduce Fig. 3 of Fried et al. (1993), plotting galactic surface density vs. distance

SQL search enabled by Python command line query tool of Tamas Budavari <budavari@jhu.edu>

Required packages:

	sqlcl
	cosmocalc
	pickle

Written by K. Willett, Aug 2011
"""

def search_sdss(blazardir='/Users/willettk/Astronomy/Research/blazars/',blazarfile='plotkin_dr8_upload.txt',savstring='',envsize=1000,ngals=1000, save=False, timing=True):
	
	import sqlcl
	import cosmocalc
	import time

	timestart = time.time()
	
	# Read in the list of blazars from combined catalogs (Plotkin, BZCAT, TeVcat)
	
	pf = open(blazardir+blazarfile,'r')
	
	# Python dictionary for the relevant data to be returned
	
	sdss_data_zns = { 'objid':[], 'ra':[], 'dec':[], 'type':[], 'nchild':[],
		'u':[], 'g':[], 'r':[], 'i':[], 'z':[],
		'p_el':[], 'p_cw':[], 'p_acw':[], 'p_edge':[], 'p_dk':[], 'p_mg':[], 'p_cs':[],
		'blazar_name':[]
		}
	
	sdss_data_zs = { 'objid':[], 'ra':[], 'dec':[], 'type':[], 'nchild':[],
		'u':[], 'g':[], 'r':[], 'i':[], 'z':[],
		'p_el':[], 'p_cw':[], 'p_acw':[], 'p_edge':[], 'p_dk':[], 'p_mg':[], 'p_cs':[],
		'blazar_name':[]
		}
	
	sdss_data_nozoo = { 'objid':[], 'ra':[], 'dec':[], 'type':[], 'nchild':[],
		'u':[], 'g':[], 'r':[], 'i':[], 'z':[],
		'blazar_name':[]
		}
	
	errorarr = {'type':[],'name':[]}
	
	for line in pf:
		#
		blazar_name, blazar_ra, blazar_dec, blazar_z, blazar_type = line.split()
		blazar_z = float(blazar_z)
		#
		# Compute the angular size of 1000 kpc
		#
		platescale = cosmocalc.cosmocalc(z=blazar_z, H0=71, WM=0.27)['PS_kpc'] # kpc per arcsec
		angsize = envsize / platescale / 60.
		#
		# Query the SDSS DR8 catalog
		#
		blazar_query_zns = """ SELECT top %i p.objID, p.ra, p.dec, p.type, p.nchild, p.u, p.g, p.r, p.i, p.z, zns.p_el, zns.p_cw, zns.p_acw, zns.p_edge, zns.p_dk, zns.p_mg, zns.p_cs FROM fGetNearbyObjEq(%f,%f,%f) n JOIN PhotoObj p on p.objid = n.objid JOIN zooNoSpec zns on zns.objid = n.objid ORDER BY p.objID """ % (ngals, float(blazar_ra), float(blazar_dec), angsize)
		blazar_query_zs = """ SELECT top %i p.objID, p.ra, p.dec, p.type, p.nchild, p.u, p.g, p.r, p.i, p.z, zs.p_el, zs.p_cw, zs.p_acw, zs.p_edge, zs.p_dk, zs.p_mg, zs.p_cs FROM fGetNearbyObjEq(%f,%f,%f) n JOIN PhotoObj p on p.objid = n.objid JOIN zooSpec zs on zs.objid = n.objid ORDER BY p.objID """ % (ngals, float(blazar_ra), float(blazar_dec), angsize)
		blazar_query_nozoo = """ SELECT top %i p.objID, p.ra, p.dec, p.type, p.nchild, p.u, p.g, p.r, p.i, p.z FROM fGetNearbyObjEq(%f,%f,%f) n JOIN PhotoObj p on p.objid = n.objid ORDER BY p.objID """ % (ngals, float(blazar_ra), float(blazar_dec), angsize)
		#
                bqtime = time.time()
		#
		queryreturn_zns = sqlcl.query(blazar_query_zns, fmt='csv')
		queryreturn_zs = sqlcl.query(blazar_query_zs, fmt='csv')
		queryreturn_nozoo = sqlcl.query(blazar_query_nozoo, fmt='csv')
		#
		lqtime = time.time()
		#
		# Make sure queries do not exceed 60 per minute (intrinsic limit of database)
		#
		if (lqtime - bqtime) < 3.1:
		        time.sleep(3.1 - (lqtime - bqtime))
		#
		qr_read_zns = queryreturn_zns.read()
		qr_read_zs = queryreturn_zs.read()
		qr_read_nozoo = queryreturn_nozoo.read()
		#
		# Print detections to the screen
		#
		if (qr_read_zns[:26] != 'No objects have been found') and (qr_read_zns[:5] != 'ERROR'):
			print blazar_name
			#
			qrsplit = qr_read_zns.split()
			print qrsplit[1:]
			for n in arange(len(qrsplit[1:]))+1:
				neighbor_data = qrsplit[n].split(',')
				sdss_data_zns['objid'].append(neighbor_data[0])
				sdss_data_zns['ra'].append(neighbor_data[1])
				sdss_data_zns['dec'].append(neighbor_data[2])
				sdss_data_zns['type'].append(neighbor_data[3])
				sdss_data_zns['nchild'].append(neighbor_data[4])
				sdss_data_zns['u'].append(neighbor_data[5])
				sdss_data_zns['g'].append(neighbor_data[6])
				sdss_data_zns['r'].append(neighbor_data[7])
				sdss_data_zns['i'].append(neighbor_data[8])
				sdss_data_zns['z'].append(neighbor_data[9])
				sdss_data_zns['p_el'].append(neighbor_data[10])
				sdss_data_zns['p_cw'].append(neighbor_data[11])
				sdss_data_zns['p_acw'].append(neighbor_data[12])
				sdss_data_zns['p_edge'].append(neighbor_data[13])
				sdss_data_zns['p_dk'].append(neighbor_data[14])
				sdss_data_zns['p_mg'].append(neighbor_data[15])
				sdss_data_zns['p_cs'].append(neighbor_data[16])
				sdss_data_zns['blazar_name'].append(blazar_name)
	
	
		if (qr_read_zs[:26] != 'No objects have been found') and (qr_read_zs[:5] != 'ERROR'):
			print blazar_name
			#
			qrsplit = qr_read_zs.split()
			print qrsplit[1:]
			for n in arange(len(qrsplit[1:]))+1:
				neighbor_data = qrsplit[n].split(',')
				sdss_data_zs['objid'].append(neighbor_data[0])
				sdss_data_zs['ra'].append(neighbor_data[1])
				sdss_data_zs['dec'].append(neighbor_data[2])
				sdss_data_zs['type'].append(neighbor_data[3])
				sdss_data_zs['nchild'].append(neighbor_data[4])
				sdss_data_zs['u'].append(neighbor_data[5])
				sdss_data_zs['g'].append(neighbor_data[6])
				sdss_data_zs['r'].append(neighbor_data[7])
				sdss_data_zs['i'].append(neighbor_data[8])
				sdss_data_zs['z'].append(neighbor_data[9])
				sdss_data_zs['p_el'].append(neighbor_data[10])
				sdss_data_zs['p_cw'].append(neighbor_data[11])
				sdss_data_zs['p_acw'].append(neighbor_data[12])
				sdss_data_zs['p_edge'].append(neighbor_data[13])
				sdss_data_zs['p_dk'].append(neighbor_data[14])
				sdss_data_zs['p_mg'].append(neighbor_data[15])
				sdss_data_zs['p_cs'].append(neighbor_data[16])
				sdss_data_zs['blazar_name'].append(blazar_name)
	
	
		if (qr_read_nozoo[:26] != 'No objects have been found') and (qr_read_nozoo[:5] != 'ERROR'):
			print blazar_name
			#
			qrsplit = qr_read_nozoo.split()
			print qrsplit[1:]
			for n in arange(len(qrsplit[1:]))+1:
				neighbor_data = qrsplit[n].split(',')
				sdss_data_nozoo['objid'].append(neighbor_data[0])
				sdss_data_nozoo['ra'].append(neighbor_data[1])
				sdss_data_nozoo['dec'].append(neighbor_data[2])
				sdss_data_nozoo['type'].append(neighbor_data[3])
				sdss_data_nozoo['nchild'].append(neighbor_data[4])
				sdss_data_nozoo['u'].append(neighbor_data[5])
				sdss_data_nozoo['g'].append(neighbor_data[6])
				sdss_data_nozoo['r'].append(neighbor_data[7])
				sdss_data_nozoo['i'].append(neighbor_data[8])
				sdss_data_nozoo['z'].append(neighbor_data[9])
				sdss_data_nozoo['blazar_name'].append(blazar_name)
	
	
		if qr_read_zns[:5] == 'ERROR':
			errorarr['type'].append('zns')
			errorarr['name'].append(blazar_name)
	
		if qr_read_zs[:5] == 'ERROR':
			errorarr['type'].append('zs')
			errorarr['name'].append(blazar_name)
	
		if qr_read_nozoo[:5] == 'ERROR':
			errorarr['type'].append('nozoo')
			errorarr['name'].append(blazar_name)

	if save is True:

		import pickle

		# Save dictionary to file so it doesn't have to be rerun
		
		output = open('sdss_'+str(envsize)+'kpc_zns'+str(savstring)+'.pkl', 'w')
		pickle.dump(sdss_data_zns, output)
		output.close()
		
		output = open('sdss_'+str(envsize)+'kpc_zs'+str(savstring)+'.pkl', 'w')
		pickle.dump(sdss_data_zs, output)
		output.close()
		
		output = open('sdss_'+str(envsize)+'kpc_nozoo'+str(savstring)+'.pkl', 'w')
		pickle.dump(sdss_data_nozoo, output)
		output.close()

		print "Saved dictionaries to file"
	
	endtime = time.time()

	print 'Total time elapsed is %5.2f min' % ((endtime - starttime)*60)

#################################################################

def plot_results(blazardir='/Users/willettk/Astronomy/Research/blazars/',loaddata=True, envsize = 1000):

	from pylab import *
	import pickle
	
	# Retrieve data from saved file

	if loaddata is True:

		pkl_file = open(blazardir+'sdss_'+str(envsize)+'kpc_zns.pkl', 'r')
		sdss_data_zns = pickle.load(pkl_file)
		pkl_file.close()
	
		pkl_file = open(blazardir+'sdss_'+str(envsize)+'kpc_zs.pkl', 'r')
		sdss_data_zs = pickle.load(pkl_file)
		pkl_file.close()
	
		pkl_file = open(blazardir+'sdss_'+str(envsize)+'kpc_nozoo.pkl', 'r')
		sdss_data_nozoo = pickle.load(pkl_file)
		pkl_file.close()

		print "Loaded data from saved files."

	
	# Plot results

	fig = figure(1)
	clf()
	#
	subplot(221)
	cm_nozoo = plot(array(sdss_data_nozoo['r'],dtype='f'), array(sdss_data_nozoo['u'],dtype='f') - array(sdss_data_nozoo['r'],dtype='f'), 'y*')
	xlabel('r')
	ylabel('u - r')
	xlim(5,35)
	ylim(-17,17)
	#
	subplot(222)
	hist(array(sdss_data_nozoo['type'],dtype='i'),color='g')
	xlabel('type')
	ylabel(r'$N_{gal}$')
	#
	subplot(223)
	cm_zns = plot(array(sdss_data_zns['r'],dtype='f'), array(sdss_data_zns['u'],dtype='f') - array(sdss_data_zns['r'],dtype='f'), 'b*')
	cm_zs = plot(array(sdss_data_zs['r'],dtype='f'), array(sdss_data_zs['u'],dtype='f') - array(sdss_data_zs['r'],dtype='f'), 'r*')
	xlabel('r')
	ylabel('u - r')
	xlim(10,22)
	ylim(-3,12)

	legend([cm_zns,cm_zs],["ZooNoSpec","ZooSpec"],loc=2,numpoints=1)
	leg = gca().get_legend()
	ltext  = leg.get_texts()
	setp(ltext, fontsize='small')
	#
	subplot(224)
	hist(array(sdss_data_zs['p_el'],dtype='f'),color='r')
	hist(array(sdss_data_zns['p_el'],dtype='f'),color='b',histtype='bar')
	xlabel(r'$p_{el}$')
	ylabel(r'$N_{gal}$')
	
	subplots_adjust(wspace=0.4)

	# Reproduce Figure 1 in Fried et al. (1993), determining whether galaxy surface density for BL Lacs changes with galactocentric radius

	

#################################################################

#################################################################

def fried_plot(blazardir='/Users/willettk/Astronomy/Research/blazars/',blazarfile='plotkin_dr8_upload.txt',loaddata=True, envsize = 1000):

	from pylab import *
	import pickle
	import coords
	import cosmocalc

	# Retrieve data from saved file

	if loaddata is True:

		pkl_file = open(blazardir+'sdss_'+str(envsize)+'kpc_zns.pkl', 'r')
		sdss_data_zns = pickle.load(pkl_file)
		pkl_file.close()
	
		pkl_file = open(blazardir+'sdss_'+str(envsize)+'kpc_zs.pkl', 'r')
		sdss_data_zs = pickle.load(pkl_file)
		pkl_file.close()
	
		pkl_file = open(blazardir+'sdss_'+str(envsize)+'kpc_nozoo.pkl', 'r')
		sdss_data_nozoo = pickle.load(pkl_file)
		pkl_file.close()

		print "Loaded data from saved files."

	nz = sdss_data_nozoo
	zns = sdss_data_zns
	zs = sdss_data_zs
	
	# Reproduce Figure 1 in Fried et al. (1993), determining whether galaxy surface density for BL Lacs changes with galactocentric radius

	# Bin by redshift; start with the Fried+ bins at median z of 0.28, 0.65, 0.97 (bin widths of z ~ 0.30)

	# Quantities to plot: 
	#	- distance from host [kpc]

	# Dictionaries have the Python ID number, RA, and dec
	# Reload the original list of blazars; assume closest object to the center (choose cone radius and redshift) is the host galaxy. 
	# For the remainder, convert angular distance to linear. Still assuming no redshift dependence?

	# Make a cut in apparent magnitude at r < 22 (identical to Fried paper)

	r = array(nz['r'],dtype='f')
	rind = where(r < 22)

	sdss_name = array(nz['blazar_name'])[rind[0]]
	sdss_ra = array(nz['ra'])[rind[0]]
	sdss_dec = array(nz['dec'])[rind[0]]

	angseparr = []
	physseparr = []

	angseparr_low = []
	angseparr_mid = []
	angseparr_hi  = []
	physseparr_low = []
	physseparr_mid = []
	physseparr_hi  = []

	b_ra = []
	b_dec = []

	z_low = []
	z_mid = []
	z_hi = []

	pf = open(blazardir+blazarfile,'r')

	for line in pf:
		#
		blazar_name, blazar_ra, blazar_dec, blazar_z = line.split()

		if blazar_name != 'Plotkin_BLLac637':
			blazar_z = float(blazar_z)

			ind1 = list(sdss_name).index(blazar_name)
			ind2 = list(sdss_name).index(blazar_name[:-3]+'%(var)03d' % {'var':float(blazar_name[-3:])+1})

			ra = sdss_ra
			dec = sdss_dec

			b_ra.append(blazar_ra)
			b_dec.append(blazar_dec)

			b_radec = (float(blazar_ra),float(blazar_dec))

			for galaxy in arange(ind1,ind2):

				s_radec = (float(ra[galaxy]),float(dec[galaxy]))

				sdss_coords = coords.Position(s_radec)
				blazar_coords = coords.Position(b_radec)

				sep = sdss_coords.angsep(blazar_coords)
				platescale = cosmocalc.cosmocalc(z=blazar_z, H0=71, WM=0.27)['PS_kpc'] # kpc per arcsec

				sep_arcmin = sep.arcsec() / 60
				sep_kpc = platescale*sep.arcsec()

				angseparr.append(sep_arcmin)
				physseparr.append(sep_kpc)

				if blazar_z < 0.46:
					angseparr_low.append(sep_arcmin)
					physseparr_low.append(sep_kpc)
					z_low.append(blazar_z)
				elif (blazar_z > 0.46) and (blazar_z < 0.81):
					angseparr_mid.append(sep_arcmin)
					physseparr_mid.append(sep_kpc)
					z_mid.append(blazar_z)
				elif (blazar_z > 0.81 and blazar_z < 1.048):
					angseparr_hi.append(sep_arcmin)
					physseparr_hi.append(sep_kpc)
					z_hi.append(blazar_z)


	#	- galaxy surface density [N_gal / sq. arcmin]
	# Bin galaxies into bins of 100 kpc each (galactocentric distance)
	# For the galaxies at each distance from the host object, find the average number of galaxies per square arcminute in that field

	bins_kpc = (arange(10)+1)*100
	bins_arcmin_low = zeros(10)
	bins_arcmin_mid = zeros(10)
	bins_arcmin_hi  = zeros(10)
	
	for b in arange(10):
		bins_arcmin_low[b] = (bins_kpc[b]/cosmocalc.cosmocalc(z=0.28,H0=71, WM=0.27)['PS_kpc'] / 60)
		bins_arcmin_mid[b] = (bins_kpc[b]/cosmocalc.cosmocalc(z=0.65,H0=71, WM=0.27)['PS_kpc'] / 60)
		bins_arcmin_hi[b] = (bins_kpc[b]/cosmocalc.cosmocalc(z=0.97,H0=71, WM=0.27)['PS_kpc'] / 60)
	
	
	hist_low = np.histogram(angseparr_low,bins=bins_arcmin_low)
	hist_mid = np.histogram(angseparr_mid,bins=bins_arcmin_mid)
	hist_hi  = np.histogram(angseparr_hi ,bins=bins_arcmin_hi )
	
	# Assume the plate scale is for the median redshift in each distance class
	
	area_low = np.pi * ((array(bins_arcmin_low[1:10]))**2 - (array(bins_arcmin_low[0:9]))**2)
	area_mid = np.pi * ((array(bins_arcmin_mid[1:10]))**2 - (array(bins_arcmin_mid[0:9]))**2)
	area_hi  = np.pi * ((array(bins_arcmin_hi[1:10]))**2  - (array(bins_arcmin_hi[0:9]))**2)

	density_low = hist_low[0] / area_low
	density_mid = hist_mid[0] / area_mid
	density_hi = hist_hi[0] / area_hi

	fig = figure(2)
	clf()

	###############
	# Sort r-band magnitudes into redshift bins
	###############

	subplot(321)
	hist_physsep = hist(r[:ind2],bins=arange(10,30),color='k')
	hist_physsep = hist(r[rind][:ind2],bins=arange(10,30),color='w')
	xlabel(r'$m_r$')
	ylabel(r'$N_{gal}$')

	subplot(323)
	hist_physsep_low = hist(physseparr_low,color='r')
	hist_physsep_mid = hist(physseparr_mid,color='g')
	hist_physsep_hi  = hist(physseparr_hi,color='b')
	xlabel('Physical separation [kpc]')
	ylabel(r'$N_{gal}$')

	subplot(325)
	plot(b_ra, b_dec, 'kx')
	xlabel('RA')
	ylabel('Dec')
	ylim(-30,90)

	subplot(322)
	errorbar(bins_kpc[1:10], density_low/len(density_low), yerr=sqrt(density_low), color='r', fmt='o')
	xlabel('Physical separation [kpc]')
	ylabel(r'$N_{gal}/sq. arcmin$')
	xlim(0,1000)
	title(r'$\langle z\rangle = %f$' % mean(z_low))
	axhline(y=mean(density_low),color='r',linewidth=2)

	subplot(324)
	errorbar(bins_kpc[1:10], density_mid/len(density_mid), yerr=sqrt(density_mid), color='g', fmt='o')
	xlabel('Physical separation [kpc]')
	ylabel(r'$N_{gal}/sq. arcmin$')
	xlim(0,1000)
	title(r'$\langle z\rangle = %f$' % mean(z_mid))
	axhline(y=mean(density_mid),color='g',linewidth=2)

	subplot(326)
	errorbar(bins_kpc[1:10], density_hi/len(density_hi), yerr = sqrt(density_hi), color='b', fmt='o')
	xlabel('Physical separation [kpc]')
	ylabel(r'$N_{gal}/sq. arcmin$')
	xlim(0,1000)
	title(r'$\langle z\rangle = %f$' % mean(z_hi))
	axhline(y=mean(density_hi),color='b',linewidth=2)

	subplots_adjust(hspace=0.5)

#################################################################

