import numpy as np
import csv

topdir='/Users/willettk/Astronomy/Research/blazars/blazars_all/'
envsize = '1000'

total = 0

for i in np.arange(1,7):

	blazardir = topdir+'all'+str(i)+'/'

	for j in np.arange(1,6):

		setno = 500*(i-1) + 100*j

		filezns = blazardir+'sdss_'+str(envsize)+'kpc_znscat'+str(setno)+'.pkl'
		filezs = blazardir+'sdss_'+str(envsize)+'kpc_zscat'+str(setno)+'.pkl'
		filenz = blazardir+'sdss_'+str(envsize)+'kpc_nozoocat'+str(setno)+'.pkl'

		pkl_file = open(filezns, 'r')
		sdss_data_zns = pickle.load(pkl_file)
		pkl_file.close()
		
		pkl_file = open(filezs, 'r')
		sdss_data_zs = pickle.load(pkl_file)
		pkl_file.close()
		
		pkl_file = open(filenz, 'r')
		sdss_data_nozoo = pickle.load(pkl_file)
		pkl_file.close()

		if (i == 1) and (j == 1):

			dict_nozoo = sdss_data_nozoo
			dict_zns = sdss_data_zns
			dict_zs = sdss_data_zs

			keys_nozoo = dict_nozoo.keys()
			keys_zns = dict_zns.keys()
			keys_zs = dict_zs.keys()

		else:
			for m in keys_nozoo:
				dict_nozoo[m].extend(sdss_data_nozoo[m])
			for m in keys_zns:
				dict_zns[m].extend(sdss_data_zns[m])
			for m in keys_zs:
				dict_zs[m].extend(sdss_data_zs[m])

		ngals = len(sdss_data_nozoo['z'])
		total += ngals
		print total,' ',setno

dw=csv.DictWriter(open('/Users/willettk/Astronomy/Research/blazars/blazars_all_nozoo.csv','wb'), fieldnames=keys_nozoo)
dw.writer.writerow(dw.fieldnames)
for l in arange(len(dict_nozoo[dw.fieldnames[0]])):
	dw.writer.writerow([str(dict_nozoo[k][l]) for k in dw.fieldnames])

dw=csv.DictWriter(open('/Users/willettk/Astronomy/Research/blazars/blazars_all_zns.csv','wb'), fieldnames=keys_zns)
dw.writer.writerow(dw.fieldnames)
for l in arange(len(dict_zns[dw.fieldnames[0]])):
	dw.writer.writerow([str(dict_zns[k][l]) for k in dw.fieldnames])

dw=csv.DictWriter(open('/Users/willettk/Astronomy/Research/blazars/blazars_all_zs.csv','wb'), fieldnames=keys_zs)
dw.writer.writerow(dw.fieldnames)
for l in arange(len(dict_zs[dw.fieldnames[0]])):
	dw.writer.writerow([str(dict_zs[k][l]) for k in dw.fieldnames])


print ' '
print total

