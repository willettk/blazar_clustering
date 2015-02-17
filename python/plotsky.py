from matplotlib import pyplot as plt
from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u

blazardir = '/Users/willettk/Astronomy/Research/blazars/'
with fits.open('%s/newdata/blazars_unique.fits' % blazardir) as f:
    data = f[1].data

ra = coord.Angle(data['RA2000']*u.degree)
ra = ra.wrap_at(180*u.degree)
dec = coord.Angle(data['DEC2000']*u.degree)

fig = plt.figure(figsize=(8,6))

ax1 = fig.add_subplot(111,projection='mollweide',axisbg='0.90')
ax1.scatter(ra.radian,dec.radian,s=10)
ax1.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
ax1.grid(True)

plt.show()

fig.savefig('%s/map.pdf' % blazardir)
