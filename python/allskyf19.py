from kapteyn import maputils, tabarray
from matplotlib import pyplot as plt
import numpy
import sys

##################################################

def cylrange(epsilon):
   X = numpy.arange(0,400.0,30.0); 
   # Replace last two (dummy) values by two values around 180 degrees   
   X[-1] = 180.0 - epsilon
   X[-2] = 180.0 + epsilon
   return X

##################################################

def doplot(frame, annim, grat, title, 
           lon_world=None, lat_world=None, 
           lon_constval=None, lat_constval=None,
           lon_fmt=None, lat_fmt=None,
           markerpos=None, 
           plotdata=False, perimeter=None, drawgrid=None, 
           smallversion=False, addangle0=0.0, addangle1=0.0, 
           framebgcolor=None, deltapx0=0.0, deltapy0=0.0,
           deltapx1=0.0, deltapy1=0.0,
           labkwargs0={'color':'r'}, labkwargs1={'color':'b'}):
# Apply some extra settings
   
   if framebgcolor != None:
      frame.set_axis_bgcolor(framebgcolor)

   if lon_constval == None:
      lon_constval = 0.0    # Reasonable for all sky plots
   if lat_constval == None:
      lat_constval = 0.0    # Reasonable for all sky plots
   if lon_fmt == None:
      lon_fmt = 'Dms'
   if lat_fmt == None:
      lat_fmt = 'Dms'
   # Plot labels inside graticule if required
   ilabs1 = grat.Insidelabels(wcsaxis=0, 
                        world=lon_world, constval=lat_constval, 
                        deltapx=deltapx0, deltapy=deltapy0, 
                        addangle=addangle0, fmt=lon_fmt, **labkwargs0)
   ilabs2 = grat.Insidelabels(wcsaxis=1, 
                        world=lat_world, constval=lon_constval, 
                        deltapx=deltapx1, deltapy=deltapy1, 
                        addangle=addangle1, fmt=lat_fmt, **labkwargs1)

   # Plot just 1 pixel c.q. marker
   if markerpos != None:
      annim.Marker(pos=markerpos, marker='o', color='red')

   if drawgrid:
      pixellabels = annim.Pixellabels(plotaxis=(2,3))

   # Plot the title
   if smallversion:
      t = frame.set_title(title, color='g', fontsize=10)
   else:
      t = frame.set_title(title, color='g', fontsize=13, linespacing=1.5)
   titlepos = 1.02
   t.set_y(titlepos)
   annim.plot()
   
   # Plot alternative borders. Do this after the graticule is plotted
   # Only then you know the frame of the graticule and plotting in that
   # frame will overwrite graticule lines so that the borders look better
   if perimeter != None:
      p = plt.Polygon(perimeter, facecolor='#d6eaef', lw=2)
      frame.add_patch(p)     # Must be in frame specified by user
      Xp, Yp = zip(*perimeter)
      grat.frame.plot(Xp, Yp, color='r')

   annim.interact_toolbarinfo()
   plt.show()

###################################################################

# Set defaults 

def make_allskyf19():

	titlepos = 1.02
	drawgrid = False
	grat = None
	plotbox = (0.1,0.05,0.8,0.8)
	epsilon = 0.0000000001
	
	fig = plt.figure(2,figsize=(4,4))
	fig.clf()
	frame = fig.add_axes(plotbox)
	title = "BzCaT BL Lac sky distribution"
	header = {'NAXIS'  : 2, 
	          'NAXIS1' : 100, 
	          'NAXIS2' : 80,
	          'CTYPE1' : 'RA---AIT',
	          'CRVAL1' : 0.0, 
		  'CRPIX1' : 50, 
		  'CUNIT1' : 'deg', 
		  'CDELT1' : -4.0,
	          'CTYPE2' : 'DEC--AIT',
	          'CRVAL2' : 0.0, 
		  'CRPIX2' : 40, 
		  'CUNIT2' : 'deg', 
		  'CDELT2' : 4.0
	         }
	X = cylrange(epsilon)
	Y = numpy.arange(-60,90,30.0)
	f = maputils.FITSimage(externalheader=header)
	annim = f.Annotatedimage(frame)
	grat = annim.Graticule(axnum= (1,2), wylim=(-90,90.0), wxlim=(0,360),
	                       startx=X, starty=Y)
	grat.setp_lineswcs0(0, lw=2)
	grat.setp_lineswcs1(0, lw=2)
	lat_world = [-60,-30, 0, 30, 60]
	
	# Remove the left 180 deg and print the right 180 deg instead
	
	w1 = numpy.arange(0,151,30.0)
	w2 = numpy.arange(210,360,30.0)
	lon_world = numpy.concatenate((w1, w2))
	labkwargs0 = {'color':'r', 'va':'bottom', 'ha':'right'}
	labkwargs1 = {'color':'b', 'va':'bottom', 'ha':'right'}
	
	# Load the positions of blazars from text file
	
	fn = '/Users/willettk/Astronomy/Research/blazars/bzcat_bllac_dec_upload.txt'
	fn = '/Users/willettk/Astronomy/Research/blazars/allblazars_remdup.cat'
	xp, yp = annim.positionsfromfile(fn, 's', cols=[1,2])
	annim.Marker(x=xp, y=yp, mode='pixels', marker='x', color='g')
	
	# Set font sizes and colors of the grid labels
	
	labkwargs0={'color':'r', 'fontsize':16, 'va':'baseline', 'ha':'right'}
	labkwargs1={'color':'b', 'fontsize':16, 'va':'center', 'ha':'right'}
	
	# Make the plot
	
	doplot(frame, annim, grat, title,
	       lon_world=lon_world, lat_world=lat_world,
	       labkwargs0=labkwargs0, labkwargs1=labkwargs1)
	
