"""
This program aims to overplot images from different surveys and/or from our work.

HISTORY:
	2016.08.30:
        - first shot

RUN:

"""

import os
import numpy as np
from astroquery.skyview import SkyView

from astropy import units as u 
from astropy import coordinates
from astropy.io import fits

from astroquery.alma import Alma, utils

import matplotlib.pyplot as plt
import aplpy


__author__="Ridlo W. Wibowo @ ITB"
__version__="0.1"



class CompareImage:
	def __init__(self):
		pass

	def read_data_points(self, ifile):
		points = np.genfromtxt(ifile, names=True)
		ra = points['alpha']
		dec = points['dec']
		fwhm = points['FWHM']

		return ra, dec, fwhm


	def get_images(self, pos, radius, surveys, savefits=True):
		images = SkyView.get_images(position=pos, coordinates='icrs', survey=surveys, radius=radius*u.arcmin, pixels=1000)
		if savefits:
			for i, s in enumerate(surveys):
				images[i].writeto(s+'.fits', clobber=True) # overwrite existing file

		return images


	def overplot_points_to_image(self, data, image, show=True, savefig=True, ofile='overplot.png'):
		ra, dec, fwhm = data
		#image = fits.open(fits_file)

		fig = aplpy.FITSFigure(image)
		
		fig.show_grayscale()
		#fig.show_colorscale(cmap='gist_heat')
		#fig.show_contour(colors='white')

		fig.show_markers(ra, dec, edgecolor='blue', facecolor='none', marker='o', s=fwhm*200, alpha=1.0)
		if savefig:
			fig.save(ofile)

		if show:
			plt.show(); plt.close()

		image.close()



if __name__ == "__main__":
	ifile = 'point_source.dat'
	position = coordinates.SkyCoord('05h22m57.98465s', '-36d27m30.8509s', frame='icrs') 
	
	CI = CompareImage()
	
	points = CI.read_data_points(ifile)

	surveys = ['DSS', '2MASS-K', '2MASS-J', 'WISE 22', 'WISE 3.4']
	
	images = CI.get_images(position, 1.0, surveys)

	for i, s in enumerate(surveys):
		# images = fits.open(s+'.fits')
		CI.overplot_points_to_image(points, images[i], ofile=s+'_and_pointSources.png')



	# surv = ['Fermi 5', 'HRI', 'DSS']
	# etaimages = SkyView.get_images(position='Eta Carinae', survey=surv)
	# for i, s in enumerate(surv):
	# 	etaimages[i].writeto(s+'.fits', clobber=True)
	# 	fig = aplpy.FITSFigure(etaimages[i])
	# 	fig.show_grayscale()
	# 	fig.save('eta_'+s+'.png')