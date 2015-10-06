#
# python program to split up all the INT spectra into their own folders
#

import os
import glob as g
from astropy.io import fits
import matplotlib.pyplot as pl

sort_yn = raw_input("Sort? (e.g. y): ")
if str(sort_yn) == 'y':
	
	t=g.glob('*tn.ms.fits')
	for i in range(0,len(t)):		
		object=fits.open(t[i])[0].header['CAT-NAME']
		print object
		if os.path.exists(object) == False:
			os.mkdir(object)	
		line="%s/%s" % (object,t[i])	
		os.rename(t[i],line)
		print line
	
compare_yn = raw_input("Compare? (e.g. y): ")
if str(compare_yn) == 'y':

	t=g.glob('*tn.ms.fits')
	h=fits.open(t[0])
	ref_spec=h[0].data[0][0]
	
	for i in range(1,len(t)):
		h2=fits.open(t[i])
		check_spec=h2[0].data[0][0]
		
		dif_spec=ref_spec-check_spec		
		pl.title('Ref Spec - Check Spec')
		pl.plot(dif_spec)
		pl.show()
		

