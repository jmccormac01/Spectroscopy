# script to wade through spectra an pull out blends
# use createTable_EBLM to make the empty table

import pymysql
from astropy.io import fits
import glob as g

t=g.glob('*.fits')

for i in t:
	h=fits.open(i)[0].data
	print h.shape
	