# script to wade through spectra and log all the vital info
# in the database in preparation for FXCOR. 
# use createTable_EBLM to make the empty table

import pymysql
from astropy.io import fits
import glob as g
import matplotlib.pyplot as pl

# globals
ref_star_name="HD168009"
db_tab="eblm_ids"

# set up the database
db=pymysql.connect(host='localhost',db='eblm')
cur=db.cursor()

# get the list of spectra
t=g.glob('*_tn.ms.fits')

# find the reference spectrum
obj=[]
for i in t:
	obj.append(fits.open(i)[0].header['OBJECT'])
try:
	template_loc=obj.index(ref_star_name)
except ValueError:
	print "Reference star %s not found, exiting..." % (ref_star_name)
	sys.exit(1) 
template_image_id=t[template_loc]

# loop over spectra and log them
for i in t:
	h=fits.open(i)
	hdr=h[0].header
	d=h[0].data
	
	image_id=i
	object_name=hdr['OBJECT']
	hjd_mid=hdr['HJD']
	n_traces=d.shape[1]
	utmiddle=hdr['UT-M_E']
	
	qry="INSERT INTO %s (image_id,object_name,template_image_id,ref_star_name,hjd_mid,n_traces,utmiddle) VALUES ('%s','%s','%s','%s',%f,%d,'%s')" % (db_tab,image_id,object_name,template_image_id,ref_star_name,hjd_mid,n_traces,utmiddle.replace("T"," "))
	cur.execute(qry)
	db.commit()
	
        