# script to wade through spectra and log all the vital info
# in the database in preparation for FXCOR. 
# use createTable_EBLM to make the empty table

import pymysql
from astropy.io import fits
import glob as g

# globals
db_tab = "eblm_ids_new"

# set up the database
db = pymysql.connect(host='localhost', db='eblm')
cur = db.cursor()

# get the list of spectra
t = g.glob('*_tn.ms.fits')

# find the reference spectrum
obj = []
for i in t:
    obj.append(fits.open(i)[0].header['OBJECT'])

# loop over spectra and log them
for i in t:
    h = fits.open(i)
    hdr = h[0].header
    d = h[0].data
    image_id = i
    object_name = hdr['OBJECT']
    hjd_mid = hdr['HJD']
    n_traces = d.shape[1]
    pa = hdr['ROTSKYPA']
    utmiddle = hdr['UT-M_E']
    vhelio = hdr['VHELIO']
    vlsr = hdr['VLSR']
    vsun = hdr['VSUN']
    qry = "INSERT INTO %s (image_id, object_name, hjd_mid, n_traces, sky_pa, utmiddle, vhelio, vlsr, vsun) VALUES ('%s', '%s', '%s', %d, %f, '%s', %f, %f, '%s')" % (db_tab, image_id, object_name, hjd_mid, n_traces, pa, utmiddle.replace("T", " "), vhelio, vlsr, vsun)
    print(qry)
    cur.execute(qry)
    db.commit()
    
