# utility script to fill in the new sky_pa column in the database.
# can be modified to fill in any missing column using header values.
import pymysql
from astropy.io import fits
import glob as g

# globals
db_tab="eblm_ids"

# set up the database
db=pymysql.connect(host='localhost',db='eblm')
cur=db.cursor()

# get the list of spectra
t=g.glob('*_tn.ms.fits')

# loop over spectra and log the missing PA value
for i in t:
	pa=fits.open(i)[0].header['ROTSKYPA']
	qry="UPDATE %s SET sky_pa=%f WHERE image_id='%s'" % (db_tab,pa,i)
	print(qry)
	cur.execute(qry)
	db.commit()
