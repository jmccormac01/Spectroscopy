#
# IRAF- FXCOR in python
# 
# to do:
#	check the names are ok in the database - go through the logs again :(
#

from pyraf import iraf
import glob as g
from astropy.io import fits
import argparse as ap
import sys, pymysql

# globals
ref_star_name="HD168009"
db_tab="eblm_ids"

# set up the database
db=pymysql.connect(host='localhost',db='eblm')
cur=db.cursor()

def argParse():
	parser=ap.ArgumentParser()
	parser.add_argument('blends',help='Analyse BLENDs (yblends | nblends)', choices=['yblends','nblends'])
	args=parser.parse_args()
	return args

def ImportPackages():
	iraf.noao(_doprint=0)
	iraf.rv(_doprint=0)
	iraf.imred(_doprint=0)
	iraf.kpnoslit(_doprint=0)	
	iraf.ccdred(_doprint=0)
	iraf.astutil(_doprint=0)
	
	iraf.keywpars.setParam('ra','CAT-RA') 
	iraf.keywpars.setParam('dec','CAT-DEC')
	iraf.keywpars.setParam('ut','UT')
	iraf.keywpars.setParam('utmiddl','UT-M_E')
	iraf.keywpars.setParam('exptime','EXPTIME')
	iraf.keywpars.setParam('epoch','CAT-EPOC')
	iraf.keywpars.setParam('date_ob','DATE-OBS')
	iraf.keywpars.setParam('hjd','HJD')
	iraf.keywpars.setParam('mjd_obs','MJD-OBS')
	iraf.keywpars.setParam('vobs','VOBS')
	iraf.keywpars.setParam('vrel','VREL')
	iraf.keywpars.setParam('vhelio','VHELIO')
	iraf.keywpars.setParam('vlsr','VLSR')
	iraf.keywpars.setParam('vsun','VSUN')
	iraf.keywpars.setParam('mode','ql')
	iraf.fxcor.setParam('continu','both')
	iraf.fxcor.setParam('filter','none')
	iraf.fxcor.setParam('rebin','smallest')
	iraf.fxcor.setParam('pixcorr','no')
	iraf.fxcor.setParam('apodize','0.2')
	iraf.fxcor.setParam('function','gaussian')
	iraf.fxcor.setParam('width','INDEF')
	iraf.fxcor.setParam('height','0.')
	iraf.fxcor.setParam('peak','no')
	iraf.fxcor.setParam('minwidt','3.')
	iraf.fxcor.setParam('maxwidt','21.')
	iraf.fxcor.setParam('weights','1.')
	iraf.fxcor.setParam('backgro','0.')
	iraf.fxcor.setParam('window','INDEF')
	iraf.fxcor.setParam('wincent','INDEF')
	iraf.fxcor.setParam('verbose','long')
	iraf.fxcor.setParam('imupdat','no')
	iraf.fxcor.setParam('graphic','stdgraph')
	iraf.fxcor.setParam('interac','yes')
	iraf.fxcor.setParam('autowri','yes')
	iraf.fxcor.setParam('autodra','yes')
	iraf.fxcor.setParam('ccftype','image')
	iraf.fxcor.setParam('observa','lapalma')
	iraf.fxcor.setParam('mode','ql')
	return 0

# get the number of traces
def getNTraces(imlist):
	n_traces=[]
	for i in range(0,len(imlist)):
		qry="SELECT n_traces FROM %s WHERE image_id='%s'" % (db_tab,imlist[i])
		cur.execute(qry)
		for row in cur:
			n_traces.append(row[0])
	# check that we got 1 n_trace for each image
	assert len(n_traces) == len(imlist)	
	return n_traces

# log the out of FXCOR to the database
def logOutput(logfile,image_id,com):
	f=open(logfile).readlines()
	peak_shift_pix=f[35].split()[4]
	correlation_height=f[36].split()[3]
	fwhm_peak_pix=f[37].split()[4]
	fwhm_peak_kms=f[37].split()[6].split('(=')[1]
	relative_velocity_kms=f[39].split()[5]
	observed_velocity_kms=f[40].split()[3]
	helio_velocity_kms=f[41].split()[3]
	qry="UPDATE %s SET peak_shift_pix=%s, correlation_height=%s, fwhm_peak_pix=%s, fwhm_peak_kms=%s, relative_velocity_kms=%s, observed_velocity_kms=%s, helio_velocity_kms=%s, comment='%s' WHERE image_id='%s'" % (db_tab,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,com,image_id)
	cur.execute(qry)
	db.commit()
	return 0

############
### MAIN ###
############

args=argParse()

# get list of directory
t=g.glob('*_tn.ms.fits')

# accepted comments
coms=['SP','DP','BP','NP']

# get n_traces to filter those we are interested in now
n_traces=getNTraces(t)

# import packages from IRAF
imported=ImportPackages()
if str(imported)!='0':
	print "Problem importing IRAF packages, exiting!"
	exit()

# find the reference star
obj=[]
for i in t:
	obj.append(fits.open(i)[0].header['OBJECT'])
try:
	template_loc=obj.index(ref_star_name)
except ValueError:
	print "Reference star %s not found, exiting..." % (ref_star_name)
	sys.exit(1) 

template_image_id=t[template_loc]
print "Ref star %s found as spectrum %s" % (ref_star_name,template_image_id)

# loop over the spectra
for i in range(0,len(t)):
	if args.blends == 'nblends' and n_traces[i] == 1:
		image_id=t[i]
		outfile="%s.out" % (image_id)
		iraf.fxcor(objects=image_id,template=template_image_id,osample='6600-6800',rsample='6600-6800',output=outfile)
		logfile=str(outfile)+".log"
	
		com=''
		while (com.upper() not in coms):
			com=raw_input('Single Peaked = SP - Double Peaked = DP - Broad Peak = BP - No Peak = NP\nEnter comment: ')
			
		done=logOutput(logfile,image_id,com)
	
	# code up the BLENDs later when analysed the first time
	# if args.blends=='yblends' and n_traces[i] > 1:
	#
