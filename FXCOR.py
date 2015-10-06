

# IRAF- FXCOR in python

from pyraf import iraf
import glob as g
from astropy.io import fits
import sys

def ImportPackages():
	iraf.noao(_doprint=0)
	iraf.rv(_doprint=0)
	iraf.imred(_doprint=0)
	iraf.kpnoslit(_doprint=0)	
	iraf.ccdred(_doprint=0)
	iraf.astutil(_doprint=0)
	return 0

# import packages from IRAF
imported=ImportPackages()

if str(imported)!='0':
	print "Problem importing IRAF packages, exiting!"
	exit()

# get list of directory
t=g.glob('*_tn.ms.fits')

# find the reference star
ref_star="HD168009"
object_id=[]
for i in t:
	h=fits.open(i)[0].header['OBJECT']
	object_id.append(h)

try:
	template_loc=object_id.index(ref_star)
except ValueError:
	print "Reference star %s not found, exiting..." % (ref_star)
	sys.exit(1) 

template=t[template_loc]
print "Ref star %s found as spectrum %s" % (ref_star,template)

iraf.keywpars.setParam('ra','RA') 
iraf.keywpars.setParam('dec','DEC')
iraf.keywpars.setParam('ut','UT')
iraf.keywpars.setParam('utmiddl','UTMIDDLE')
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


# need to write a function to get the following numbers from the log file:
# peak_shift_pix       
# correlation_height   
# fwhm_peak_pix        
# fwhm_peak_kms        
# relative_velocity_kms
# observed_velocity_kms
# helio_velocity_kms   
# and log them to the database where image_id = t[i]

for i in range(0,len(t)):
	object=t[i]
	outfile="%s.out" % (object)
	iraf.fxcor(objects=object,template=template,osample='6600-6800',rsample='6600-6800',output=outfile)
	logfile=str(outfile)+".log"
	f=file(logfile,'r').readlines()
	vshift=f[-9].split('=')[1].split(' ')[1]
	print "Target: %s\tMeasured shift: %s km/s\n" % (object_id[i],vshift)

