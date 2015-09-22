

# IRAF- FXCOR in python

from pyraf import iraf
import commands
import pyfits

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
templist=commands.getoutput('ls *_tn.ms.fits').split('\n')

if str(templist[0]) == 'ls: *_tn.ms.fits: No such file or directory':
	templist=commands.getoutput('ls *_tn.ms.fits').split('\n')

templist2=commands.getoutput('ls *_tn-S.ms.fits').split('\n')

template=templist2[0]

hdulist=pyfits.open(template)
prihdr=hdulist[0].header

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


for i in range(0,len(templist)):
	object=templist[i]
	
	target=str(templist[i])
	
	outfile=str(target)+"-"+str(i)
	
	iraf.fxcor(objects=object,template=template,osample='6600-6800',rsample='6600-6800',output=outfile)
	
	logfile=str(outfile)+".log"
	
	f=file(logfile,'r')
	s=f.readlines()
	f.close()
	
	vshift=s[-9].split('=')[1].split(' ')[1]

	print "Target: " + str(target) + "\tMeasured shift: " + str(vshift) + " Km/sec\n"

