from pyraf import iraf

def ImportPackages():
	iraf.noao(_doprint=0)
	iraf.rv(_doprint=0)
	iraf.keywpars.setParam('ra','CAT-RA') 
	iraf.keywpars.setParam('dec','CAT-DEC')
	# must use UT as UT mid exposure for RV correct or it 
	# messes up the HJD values by a lot
	iraf.keywpars.setParam('ut','UT-M_E')
	iraf.keywpars.setParam('utmiddl','')
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
	iraf.rvcorrect.setParam('imupdat','yes')
	iraf.rvcorrect.setParam('observa','lapalma')

	return 0

ImportPackages()
iraf.rvcorrect(images='*.fits')