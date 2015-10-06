########################################################################################
################################## SPECTOR.py ##########################################
########################################################################################
#
# 	INT/IDS Spectral Reduction v2.0: A Script for complete reduction of INT/IDS data
#		
#								James McCormac
#	
#	Before running Spector:
#		- 	Remove all files that should not be included in reduction
#			i.e. bad flats, data from other central wavelengths or gratings etc.
#			Place all this data in a folder called 'junk'. Use *.int logfile to
#			help weed out unwanted files
#		-	Adjust the settings in the **EDIT HERE** section below
#
#	Assumptions:
#		-	This code was written for a particular observing strategy:
#				Acquire, Arc, Science, Arc, REPEAT...
#			Therefore the script is looking for this structure in the data when matching 
#			arcs to spectra. This can be easily changed by modifying the relevant section
# 
# 	What Spector does:
#		- 	Prepares the working environment, making copies of original data
#		-	Loads IRAF
#		- 	Strips multi extension fits to single extension
#		- 	Renames files for easier viewing:
#				i_* = integration (spectra)
#				a_* = arc
#				b_* = bias
#				f_* = flat
#		-	Debiases and flat field corrects the data
#		- 	Tidies up calibration frames into 'calibs' folder
#		- 	Trims the spectra based on INT's heavily vignetted region, see TRIMSEC below
#		-	Sets JD and Airmass (Note JD needs updating to Eastman 2010's method)
#		- 	Interactive spectral analysis with KPNOSLIT in IRAF 	
#		- 	Normalises the spectra	
#		- 	Tidies up working directory after reductions are finished
#
#	What Spector does not do:
#		- 	Flux calibration
#
#	What you will end up with:
#		-	A set of reduced, wavelength calibrated, normalised 1D spectra with the 
#			byproducts from each step located in subdirectories within the 'reduced'
#			folder. 
#
#	Who to annoy if you have questions:
#		jmccormac001@gmail.com
#
# Changes:
# 	2015-10-06: Removed setdJD and setAirmass as they were deeply broken!
#

# modules imported
from numpy import *                       
import os, sys, string, math
import pyfits
import time
import commands

#################################
########### EDIT HERE ###########
#################################

lineList_location = "/Users/James/Documents/LineLists/cuar_cune.dat"
loginCl_location = "/Users/James/"

GAIN="1.03"
RDNOISE="3.48"
TRIMSEC='[1:360,1000:3000]'

#################################

def Makedirs():
	if os.path.isdir('original') == False:
		os.mkdir('original')
	if os.path.isdir('reduced') == False:
		os.mkdir('reduced')
	return 0

def CopyFiles():
	currentdir=os.getcwd()
	contents=os.listdir(currentdir)
	madedirs=Makedirs()
			
	for i in range(0,len(contents)):
		f1=contents[i]
		
		if str(f1) != 'junk':
			f2="original/"+str(contents[i])
			os.rename(f1, f2)
			print str(f1) + " moved to 'original' folder"	
		
	os.chdir('original/')
	templist=commands.getoutput('ls r*.fit').split('\n')
		
	for i in range(0,len(templist)):
		f3=templist[i]
		f4="../reduced/"+str(templist[i])
		os.system("cp %s %s" % (f3, f4))
		print str(f3) + " copied to reduction folder"
	
	# move the login and pyraf files into reduced
	os.system("mv pyraf/ ../")
	
	# leave off in reduced dir for working
	os.chdir('../reduced/')
	return 0

def ImportPackages():
	iraf.noao(_doprint=0)
	iraf.imred(_doprint=0)
	iraf.kpnoslit(_doprint=0)	
	iraf.ccdred(_doprint=0)
	iraf.astutil(_doprint=0)
	return 0

def StripFiles():
	os.mkdir("unstripped")
	print "Directory Created: 'unstripped'"
	
	ob,bi,fl,ar=0,0,0,0
	
	templist=commands.getoutput('ls r*').split('\n')
		
	for i in range(0,len(templist)):
		image=str(templist[i]+"[1]")
		header=str(templist[i])
		hdulist=pyfits.open(header)
		prihdr=hdulist[0].header
		object=prihdr['IMAGETYP']
		
		if str(object)=='object':
			image2=str("i_s_"+str(templist[i]))
			ob=ob+1
		if str(object)=='zero':
			image2=str("b_s_"+str(templist[i]))
			bi=bi+1
		if str(object)=='flat':
			image2=str("f_s_"+str(templist[i]))
			fl=fl+1
		if str(object)=='arc':
			image2=str("a_s_"+str(templist[i]))
			ar=ar+1
	
		iraf.imcopy(input=image,output=image2)
		print image
		print image2
		f2="unstripped/"+str(templist[i])
		os.rename(templist[i], f2)
	
	print "\n\tSpectra: " + str(ob)
	print "\tArcs: " + str(ar)
	print "\tFlats: " + str(fl)
	print "\tBiases: " + str(bi) + "\n"
	return 0

def Zerocombine():
	iraf.zerocombine(input="b_s_r*",output="Zero", combine="median", reject="minmax", ccdtype="zero", process="no", delete="no", clobber="no", scale="none", statsec="", nlow="0", nhigh="1", nkeep="1" , mclip="yes", lsigma="3.", hsigma="3.", rdnoise=RDNOISE, gain=GAIN, snoise="0.", pclip="-0.5", blank="0.")
	return 0

def ProcFlats():
	iraf.ccdproc(images="f_s_r*",output="", ccdtype="", max_cache="0", noproc="no", fixpix="no", overscan="no",trim="no", zerocor="yes", darkcor="no", flatcor="no", illumcor="no", fringecor="no",readcor="no", scancor="no", readaxis="line", fixfile="", biassec="", trimsec="",zero="Zero.fits", dark="", flat="", illum="", fringe="", minreplace="1.",scantype="shortscan", nscan="1", interactive="no", function="chebyshev", order="1", sample="*", naverage="1", niterate="1", low_reject="3.", high_reject="3.", grow="1.")
	return 0
	
def Flatcombine():	
	iraf.flatcombine(input="f_s_r*",output="Flat", combine="median", reject="avsigclip", ccdtype="flat", process="no", subsets="yes", delete="no", clobber="no", scale="mode", statsec="", nlow="1", nhigh="1", nkeep="1", mclip="yes", lsigma="3.", hsigma="3.", rdnoise=RDNOISE, gain=GAIN,snoise="0.", pclip="-0.5", blank="1.")
	return 0

def ProcArcs():
	iraf.ccdproc(images="a_s_r*",output="", ccdtype="", max_cache="0", noproc="no", fixpix="no", overscan="no", trim="no", zerocor="yes", darkcor="no", flatcor="yes", illumcor="no", fringecor="no", readcor="no", scancor="no", readaxis="line", fixfile="", biassec="", trimsec="",zero="Zero.fits", dark="", flat="Flat.fits", illum="", fringe="",minreplace="1.", scantype="shortscan", nscan="1", interactive="no", function="chebyshev", order="1", sample="*", naverage="1", niterate="1", low_reject="3.", high_reject="3.", grow="1.")
	return 0
	
def ProcSpec():
	iraf.ccdproc(images="i_s_r*",output="", ccdtype="", max_cache="0", noproc="no", fixpix="no", overscan="no",trim="no", zerocor="yes", darkcor="no", flatcor="yes", illumcor="no", fringecor="no", readcor="no", scancor="no", readaxis="line", fixfile="", biassec="", trimsec="",zero="Zero.fits", dark="", flat="Flat.fits", illum="", fringe="",minreplace="1.", scantype="shortscan", nscan="1", interactive="no",function="chebyshev", order="1", sample="*", naverage="1", niterate="1", low_reject="3.", high_reject="3.", grow="1.")
	return 0

def CleanCalibs():
	
	os.mkdir('calibs')
	
	templist=commands.getoutput('ls f_s_r*').split('\n')
	templist2=commands.getoutput('ls b_s_r*').split('\n')
	
	for i in range(0,len(templist)):
		f3="calibs/"+str(templist[i])
		os.rename(templist[i], f3)
	for i in range(0,len(templist2)):
		f4="calibs/"+str(templist2[i])
		os.rename(templist2[i], f4)
		
	os.rename("Zero.fits","calibs/Zero.fits")
	os.rename("Flat.fits","calibs/Flat.fits")
	return 0

def TrimFiles():
	
	section=TRIMSEC
	os.mkdir("untrimmed")
	print "Directory Created: 'untrimmed'"
	
	templist=commands.getoutput('ls i_s_r*').split('\n')
	templist2=commands.getoutput('ls a_s_r*').split('\n')
	
	for i in range(0,len(templist)):
		image=str(templist[i])+section
		image2=str(templist[i].split('.')[0]) + "_t"
		iraf.imcopy(input=image,output=image2)
		f2="untrimmed/"+str(templist[i])
		os.rename(templist[i], f2)
		
	for i in range(0,len(templist2)):
		image=str(templist2[i])+section
		image2=str(templist2[i].split('.')[0]) + "_t"
		iraf.imcopy(input=image,output=image2)
		f3="untrimmed/"+str(templist2[i])
		os.rename(templist2[i], f3)
	return 0

#def SetJD():
#	iraf.setjd(images="*.fits", date="date", time='utmiddle', exposur="exptime", ra="ra", dec="dec", epoch="cat-epoc", jd="jd", hjd="hjd", ljd="ljd", utdate="yes", uttime="yes", listonl="no", mode="ql")
#	return 0
	
#def SetAirmass():
#	iraf.setairmass(images="*.fits", ra="cat-ra", dec="cat-dec", equinox="cat-epoc", st="st", ut="utstart", date="date", exposur="exptime", airmass="airmass", utmiddl="utmiddle", show="yes", update="yes", overrid="yes",mode="ql")
#	return 0

def Analyse():
	# make a list of the science images to be analysed
	templist=commands.getoutput('ls i_s*').split('\n')
			
	# apall parameters
	iraf.apall.setParam('format','multispec') 	
	iraf.apall.setParam('interac','yes')     				
	iraf.apall.setParam('find','yes') 
	iraf.apall.setParam('recen','yes') 
	iraf.apall.setParam('resize','yes') 
	iraf.apall.setParam('trace','yes') 
	iraf.apall.setParam('fittrac','yes') 
	iraf.apall.setParam('extract','yes') 
	iraf.apall.setParam('extras','yes') 
	iraf.apall.setParam('review','yes') 
	
	iraf.apall.setParam('line','INDEF') 
	iraf.apall.setParam('nsum','12') 
	
	iraf.apall.setParam('lower','-6') 
	iraf.apall.setParam('upper','6') 
	
	iraf.apall.setParam('b_funct','chebyshev') 
	iraf.apall.setParam('b_order','1') 
	iraf.apall.setParam('b_sampl','-25:-15,15:25') 
	iraf.apall.setParam('b_naver','-100') 
	iraf.apall.setParam('b_niter','0') 
	iraf.apall.setParam('b_low_r','3') 
	iraf.apall.setParam('b_high','3') 
	iraf.apall.setParam('b_grow','0') 
	
	iraf.apall.setParam('width','10') 
	iraf.apall.setParam('radius','10') 
	iraf.apall.setParam('threshold','0') 
	
	iraf.apall.setParam('nfind','1') 
	
	iraf.apall.setParam('t_nsum','10') 
	iraf.apall.setParam('t_step','10')
	iraf.apall.setParam('t_nlost','3')
	iraf.apall.setParam('t_funct','spline3')
	iraf.apall.setParam('t_order','3')
	
	iraf.apall.setParam('backgro','fit')
	iraf.apall.setParam('skybox','1')
	iraf.apall.setParam('weights','variance')
	iraf.apall.setParam('pfit','fit1d')
	iraf.apall.setParam('clean','yes')
	iraf.apall.setParam('saturat','40000')
	iraf.apall.setParam('readnoi','4.38')
	iraf.apall.setParam('gain',GAIN)
	iraf.apall.setParam('lsigma','4.0')
	iraf.apall.setParam('usigma','4.0')
	iraf.apall.setParam('nsubaps','1')
	
	iraf.apall.saveParList(filename="apall.pars")
	
	# make reference arc for reidentify
	refarc=[]
	
	for i in range(0,len(templist)):
		image=templist[i]		
	
		hdulist=pyfits.open(templist[i])
		prihdr=hdulist[0].header
		object=prihdr['CAT-NAME']
	
		# extract the object spectrum
		print "Extracting spectrum of " + str(object) + " from image " + str(image)
		print "Check aperture and background. Change if required"
		print "AP: m = mark aperture, d = delete aperture"
		print "SKY: s = mark sky, t = delete sky, f = refit"
		print "q = continue"
		iraf.apall(input=image)
		print "Spectrum extracted!"
	
		# find the arcs either side of the object
		arc1="a_s_r" + str(int(templist[i].split('_')[2].split('r')[1])-1) + "_t.fits"
		arc2="a_s_r" + str(int(templist[i].split('_')[2].split('r')[1])+1) + "_t.fits"
		
		# note the reference arc
		if i==0:
			refarc.append(arc1)
		
		# predict the arc names
		print "\nPredicting arcs names..."
		print "Arc1: " + str(arc1)
		print "Arc2: " + str(arc2)
		
		# setup a reference filename for the arc conditions in database
		reffile=templist[i].split('.fits')[0]
		
		# extract the arcs
		print "\nExtracting arcs under the same conditions..."
		iraf.apall(input=arc1,reference=reffile,recente="no",trace="no",backgro="no",interac="no")
		print "Arc1 extracted"
		iraf.apall(input=arc2,reference=reffile,recente="no",trace="no",backgro="no",interac="no")
		print "Arc2 extracted"
		
		# get a list of the extracted arcs and objects
		templist2=commands.getoutput('ls a_s*.ms.fits').split('\n')[-2:]
		templist3=commands.getoutput('ls i_s*.ms.fits').split('\n')[-1]
		
		print "\nArcs extracted to:"
		print templist2[0]
		print templist2[1]
		
		if i==0:
			print "\nIdentify arc lines:"
			print "Enter the following in the splot window"
			print "\t:thres 500"
			print "\t:order 3"
			print "\tfwidth 2"
			print "Select 3-5 arc lines from line atlas"
			print "Press 'm' to mark, then enter wavelength"
			print "Then press 'l' to automatically ID the other lines"
			print "Press 'f' to fit the dispersion correction"
			print "Use 'd' to remove bad points, 'f' to refit"
			print "'q' from fit, then 'q' from identify to continue\n"
					
		if i==0:
			arcid=templist2[0]
			refarc.append(arcid)
			iraf.identify(images=arcid,coordlist=lineList_location)
			iraf.reidentify(reference=refarc[1],images=templist2[1])
		if i>0:
			print "\nReidentifying arclines from " + str(refarc[1])
			iraf.reidentify(reference=refarc[1],images=templist2[0])
			iraf.reidentify(reference=refarc[1],images=templist2[1])	
		
		# add the refspec keywords to the image header for dispcor
		refspec1=str(templist2[0].split(".fits")[0]) + " 0.5"
		refspec2=str(templist2[1].split(".fits")[0]) + " 0.5"
		
		print "\nREFSPEC1: " + str(refspec1)
		print "REFSPEC2: " + str(refspec2)
		iraf.hedit(images=templist3,fields="REFSPEC1",value=refspec1,add="yes",verify="no",show="yes",)
		iraf.hedit(images=templist3,fields="REFSPEC2",value=refspec2,add="yes",verify="no",show="yes",)
		print "Headers updated!\n"
		
		# apply the dispersion correction
		print "Applying the dispersion correction"
		iraf.dispcor(input=templist3,output=templist3,lineari="yes",databas="database",table="")
		print "Correction applied!"
		
		# normalize the spectrum using continuum
		normspec=str(templist3.split('.ms')[0]) + "n.ms.fits"
		iraf.continuum(input=templist3, output=normspec,logfile="logfile",interac="yes", functio="spline3", order="5", niterat="10", markrej="yes")
		print "\n\n"
	return 0 

def RoundUpFiles():

	templist=commands.getoutput('ls *_tn.ms.fits').split('\n')
	if os.path.exists('spectra_r') == False:
		os.mkdir('spectra_r')
	
	for i in range(0,len(templist)):
		f1="spectra_r/"+str(templist[i])
		os.rename(templist[i], f1)
	
	# tar up the results	
	os.system('tar cvzf spectra_r.tgz spectra_r/')
	return 0
		
def RoundUpArcs():

	templist=commands.getoutput('ls i*t.ms.fits').split('\n')
	templist2=commands.getoutput('ls i*t.fits').split('\n')
	templist3=commands.getoutput('ls a*t.ms.fits').split('\n')
	templist4=commands.getoutput('ls a*t.fits').split('\n')
	
	if os.path.exists('unnormalized') == False:
		os.mkdir('unnormalized')
	
	for i in range(0,len(templist)):
		f1="unnormalized/"+str(templist[i])
		os.rename(templist[i], f1)
	for i in range(0,len(templist2)):
		f2="unnormalized/"+str(templist2[i])
		os.rename(templist2[i], f2)
	for i in range(0,len(templist3)):
		f3="unnormalized/"+str(templist3[i])
		os.rename(templist3[i], f3)
	for i in range(0,len(templist4)):
		f4="unnormalized/"+str(templist4[i])
		os.rename(templist4[i], f4)
	return 0

########################################################################################
##################################### MAIN #############################################
########################################################################################

print '\n\n---------------------------------------------------------'
print '------------- Spectroscopy by Spector.py ----------------'
print '---------------------------------------------------------\n'

#get time of start
starttime=time.ctime()
print 'Program started on '+starttime

# load IRAF from directory where the login.cl file is
here=os.getcwd()
os.chdir(loginCl_location)
from pyraf import iraf
os.chdir(here)

time.sleep(2)

print "Script set up and junk removed?"
junk_yn = raw_input("(e.g. y): ")
if str(junk_yn)=='n':
	print "Edit the EDIT HERE section of the script and/or"
	print "go back and weed out the bad frames, exiting!"
	exit()
	
print "Setup work environment?"
we_yn = raw_input("(e.g y): ")
if str(we_yn) == 'y':
	# set up working environment and move into reduced 
	# folder for working
	print "Setting up work environment..."
	time.sleep(1)
	copiedfiles=CopyFiles()
	if copiedfiles != 0:
		print "Problem with the environment setup, exiting!"
		exit()
	time.sleep(2)

print "Load IRAF packages?"
pack_yn = raw_input("(e.g y): ")
if str(pack_yn) == 'y':
	# load packages
	print "Loading IRAF packages..."
	time.sleep(1)
	loaded=ImportPackages()
	if loaded != 0:
		print "Problem loading IRAF packages, exiting!"
		exit()

print "Strip and rename files?"
strip_yn = raw_input("(e.g y): ")
if str(strip_yn) == 'y':
	# strip files and rename
	# mv unstripped files to separate directory
	print "Renaming Files..."
	time.sleep(1)
	stripped=StripFiles()
	if stripped !=0:
		print "Problem striping and renaming files, exiting!"
		exit()

print "Reduce the data?"
reduce_yn = raw_input("(e.g y): ")
if str(reduce_yn) == 'y':
	# reduce the data
	print "Reducing Data..."
	print "Combining Biases..."
	time.sleep(1)
	zerocmbd=Zerocombine()
	if zerocmbd !=0:
		print "Problem zerocombining, exiting!"
		exit()

	print "Processing Flats..."
	time.sleep(1)	
	flatsprocd=ProcFlats()
	if flatsprocd !=0:
		print "Problem processing flats, exiting!"
		exit()
	
	print "Combining Flats..."
	time.sleep(1)
	flatscmbd=Flatcombine()
	if flatscmbd !=0:
		print "Problem flatcombining, exiting!"
		exit()
		
	print "Processing Arcs..."
	time.sleep(1)
	arcsprocd=ProcArcs()
	if arcsprocd !=0:
		print "Problem processing arcs, exiting!"
		exit()
	
	print "Processing Spectra..."
	time.sleep(1)
	specprocd=ProcSpec()
	if specprocd !=0:
		print "Problem processing spectra, exiting!"
		exit()

	# clean up the working directory
	print "Cleaning up calibration frames..."
	time.sleep(1)
	cleancalibs=CleanCalibs()
	if cleancalibs !=0:
		print "Problem cleaning up calibration frames, exiting!"
		exit()

print "Trim the frames?"
trim_yn = raw_input("(e.g y): ")
if str(trim_yn) == 'y':
	# Trim the vignetted parts off the files
	print "Trimming reduced frames..."
	time.sleep(1)
	trimmed=TrimFiles()
	if trimmed !=0:
		print "Problem trimming frames, exiting!"
		exit()

#print "Set AIRMASS and HJD?"
#am_yn = raw_input("(e.g y): ")
#if str(am_yn) == 'y':
#	# Set airmass
#	print "Setting Airmass..."
#	time.sleep(1)
#	airmass_set=SetAirmass()
#	if airmass_set != 0:
#		print "Problem setting airmass, exiting!"
#		exit()
#		
#	# Set JD
#	print "Setting JD..."
#	time.sleep(1)	
#	jd_set=SetJD()
#	if jd_set !=0:
#		print "Problem setting JD, exiting!"
#		exit()
	
print "Analyse the specta interactively?"
analyse_yn = raw_input("(e.g y): ")
if str(analyse_yn) == 'y':
	# begin interactive analysis
	print "Beginning Interactive Analysis..."
	time.sleep(1)
	analysed=Analyse()
	if analysed != 0:
		print "Problem during spectral analysis, exiting!"
		exit()

print "Round up reduced spectra?"
ru_yn = raw_input("(e.g y): ")
if str(ru_yn) == 'y':
	# round up reduced spectra into one folder
	# make a tar ball for emailing etc
	print "Putting reduced spectra in 'spectra_r'..."
	time.sleep(1)
	roundedup=RoundUpFiles()
	if roundedup != 0:
		print "Problem rounding up the reduced spectra, exiting!"
		exit()

print "Round up arcs?"
rua_yn = raw_input("(e.g y): ")
if str(rua_yn) == 'y':
	# round up files lift by the reduction into one folder
	print "Putting reduction left overs in 'unnormalized'..."
	time.sleep(1)
	roundedupa=RoundUpArcs()
	if roundedupa != 0:
		print "Problem rounding up the reduction left overs, exiting!"
		exit()

endtime=time.ctime()
print "\n\nSpectral reduction & analysis with SPECTOR.PY complete! "
print '\nProgram started at '+starttime
print 'Program ended at '+endtime


