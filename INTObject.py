
# python program to split up all the INT spectra into their own folders

import os,commands
import pyfits
from pylab import *

print "sort?"
sort_yn = raw_input("(e.g. y): ")
if str(sort_yn) == 'y':
	
	templist=commands.getoutput('ls *tn.ms.fits').split('\n')
	
	for i in range(0,len(templist)):
		image=templist[i]		
		
		hdulist=pyfits.open(templist[i])
		prihdr=hdulist[0].header
		object=prihdr['CAT-NAME']
		
		print object
		
		if os.path.exists(object) == False:
			os.mkdir(object)
			
		line=str(object)+"/"+str(templist[i])	
		os.rename(templist[i],line)
		
		print line
	
print "compare?"
compare_yn = raw_input("(e.g. y): ")
if str(compare_yn) == 'y':

	templist=commands.getoutput('ls *t.ms.fits').split('\n')

	ref=templist[0]
	hdulist=pyfits.open(ref)
	ref_data=hdulist[0].data
	ref_spec=ref_data[0][0]
	
	for i in range(1,len(templist)):
		check=templist[i]
		hdulist=pyfits.open(check)
		check_data=hdulist[0].data
		check_spec=check_data[0][0]
		
		dif_spec=ref_spec-check_spec
		
		plot(dif_spec)
		
		show()
		
		
		
		
		
		
		
		
	