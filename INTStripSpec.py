# Code to strip multi extension INT files and rename them according to
# what type of file they are

#! /usr/bin/env python

from numpy import * 
from numarray import *                       #modules imported
import pyfits, os, sys, string, re, math
import time
import commands
import sys
from pylab import *
import matplotlib.pyplot as plt
from pyraf import iraf

#Import the good packages from IRAF
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)	

print "Strip images?"
strip_yn = raw_input("(e.g y): ")
if str(strip_yn)=='y':

	templist=commands.getoutput('ls r*').split('\n')
	
	for i in range(0,len(templist)):
		image=str(templist[i]+"[1]")
		header=str(templist[i])
		hdulist=pyfits.open(header)
		prihdr=hdulist[0].header
		object=prihdr['IMAGETYP']
		
		if str(object)=='object':
			image2=str("i_s_"+str(templist[i]))
	
		if str(object)=='zero':
			image2=str("b_s_"+str(templist[i]))
			
		if str(object)=='flat':
			image2=str("f_s_"+str(templist[i]))
		
		if str(object)=='arc':
			image2=str("a_s_"+str(templist[i]))
	
		iraf.imcopy(input=image,output=image2)
	
		print image
		print image2

print "Trim images?"
trim_yn = raw_input("(e.g y): ")
if str(trim_yn)=='y':
	
	section='[1:251,1000:3000]'
	
	templist=commands.getoutput('ls i_s_r*').split('\n')
	
	for i in range(0,len(templist)):
		image=str(templist[i])+section
		image2=str(templist[i].split('.')[0]) + "_t"
	
		iraf.imcopy(input=image,output=image2)


	
	