# new version of INT HJD correction

import sys
import numpy as np
from math import *
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
import glob as g

output = ' '
n=965
j=1
pi = 3.141592654
el = np.zeros(55,float)
Rads = (pi/180)
#nm = np.zeros(10,string)

def GetPosition(Year, Month, day, Hour, Minutes, Seconds, Planet):
	#Compute Positions of All Planets
	y = Year
	m = Month
	zDay = day
	h = Hour
	mins = Minutes
	secs = Seconds
	h = h + mins / 60 +  secs / 3600
	p = Planet

	D=  367 * y - floor(7 * (y + floor((m + 9) / 12)) / 4) + floor(275 * m / 9) + zDay - 730531.5 + h / 24
	# the Planets
	#nm[1] = "Mercury"
	#nm[2] = "Venus"
	#nm[3] = "Sun"
	#nm[4] = "Mars"
	#nm[5] = "Jupiter"
	#nm[6] = "Saturn"
	#nm[7] = "Uranus"
	#nm[8] = "Neptune"
	#nm[9] = "Pluto"
	# Mercury
	el[1] = (7.00487 - 0.000000178797 * D) * Rads
	el[2] = (48.33167 - 0.0000033942 * D) * Rads
	el[3] = (77.45645 + 0.00000436208 * D) * Rads
	el[4] = 0.38709893 + 1.80698E-11 * D
	el[5] = 0.20563069 + 0.000000000691855 * D
	el[6] = (Rads * (252.25084 + 4.092338796 * D))
	# Venus
	el[7] = (3.39471 - 0.0000000217507 * D) * Rads
	el[8] = (76.68069 - 0.0000075815 * D) * Rads
	el[9] = (131.53298 - 0.000000827439 * D) * Rads
	el[10] = 0.72333199 + 2.51882E-11 * D
	el[11] = 0.00677323 - 0.00000000135195 * D
	el[12] = (Rads * (181.97973 + 1.602130474 * D))
	# Earth
	el[13] = (0.00005 - 0.000000356985 * D) * Rads
	el[14] = (-11.26064 - 0.00013863 * D) * Rads
	el[15] = (102.94719 + 0.00000911309 * D) * Rads
	el[16] = 1.00000011 - 1.36893E-12 * D
	el[17] = 0.01671022 - 0.00000000104148 * D
	el[18] = (Rads * (100.46435 + 0.985609101 * D))
	# Mars
	el[19] = (1.85061 - 0.000000193703 * D) * Rads
	el[20] = (49.57854 - 0.0000077587 * D) * Rads
	el[21] = (336.04084 + 0.00001187 * D) * Rads
	el[22] = 1.52366231 - 0.000000001977 * D
	el[23] = 0.09341233 - 0.00000000325859 * D
	el[24] = (Rads * (355.45332 + 0.524033035 * D))
	# Jupiter
	el[25] = (1.3053 - 0.0000000315613 * D) * Rads
	el[26] = (100.55615 + 0.00000925675 * D) * Rads
	el[27] = (14.75385 + 0.00000638779 * D) * Rads
	el[28] = 5.20336301 + 0.0000000166289 * D
	el[29] = 0.04839266 - 0.00000000352635 * D
	el[30] = (Rads * (34.40438 + 0.083086762 * D))
	# Saturn
	el[31] = (2.48446 + 0.0000000464674 * D) * Rads
	el[32] = (113.71504 - 0.0000121 * D) * Rads
	el[33] = (92.43194 - 0.0000148216 * D) * Rads
	el[34] = 9.53707032 - 0.0000000825544 * D
	el[35] = 0.0541506 - 0.0000000100649 * D
	el[36] = (Rads * (49.94432 + 0.033470629 * D))
	# Uranus
	el[37] = (0.76986 - 0.0000000158947 * D) * Rads
	el[38] = (74.22988 + 0.0000127873 * D) * Rads
	el[39] = (170.96424 + 0.0000099822 * D) * Rads
	el[40] = 19.19126393 + 0.0000000416222 * D
	el[41] = 0.04716771 - 0.00000000524298 * D
	el[42] = (Rads * (313.23218 + 0.011731294 * D))
	# Neptune
	el[43] = (1.76917 - 0.0000000276827 * D) * Rads
	el[44] = (131.72169 - 0.0000011503 * D) * Rads
	el[45] = (44.97135 - 0.00000642201 * D) * Rads
	el[46] = 30.06896348 - 0.0000000342768 * D
	el[47] = 0.00858587 + 0.000000000688296 * D
	el[48] = (Rads * (304.88003 + 0.0059810572 * D))
	# Pluto
	el[49] = (17.14175 + 0.0000000841889 * D) * Rads
	el[50] = (110.30347 - 0.0000002839 * D) * Rads
	el[51] = (224.06676 - 0.00000100578 * D) * Rads
	el[52] = 39.48168677 - 0.0000000210574 * D
	el[53] = 0.24880766 + 0.00000000177002 * D
	el[54] = (Rads * (238.92881 + 0.003931834 * D))
	
	q = 6 * (p - 1)
	ip = el[q + 1]
	op = el[q + 2]
	pp = el[q + 3]
	ap = el[q + 4]
	ep = el[q + 5]
	lp = el[q + 6]
	ie = el[13]
	oe = el[14]
	pe = el[15]
	ae = el[16]
	ee = el[17]
	le = el[18]
	
	# Get Earth's position using Kepler's equation
	Me1 = (le - pe)
	B = Me1 / (2 * pi)
	Me1 = 2 * pi * (B - floor(abs(B)))
	
	if (B<0): 
		Me1 = 2 * pi * (B + floor(abs(B)))
	if (Me1 < 0): 
		Me1 = 2 * pi + Me1
	
	e = Me1
	delta = 0.05
	
	while (abs(delta) >= pow(10,-12)):
 		delta = e - ee * sin(e) - Me1
 		e = e - delta / (1 - ee * cos(e))
 	
 	
	ve = 2 * atan(pow(((1 + ee) / (1 - ee)) , 0.5) * tan(0.5 * e))
	
	if (ve < 0):
		ve = ve + 2 * pi 
	
	re = ae * (1 - ee * ee) / (1 + ee * cos(ve))
	xe = re * cos(ve + pe)
	ye = re * sin(ve + pe)
	ze = 0
	
	
	# Get planet's position using Kepler's equation
	mp = (lp - pp)
	B = mp / (2 * pi)
	mp = 2 * pi * (B - floor(abs(B)))
	
	if (B<0): 
		mp = 2 * pi * (B + floor(abs(B)))
	
	if (mp < 0):
		mp = 2 * pi + mp
	
	e = mp
	delta = 0.05
	
	while (abs(delta) >= pow(10,-12)):
		delta = e - ep * sin(e) - mp
		e = e - delta / (1 - ep * cos(e))
		
	vp = 2 * atan(pow(((1 + ep) / (1 - ep)) , 0.5) * tan(0.5 * e))
 
	if (vp < 0):
		vp = vp + 2 * pi
		
	rp = ap * (1 - ep * ep) / (1 + ep * cos(vp))
	xh = rp * (cos(op) * cos(vp + pp - op) - sin(op) * sin(vp + pp - op) * cos(ip))
	yh = rp * (sin(op) * cos(vp + pp - op) + cos(op) * sin(vp + pp - op) * cos(ip))
	zh = rp * (sin(vp + pp - op) * sin(ip))
	xg = xh - xe
	yg = yh - ye
	zg = zh
	
	# compute RA and DEC
	ecl = 23.439292 * Rads
	xeq = xg
	yeq = yg * cos(ecl) - zg * sin(ecl)
	zeq = yg * sin(ecl) + zg * cos(ecl)
	ra = atan(yeq/ xeq)*12/pi
	
	if (xeq < 0):
		ra = ra + 12

	if (yeq < 0):
		if (xeq>0): 
			ra = ra + 24
	
	dec = 180*atan(zeq / pow((xeq * xeq + yeq * yeq),0.5))/pi
	
	# Sun Coodinates
	xeq = xe
	yeq = ye * cos(ecl) - ze * sin(ecl)
	zeq = ye * sin(ecl) + ze * cos(ecl)
	rae = 12 + atan(yeq/ xeq)*12/pi
	
	if (xe < 0):
		rae = rae + 12

	if (ye < 0):
		if (xe>0):
			rae = rae + 24
  
	dece = -180*atan(zeq / pow((xeq * xeq + yeq * yeq),0.5))/pi
  
	if (p==3):
		ra=rae
		dec=dece

	if (ra<12):
		raa=12-ra
	else:
		raa=36-ra

	xpix=32+floor(raa*965/24)
	ypix=201+floor(dec*-3)

	return ra,dec,D

def gethjd(YEAR,MONTH,DAY,HOUR,MINUTE,SEC,ORAH,ORAM,ORAS,ODECD,ODECM,ODECS):
	# Sun
	MM=MONTH     
	DD=DAY     
	YY=YEAR  
	HR=HOUR    
	MN=MINUTE    
	SC=SEC    

	ra,dec,D=GetPosition(YY,MM,DD,HR,MN,SC,3)

	rah=floor(ra)
	ram=floor((ra - floor(ra))*60)
	decd=floor(abs(dec))

	if (dec<0):
		decd=-1*decd
	
	decm=floor((abs(dec)-floor(abs(dec)))*60)
	ASDF=str(MM) + "/" + str(DD) + "/" + str(YY) + "   " + str(HR) + ":" + str(MN) + ":" + str(SC)
	ASDF=str(ASDF) + "\n" + "Sun -- RA: "+str(rah)+"h "+str(ram)+"m   DEC: "+str(decd)+" "+str(decm)+"'"

	# Earth
	dec=-1*dec
	ra=ra+12

	if (ra>24):
		ra=ra-24

	rah=floor(ra)
	ram=floor((ra -floor(ra))*60)
	decd=floor(abs(dec))

	if (dec<0):
		decd=-1*decd
	
	decm=floor((abs(dec)-floor(abs(dec)))*60)
	ASDF=str(ASDF) + "\n" + "Earth -- RA: "+str(rah)+"h "+str(ram)+"m   DEC: "+str(decd)+" "+str(decm)+"'"

	ORA = (ORAH + ORAM/60 + ORAS/3600)*15
	ODEC = abs(ODECD) + ODECM/60 + ODECS/3600

	if (ODECD<0):
		ODEC=-1*ODEC

	ASDF=str(ASDF) + "\n" + "Object -- RA: " + str(ORA) + "   DEC: " + str(ODEC) + " "

	# Earth XYZ
	cel = cos(dec * pi/180)
	earthx = cos(ra * pi/12) * cel
	earthy = sin(ra * pi/12) * cel
	earthz = sin(dec * pi/180)
	ASDF=str(ASDF) + "\n" + "Earth -- X,Y,Z: " + str(earthx) + "," + str(earthy) + "," + str(earthz)

	# Object XYZ
	cel = cos(ODEC * pi/180)
	objectx = cos(ORA * pi/180) * cel
	objecty = sin(ORA * pi/180) * cel
	objectz = sin(ODEC * pi/180)
	ASDF=str(ASDF) + "\n" + "Object -- X,Y,Z: " + str(objectx) + "," + str(objecty) + "," + str(objectz)

	# Light Time (Minutes per AU)
	ausec=8.3168775
	correction = ausec * (earthx * objectx + earthy * objecty + earthz * objectz)
	ASDF=str(ASDF) + "\n" + "Light Time: " + str(correction) + " minutes"
	D=D+2451545.
	ASDF=str(ASDF) + "\nJulian Date: " + str(D)
	D1=D+correction/(24*60)
	ASDF=str(ASDF) + "\nHeliocentric Julian Date: " + str(float(D1))

	#form.coordinates.value=ASDF
	print ASDF
	
	#	return D
	return D1,D

def ImportPackages():
	iraf.noao(_doprint=0)
	return 0	
	
t=g.glob('*.fits')
start_times,mid_times,coords=[],[],[]
HJD=np.empty(len(t))
JD_e=np.empty(len(t))

# load IRAF for updating the headers
from pyraf import iraf
imp=ImportPackages()
if imp != 0:
	print "Problem loading IRAF packages, exiting!"
	sys.exit(1)

for i in range(0,len(t)):
	hdr=fits.open(t[i])[0].header
	d=hdr['DATE-OBS']
	uts=hdr['UTSTART']
	ra=hdr['CAT-RA']
	dec=hdr['CAT-DEC']
	object=hdr['OBJECT']
	exptime=float(hdr['EXPTIME'])
	start_times.append(Time('%s %s' %(d,uts),format='iso', scale='utc'))
	mid_times.append(start_times[i]+(exptime/2.)*u.second)
	coords.append(SkyCoord('%s %s' % (ra,dec),unit=(u.hourangle,u.deg)))
		
	# call this xt to shorten next line
	xt=mid_times[i].datetime
	c_ra=coords[i].ra
	c_dec=coords[i].dec
	HJD[i],JD_e[i]=gethjd(float(xt.year),float(xt.month),float(xt.day),float(xt.hour),float(xt.minute),float(xt.second)+float(xt.microsecond/1.0E6),float(c_ra.hms.h),float(c_ra.hms.m),float(c_ra.hms.s),float(c_dec.dms.d),float(c_dec.dms.m),float(c_dec.dms.s))
	print "%s %s %s %.2f %s %s %.8f" % (object,ra,dec,exptime,start_times[i],mid_times[i],HJD[i])
	
	# ok so now we have the real JD-MID, HJD-MID and UTMIDDLE
	# update the headers of the spectra
	iraf.hedit(images=t[i],fields='HJD',value=str('%.8f' % HJD[i]),add='yes',verify='no',show='yes')
	iraf.hedit(images=t[i],fields='JD_e',value='xxx',add='yes',verify='no',show='yes')
	iraf.hedit(images=t[i],fields='JD_e',value=str('%.8f' % JD_e[i]),add='yes',verify='no',show='yes')
	iraf.hedit(image=t[i],field='UTMIDDLE',value="%s" % (mid_times[i].iso.replace(' ','T')),add='yes',ver='no',show='yes')

# ensure the airmass is ok now that utmiddle etc are fixed	
iraf.setairmass(images="*.fits", ra="cat-ra", dec="cat-dec", equinox="cat-epoc", st="st", ut="utstart", date="date", exposur="exptime", airmass="airmass", utmiddl="utmiddle", show="yes", update="yes", overrid="yes",mode="ql")	

