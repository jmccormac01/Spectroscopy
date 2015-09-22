from pyraf import iraf
import commands
import pyfits

def ImportPackages():
	iraf.noao(_doprint=0)
	iraf.onedspec(_doprint=0)

	return 0

imported=ImportPackages()
	
if imported != 0:
	print "Problem importing IRAF Packages, exiting!"
	exit()

templist=commands.getoutput('ls *_tn.ms.fits').split('\n')

print "How many spectra are to be combined?"
num=int(raw_input('(e.g 4): '))


x=templist[0:num]
y=templist[num:len(templist)+1]

if num == 2:
	xin=str(x[0]) + ', ' + str(x[1])
	yin=str(y[0]) + ', ' + str(y[1])

if num == 3:
	xin=str(x[0]) + ', ' + str(x[1]) + ', ' + str(x[2])
	yin=str(y[0]) + ', ' + str(y[1]) + ', ' + str(y[2])
	
if num == 4:
	xin=str(x[0]) + ', ' + str(x[1]) + ', ' + str(x[2]) + ', ' + str(x[3])
	yin=str(y[0]) + ', ' + str(y[1]) + ', ' + str(y[2])	+ ', ' + str(y[3])


xout=str(x[0].split('.')[0]) + 'C.' + str(x[0].split('.')[1]) + '.' + str(x[0].split('.')[2])
yout=str(y[0].split('.')[0]) + 'C.' + str(y[0].split('.')[1]) + '.' + str(y[0].split('.')[2])

print x
print y

print xin
print yin

print xout
print yout

iraf.scombine(input=xin,output=xout)
iraf.scombine(input=yin,output=yout)

