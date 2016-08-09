"""
########################################################################################
################################## SPECTOR.py ##########################################
########################################################################################
#
#   INT/IDS Spectral Reduction v2.0: A Script for complete reduction of INT/IDS data
#
#                               James McCormac
#
#   Before running Spector:
#       -   Remove all files that should not be included in reduction
#           i.e. bad flats, data from other central wavelengths or gratings etc.
#           Place all this data in a folder called 'junk'. Use *.int logfile to
#           help weed out unwanted files
#       -   Adjust the settings in the **EDIT HERE** section below
#
#   Assumptions:
#       -   This code was written for a particular observing strategy:
#               Acquire, Arc, Science, Arc, REPEAT...
#           Therefore the script is looking for this structure in the data when matching 
#           arcs to spectra. This can be easily changed by modifying the relevant section
#
#   What Spector does:
#       -   Prepares the working environment, making copies of original data
#       -   Loads IRAF
#       -   Strips multi extension fits to single extension
#       -   Renames files for easier viewing:
#               i_* = integration (spectra)
#               a_* = arc
#               b_* = bias
#               f_* = flat
#       -   Debiases and flat field corrects the data
#       -   Tidies up calibration frames into 'calibs' folder
#       -   Trims the spectra based on INT's heavily vignetted region, see TRIMSEC below
#       -   Interactive spectral analysis with KPNOSLIT in IRAF
#       -   Normalises the spectra
#       -   Tidies up working directory after reductions are finished
#
#   What Spector does not do:
#       -   Flux calibration
#
#   What you will end up with:
#       -   A set of reduced, wavelength calibrated, normalised 1D spectra with the
#           byproducts from each step located in subdirectories within the 'reduced'
#           folder.
#
#   Who to annoy if you have questions:
#       jmccormac001@gmail.com
#
# Changes:
#   2015-10-06: Removed setdJD and setAirmass as they were deeply broken!
#
"""
import sys
import os
import re
import time
import argparse as ap
from datetime import (
    datetime,
    timedelta
    )
from collections import defaultdict
import numpy as np
from ccdproc import (
    CCDData,
    ImageFileCollection,
    combine,
    subtract_bias,
    subtract_dark,
    flat_correct,
    trim_image
    )
from astropy.io import fits
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    SkyCoord
    )
import matplotlib.pyplot as plt
import glob as g
#import pyds9

# set up some exceptions to work cross Python2 and 3
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

#################################
########### EDIT HERE ###########
#################################

# file locations
lineList_location = "/Users/James/Documents/LineLists/cuar_cune.dat"
loginCl_location = "/Users/James/"
# CCD params
GAIN = "1.03"
RDNOISE = "3.48"
SATURATION = "40000"
# header keyword map
BIAS_KEYWORD = 'zero'
FLAT_KEYWORD = 'flat'
ARC_KEYWORD = 'arc'
SCIENCE_KEYWORD = 'object'
EXPTIME_KEYWORD = 'EXPTIME'
UTSTART_KEYWORD = 'UTSTART'
DATEOBS_KEYWORD = 'DATE-OBS'
RA_KEYWORD = 'CAT-RA'
DEC_KEYWORD = 'CAT-DEC'

# Observatory
# this fails when there is no internet, how stupid!
# do it manually so I can work on the bus
#OBSERVATORY = EarthLocation.of_site('Roque de los Muchachos')
olat = 28.+(45./60.)-(37./3600.)
olon = -17.-(52./60.)-(46./3600.)
elev = 2332.
OBSERVATORY = EarthLocation(lat=olat*u.deg,lon=olon*u.deg,height=elev*u.m)

#################################

def makeDirs():
    """
    Make directories for original and reduced data
    """
    if os.path.isdir('original') == False:
        os.mkdir('original')
    if os.path.isdir('reduced') == False:
        os.mkdir('reduced')

def copyFiles():
    """
    Copy the data
    """
    currentdir = os.getcwd()
    contents = os.listdir(currentdir)
    if not os.path.isdir('original'):
        os.mkdir('original')
    else:
        print('original directory exists, it shouldn\'t, exiting...')
        sys.exit(1)
    for i in contents:
        os.system('mv {} original/'.format(i))
        print("Moved {} to original folder".format(i))
    if not os.path.isdir('reduced'):
        os.mkdir('reduced')
    else:
        print('reduced directory exists, it shouldn\'t, exiting...')
        sys.exit(1)
    os.chdir('original/')
    for i in g.glob('r*.fit*'):
        os.system("cp {} ../reduced/".format(i))
        print("Copied {} to reduction folder".format(i))
    os.chdir('../reduced/')

def getImageList():
    """
    Return ImageFileCollection object for the current
    working directory
    """
    return ImageFileCollection('.')

def renameFiles(images):
    """
    Rename the files so they are easier to scan by eye
    """
    for f in images.files_filtered(imagetyp=FLAT_KEYWORD):
        os.rename(f, "f_s_{}".format(f))
    for f in images.files_filtered(imagetyp=BIAS_KEYWORD):
        os.rename(f, "b_s_{}".format(f))
    for f in images.files_filtered(imagetyp=ARC_KEYWORD):
        os.rename(f, "a_s_{}".format(f))
    for f in images.files_filtered(imagetyp=SCIENCE_KEYWORD):
        os.rename(f, "i_s_{}".format(f))

def makeMasterBias(images):
    """
    Make a master bias using all biases found in
    images object

    TODO: Finish docstring
    """
    try:
        master_bias = CCDData.read('master_bias.fits', unit=u.adu)
        return master_bias
    except FileNotFoundError:
        bias_list = []
        for f in images.files_filtered(imagetyp=BIAS_KEYWORD):
            print(f)
            ccd = CCDData.read(f, unit=u.adu)
            bias_list.append(ccd)
        try:
            master_bias = combine(bias_list, method='median')
            master_bias.write('master_bias.fits', clobber=True)
            return master_bias
        except IndexError:
            return  None

def makeMasterFlat(images, master_bias):
    """
    Flats are corrected for their bias level (if
    master_bias)

    TODO: Finish docstring
    """
    try:
        fitsfile = 'master_flat.fits'
        master_flat = CCDData.read(fitsfile, unit=u.adu)
        return master_flat
    except FileNotFoundError:
        # empty list for the flats
        flat_list = []
        # create the master flat field
        print('Reducing flats')
        for f in images.files_filtered(imagetyp=FLAT_KEYWORD):
            print(f)
            with fits.open(f) as fitsfile:
                data_exp = fitsfile[0].header[EXPTIME_KEYWORD]
            ccd = CCDData.read(f, unit=u.adu)
            if master_bias:
                ccd = subtract_bias(ccd, master_bias)
            else:
                print('No master bias, skipping correction...')
            flat_list.append(ccd)
        try:
            master_flat = combine(flat_list, method='median')
            master_flat.write('master_flat.fits', clobber=True)
            return master_flat
        except IndexError:
            print('There are no flats, skipping...')
            master_flat = None

def getLightTravelTimes(ra, dec, time2corr):
    """
    Get the light travel times to the helio and
    bary centres

    TODO: Finish docstring
    """
    target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    ltt_bary = time2corr.light_travel_time(target)
    ltt_helio = time2corr.light_travel_time(target, 'heliocentric')
    return ltt_bary, ltt_helio

def correctData(filename, master_bias, master_flat, filetype):
    """
    Correct a science image using the available
    master calibrations. Skip a calibration step if the
    master frame does not exist.

    No reduced file is written in this new scheme.
    Instead, the corrected data is passed directly
    to the phot() routine, photometry is done as per
    the configuration and the photometry is written out
    only.

    TODO: Finish docstring
    """
    print('Reducing {0:s}...'.format(filename))
    with fits.open(filename) as fitsfile:
        # correct times for science spectra,
        # don't bother for arcs
        hdr = fitsfile[0].header
        if filetype == 'science':
            half_exptime = hdr[EXPTIME_KEYWORD]/2.
            utstart = hdr[UTSTART_KEYWORD]
            dateobs = hdr[DATEOBS_KEYWORD]
            ra = hdr[RA_KEYWORD]
            dec = hdr[DEC_KEYWORD]
            time_start = Time('{}T{}'.format(dateobs, utstart),
                                             scale='utc',
                                             format='isot',
                                             location=OBSERVATORY)
            # correct to mid exposure time
            jd_mid = time_start + half_exptime*u.second
            ltt_bary, ltt_helio = getLightTravelTimes(ra, dec, jd_mid)
            time_bary = jd_mid.tdb + ltt_bary
            time_helio = jd_mid.utc + ltt_helio
            hdr['BJD-MID'] = time_bary.jd
            hdr['HJD-MID'] = time_helio.jd
            hdr['JD-MID'] = jd_mid.jd
    ccd = CCDData.read(filename, unit=u.adu)
    if master_bias:
        ccd = subtract_bias(ccd, master_bias)
    else:
        print('No master bias, skipping correction...')
    if master_flat:
        ccd = flat_correct(ccd, master_flat)
    else:
        print('No master flat, skipping correction...')
    # after calibrating we get np.float64 data
    # if there are no calibrations we maintain dtype = np.uint16
    # sep weeps
    # fix this by doing the following
    if isinstance(ccd.data[0][0], np.uint16):
        ccd.data = ccd.data.astype(np.float64)
    # trim the data
    ccd_trimmed = trim_image(ccd[1000:3001,:])
    # write out the trimmed file and the updated header
    #ccd_trimmed.write(filename, hdr, clobber=True)
    trimmed_filename = '{}_t.fits'.format(filename.split('.')[0])
    fits.writeto(trimmed_filename, ccd_trimmed.data, hdr)
    # remove the old untrimmed data
    os.system('rm {}'.format(filename))

def cleanCalibs():
    """
    Tidy up the calibration frames
    """
    os.mkdir('calibs')
    for i in g.glob('f_s_r*'):
        os.system('mv {} calibs/'.format(i))
    for i in g.glob('b_s_r*'):
        os.system('mv {} calibs/'.format(i))
    os.system('mv master_bias.fits calibs/')
    os.system('mv master_flat.fits calibs/')

def extractSpectra():
    """
    Extract 1D spectra using IRAF interactively

    Interpolate across the two arcs either side to
    get the most accurate wavelength solution

    TODO: Finish docstring
          Add method of using super arc for inital 
          identify
    """
    # load IRAF from the location of the login.cl file
    here=os.getcwd()
    os.chdir(loginCl_location)
    from pyraf import iraf
    os.chdir(here)
    time.sleep(2)

    # make a list of the science images to be analysed
    templist = g.glob('i_s*')
    # import IRAF packages for spectroscopy
    iraf.imred(_doprint=0)
    iraf.kpnoslit(_doprint=0)
    # apall parameters
    iraf.apall.setParam('format', 'multispec')
    iraf.apall.setParam('interac', 'yes')
    iraf.apall.setParam('find', 'yes')
    iraf.apall.setParam('recen', 'yes')
    iraf.apall.setParam('resize', 'yes')
    iraf.apall.setParam('trace', 'yes')
    iraf.apall.setParam('fittrac', 'yes')
    iraf.apall.setParam('extract', 'yes')
    iraf.apall.setParam('extras', 'yes')
    iraf.apall.setParam('review', 'yes')
    iraf.apall.setParam('line', 'INDEF')
    iraf.apall.setParam('nsum', '12')
    iraf.apall.setParam('lower', '-6')
    iraf.apall.setParam('upper', '6')
    iraf.apall.setParam('b_funct', 'chebyshev')
    iraf.apall.setParam('b_order', '1')
    iraf.apall.setParam('b_sampl', '-25:-15,15:25')
    iraf.apall.setParam('b_naver', '-100')
    iraf.apall.setParam('b_niter', '0')
    iraf.apall.setParam('b_low_r', '3')
    iraf.apall.setParam('b_high', '3')
    iraf.apall.setParam('b_grow', '0')
    iraf.apall.setParam('width', '10')
    iraf.apall.setParam('radius', '10')
    iraf.apall.setParam('threshold', '0')
    iraf.apall.setParam('nfind', '1')
    iraf.apall.setParam('t_nsum', '10')
    iraf.apall.setParam('t_step', '10')
    iraf.apall.setParam('t_nlost', '3')
    iraf.apall.setParam('t_funct', 'spline3')
    iraf.apall.setParam('t_order', '3')
    iraf.apall.setParam('backgro', 'fit')
    iraf.apall.setParam('skybox', '1')
    iraf.apall.setParam('weights', 'variance')
    iraf.apall.setParam('pfit', 'fit1d')
    iraf.apall.setParam('clean', 'yes')
    iraf.apall.setParam('saturat', SATURATION)
    iraf.apall.setParam('readnoi', RDNOISE)
    iraf.apall.setParam('gain', GAIN)
    iraf.apall.setParam('lsigma', '4.0')
    iraf.apall.setParam('usigma', '4.0')
    iraf.apall.setParam('nsubaps', '1')
    iraf.apall.saveParList(filename="apall.pars")

    # TODO
    # Make a reference using a long exposure arc
    # TODO

    # make reference arc for reidentify
    refarc=None
    for i in range(0,len(templist)):
        hdulist = fits.open(templist[i])
        prihdr = hdulist[0].header
        target_id = prihdr['CAT-NAME']
        # extract the object spectrum
        print("Extracting spectrum of {} from image {}".format(target_id, templist[i]))
        print("Check aperture and background. Change if required")
        print("AP: m = mark aperture, d = delete aperture")
        print("SKY: s = mark sky, t = delete sky, f = refit")
        print("q = continue")
        iraf.apall(input=image)
        print("Spectrum extracted!")
        # find the arcs either side of the object
        arc1 = "a_s_r{0:d}_t.fits".format(int(templist[i].split('_')[2].split('r')[1])-1)
        arc2 = "a_s_r{0:d}_t.fits".format(int(templist[i].split('_')[2].split('r')[1])+1)
        # predict the arc names
        print("\nPredicting arcs names...")
        print("Arc1: {}".format(arc1))
        print("Arc2: {}".format(arc2))
        # setup a reference filename for the arc conditions in database
        reffile = templist[i].split('.fits')[0]
        # extract the arcs
        print("\nExtracting arcs under the same conditions...")
        iraf.apall(input=arc1,
                   reference=reffile,
                   recente="no",
                   trace="no",
                   backgro="no",
                   interac="no")
        print("Arc1 extracted")
        iraf.apall(input=arc2,
                   reference=reffile,
                   recente="no",
                   trace="no",
                   backgro="no",
                   interac="no")
        print("Arc2 extracted")
        # get a list of the extracted arcs and objects
        templist2 = g.glob('a_s*.ms.fits')[-2:]
        templist3 = g.glob('i_s*.ms.fits')[-1]
        print("\nArcs extracted to:")
        print(templist2[0])
        print(templist2[1])
        if i==0:
            print("\nIdentify arc lines:")
            print("Enter the following in the splot window")
            print("\t:thres 500")
            print("\t:order 3")
            print("\tfwidth 2")
            print("Select 3-5 arc lines from line atlas")
            print("Press 'm' to mark, then enter wavelength")
            print("Then press 'l' to automatically ID the other lines")
            print("Press 'f' to fit the dispersion correction")
            print("Use 'd' to remove bad points, 'f' to refit")
            print("'q' from fit, then 'q' from identify to continue\n")
            refarc = templist2[0]
            iraf.identify(images=refarc, coordlist=lineList_location)
            iraf.reidentify(reference=refarc, images=templist2[1])
        if i>0:
            print ("\nReidentifying arclines from {}".format(refarc))
            iraf.reidentify(reference=refarc, images=templist2[0])
            iraf.reidentify(reference=refarc, images=templist2[1])
        # add the refspec keywords to the image header for dispcor
        refspec1 = "{} 0.5".format(templist2[0].split(".fits")[0])
        refspec2 = "{} 0.5".format(templist2[1].split(".fits")[0])
        print("\nREFSPEC1: {}".format(refspec1))
        print("REFSPEC2: {}".format(refspec2))
        iraf.hedit(images=templist3,
                   fields="REFSPEC1",
                   value=refspec1,
                   add="yes",
                   verify="no",
                   show="yes")
        iraf.hedit(images=templist3,
                   fields="REFSPEC2",
                   value=refspec2,
                   add="yes",
                   verify="no",
                   show="yes")
        print("Headers updated!\n")
        # apply the dispersion correction
        print("Applying the dispersion correction")
        iraf.dispcor(input=templist3,
                     output=templist3,
                     lineari="yes",
                     databas="database",
                     table="")
        print("Correction applied!")
        # normalize the spectrum using continuum
        normspec = "{}n.ms.fits".format(templist3.split('.ms')[0])
        iraf.continuum(input=templist3,
                       output=normspec,
                       logfile="logfile",
                       interac="yes",
                       functio="spline3",
                       order="5",
                       niterat="10",
                       markrej="yes")
        print("\n\n")

def roundUpSpectra():
    """
    Gather up reduced spectra
    """
    if os.path.exists('spectra_r') == False:
        os.mkdir('spectra_r')
    for i in g.glob('*_tn.ms.fits'):
        os.system('mv {} spectra_r/'.format(i))

def roundUpArcs():
    """
    Gather up all the reduced arcs and put them 
    somewhere out of the way
    """
    templist=g.glob('i*t.ms.fits')
    templist2=g.glob('i*t.fits')
    templist3=g.glob('a*t.ms.fits')
    templist4=g.glob('a*t.fits')
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

if __name__ == '__main__':
    print '\n\n---------------------------------------------------------'
    print '------------- Spectroscopy by Spector.py ----------------'
    print '---------------------------------------------------------\n'

    junk_yn = raw_input("Has the junk been removed from current folder? (y/n): ")
    if junk_yn != 'y':
        print('Remove all files that should not be analysed, then restart')
        sys.exit(1)

    # make a local backup copy of the data
    copyFiles()
    # get a list of images in current directory
    images = getImageList()
    # rename the images to make them easier to read
    renameFiles(images)
    # get new filenames
    images = getImageList()
    # make a master bias
    master_bias = makeMasterBias(images)
    # make a master flat
    master_flat = makeMasterFlat(images, master_bias)
    # correct the arcs
    for filename in images.files_filtered(imagetyp=ARC_KEYWORD):
        correctData(filename, master_bias, master_flat, 'arc')
    # correct the science spectra
    for filename in images.files_filtered(imagetyp=SCIENCE_KEYWORD):
        correctData(filename, master_bias, master_flat, 'science')
    cleanCalibs()
    # extract wavelength calibrated and continuum normalised spectra
    extractSpectra()
    # round up the reduced spectra
    roundUpSpectra()

    cont = 0
    if cont > 0:

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

