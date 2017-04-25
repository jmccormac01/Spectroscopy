"""
Spector - A tool for extracting 1D spectra from INT/IDS

TODO:
    Implement a way to do the blends at the same time
    Fix the flat field normalisation

"""
import sys
import os
import time
import glob as g
import argparse as ap
from datetime import (
    datetime,
    timedelta
    )
import numpy as np
from ccdproc import (
    CCDData,
    ImageFileCollection,
    combine,
    subtract_bias,
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
import pymysql

# set up some exceptions to work cross Python2 and 3
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

# pylint: disable=invalid-name
# pylint: disable=redefined-builtin
# pylint: disable=no-member
# pylint: redefined-outer-name

#################################
########### EDIT HERE ###########
#################################

# file locations
lineList_location = "/Users/jmcc/Dropbox/LineLists/cuar_cune.dat"
loginCl_location = "/Users/jmcc/"
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
olat = 28.+(45./60.)-(37./3600.)
olon = -17.-(52./60.)-(46./3600.)
elev = 2332.
OBSERVATORY = EarthLocation(lat=olat*u.deg, lon=olon*u.deg, height=elev*u.m)
# this below is useful but fails when there is no internet
# OBSERVATORY = EarthLocation.of_site('Roque de los Muchachos')

#################################

def argParse():
    parser = ap.ArgumentParser()
    parser.add_argument('refarc', help='image_id of the long arc used for master solution')
    parser.add_argument('--ready',
                        help='set this to show all the junk is removed',
                        action='store_true')
    parser.add_argument('--log_only',
                        help='use this to log spectra to db only',
                        action='store_true')
    parser.add_argument('--ds9',
                        help='use this to display spectra in DS9',
                        action='store_true')
    return parser.parse_args()

def makeDirs():
    """
    Make directories for original and reduced data
    """
    if not os.path.isdir('original'):
        os.mkdir('original')
    if not os.path.isdir('reduced'):
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
            return None

def estimateSkyLevel(data):
    """
    Function to interatively sigma clip the sky background
    to estimate the sky level without the influence of stars
    """
    mean_diff = 1E6
    mean_diff_limit = 1E-6
    sigma = 3
    # create a masked array where nothing is masked
    data = np.ma.masked_where(data < -1E6, data)
    i = 0
    while mean_diff > mean_diff_limit:
        mean = np.ma.average(data)
        rms = np.ma.std(data)
        masked_data = np.ma.masked_where(((data > mean+sigma*rms) | (data < mean-sigma*rms)), data)
        new_mean = np.ma.average(masked_data)
        new_rms = np.ma.std(masked_data)
        print('Sky level: {}, RMS: {}'.format(new_mean, new_rms))
        data = masked_data
        mean_diff = abs(new_mean-mean)/new_mean
    return new_mean, new_rms

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
            sky_level, sky_rms = estimateSkyLevel(ccd.data)
            ccd.data = ccd.data/sky_level
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
            hdr['UT-MID'] = jd_mid.isot
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
    ccd_trimmed = trim_image(ccd[1000:3001, :])
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
    here = os.getcwd()
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
    iraf.apall.setParam('t_niter', '7')
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
    iraf.identify.setParam('fwidth', '2.')
    iraf.identify.setParam('order', '4')
    iraf.identify.setParam('niterate', '7')
    iraf.apall.saveParList(filename="identify.pars")

    # make reference arc for reidentify
    if '.' in args.refarc:
        args.refarc = args.refarc.split('.')[0]
    refarc = "a_s_{}_t.fits".format(args.refarc)
    refarc_out = "a_s_{}_t.ms.fits".format(args.refarc)

    # loop over all the the spectra
    for i in range(0, len(templist)):
        hdulist = fits.open(templist[i])
        prihdr = hdulist[0].header
        target_id = prihdr['CAT-NAME']
        spectrum_id = int(templist[i].split('_')[2].split('r')[1])
        if args.ds9:
            os.system('xpaset SpectorDS9 fits < {}'.format(templist[i]))
        # extract the object spectrum
        print("[{}/{}] Extracting spectrum of {} from image {}".format(i+1, len(templist), target_id, templist[i]))
        print("[{}/{}] Check aperture and background. Change if required".format(i+1, len(templist)))
        print("[{}/{}] AP: m = mark aperture, d = delete aperture".format(i+1, len(templist)))
        print("[{}/{}] SKY: s = mark sky, t = delete sky, f = refit".format(i+1, len(templist)))
        print("[{}/{}] q = continue".format(i+1, len(templist)))
        iraf.apall(input=templist[i])
        print("Spectrum extracted!")
        # find the arcs either side of the object
        arclist = []
        arc1 = "a_s_r{0:d}_t.fits".format(spectrum_id-1)
        arc2 = "a_s_r{0:d}_t.fits".format(spectrum_id+1)
        arc1_out = "a_s_r{0:d}_t.ms.fits".format(spectrum_id-1)
        arc2_out = "a_s_r{0:d}_t.ms.fits".format(spectrum_id+1)
        # predict the arc names
        print("\nPredicting arcs names...")
        print("Arc1: {}".format(arc1))
        print("Arc2: {}".format(arc2))
        # setup a reference filename for the arc conditions in database
        reffile = templist[i].split('.fits')[0]
        # extract the arcs
        print("\nExtracting arcs under the same conditions...")
        if os.path.exists(arc1):
            iraf.apall(input=arc1,
                       reference=reffile,
                       recente="no",
                       trace="no",
                       backgro="no",
                       interac="no")
            print("Arc1 {} extracted".format(arc1))
            arclist.append(arc1_out)
        else:
            print("\n\nArc1 {} FILE NOT FOUND\n\n".format(arc1))
        if os.path.exists(arc2):
            iraf.apall(input=arc2,
                       reference=reffile,
                       recente="no",
                       trace="no",
                       backgro="no",
                       interac="no")
            print("Arc2 {} extracted".format(arc2))
            arclist.append(arc2_out)
        else:
            print("\n\nArc2 {} FILE NOT FOUND\n\n".format(arc2))
        # get a list of the extracted arcs and objects
        spectrum_out = "i_s_r{0:d}_t.ms.fits".format(spectrum_id)
        if i == 0:
            # extract the master reference arc
            print("\nExtracting master arc {} under the same conditions...".format(refarc))
            iraf.apall(input=refarc,
                       reference=reffile,
                       recente="no",
                       trace="no",
                       backgro="no",
                       interac="no")
            print("Reference arc {} extracted".format(refarc))
            # identify the lines in it
            print("\nIdentify arc lines:")
            print("Enter the following in the splot window")
            print("\t:thres 500")
            print("\t:order 4, max = 5")
            print("\tfwidth 2")
            print("Select 3-5 arc lines from line atlas")
            print("Press 'm' to mark, then enter wavelength")
            print("Then press 'l' to automatically ID the other lines")
            print("Press 'f' to fit the dispersion correction")
            print("Use 'd' to remove bad points, 'f' to refit")
            print("'q' from fit, then 'q' from identify to continue\n")
            iraf.identify(images=refarc_out, coordlist=lineList_location)
        # use the refarc to ID all the subsequent arcs
        for arc in arclist:
            print("\nReidentifying arclines from {}".format(arc))
            iraf.reidentify(reference=refarc_out, images=arc)
        # add the refspec keywords to the image header for dispcor
        # refspec_factor tells IRAF how to interpolate the arcs
        refspec_factor = round((1./len(arclist)), 1)
        for i in range(0, len(arclist)):
            refspec = "{} {}".format(arclist[i].split(".fits")[0], refspec_factor)
            print("REFSPEC{}: {}".format(i+1, refspec))
            iraf.hedit(images=spectrum_out,
                       fields="REFSPEC{}".format(i+1),
                       value=refspec,
                       add="yes",
                       verify="no",
                       show="yes")
        print("Headers updated!\n")
        # apply the dispersion correction
        print("Applying the dispersion correction")
        iraf.dispcor(input=spectrum_out,
                     output=spectrum_out,
                     lineari="yes",
                     databas="database",
                     table="")
        print("Correction applied!")
        # normalize the spectrum using continuum
        normspec_out = "{}n.ms.fits".format(spectrum_out.split('.ms')[0])
        iraf.continuum(input=spectrum_out,
                       output=normspec_out,
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
    if not os.path.exists('spectra_r'):
        os.mkdir('spectra_r')
    for i in g.glob('*_tn.ms.fits'):
        os.system('mv {} spectra_r/'.format(i))

def eraseIntermediateProducts():
    """
    Delete all the partial steps, I never look at this crap
    and always start from scratch if there is an issue.

    Waste of disc space
    """
    for i in g.glob('i*.fits'):
        os.system('rm {}'.format(i))
    for i in g.glob('a*.fits'):
        os.system('rm {}'.format(i))

def utmiddleToNight(utmid):
    """
    Function to calculate night using the utmiddle
    """
    utmid = datetime.strptime(utmid, '%Y-%m-%dT%H:%M:%S.%f')
    if utmid.hour < 12:
        td = 1
    else:
        td = 0
    night = (utmid-timedelta(days=td)).strftime('%Y%m%d')
    return night

def logSpectraToDb():
    """
    Merged from LogSpectraToDB.py
    """
    db = pymysql.connect(host='localhost', db='eblm', password='mysqlpassword')
    os.chdir('spectra_r/')
    # get the list of spectra
    t = g.glob('*_tn.ms.fits')
    # loop over spectra and log them
    for i in t:
        h = fits.open(i)
        hdr = h[0].header
        d = h[0].data
        image_id = i
        object_name = hdr['OBJECT']
        hjd_mid = hdr['HJD-MID']
        bjd_mid = hdr['BJD-MID']
        jd_mid = hdr['JD-MID']
        utmiddle = hdr['UT-MID']
        pa = hdr['ROTSKYPA']
        n_traces = d.shape[1]
        night = utmiddleToNight(utmiddle)
        qry = """REPLACE INTO eblm_ids_newest
                (image_id,
                object_name,
                bjd_mid,
                hjd_mid,
                jd_mid,
                utmiddle,
                n_traces,
                sky_pa,
                night,
                analyse)
                VALUES
                ('{}', '{}', {}, {}, {}, '{}', {}, {}, '{}', 1)
                """.format(image_id,
                           object_name,
                           bjd_mid,
                           hjd_mid,
                           jd_mid,
                           utmiddle.replace("T", " "),
                           n_traces,
                           pa,
                           night)
        print(qry)
        with db.cursor() as cur:
            cur.execute(qry)
            db.commit()

def checkForDs9():
    """
    Find SpectorDS9
    """
    processes = os.popen('ps aux | grep SpectorDS9').readlines()
    for process in processes:
        if "ds9 -title SpectorDS9" in process:
            print('Hurray!')
            return True
    else:
        print('Boo!')
        return False

if __name__ == '__main__':
    print('\n\n---------------------------------------------------------')
    print('------------- Spectroscopy by Spector.py ----------------')
    print('---------------------------------------------------------\n')
    # get command line args
    args = argParse()
    if args.log_only:
        logSpectraToDb()
        sys.exit()
    if not os.path.exists(args.refarc):
        print('{} not found, quitting'.format(args.refarc))
        sys.exit(1)
    if args.ready:
        if args.ds9:
            # pyds9 was being a shit in astroconda env
            # beat it into submission with direct access via xpa
            # check for window
            already_ds9 = checkForDs9()
            if not already_ds9:
                os.system('ds9 -title SpectorDS9 &')
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
        # remove all the intermediate data products
        eraseIntermediateProducts()
        # log the spectra to the database
        logSpectraToDb()
    else:
        print('The --ready flag is needed to show that all\njunk images have been removed')
