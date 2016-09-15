"""
Code to measure radial velocites using iSpec 
and atomic line lists. Corrects for barycentric
velocity and telluric features before doing the
cross correlation

TODO:
    USE THE RIGHT MASK FOR EACH OBJECT!!!!
    READ THEM ALL THEN USE THE CLOSEST SPEC TYPE!!!!
    FIX THIS ASAP!!!

    Save out plots of:
        1. CCF for each spectrum
        2. Spectrum (various steps of reduction)
        3. Radial velocity summary for the object
"""
import os
import sys
import glob as g
import argparse as ap
import numpy as np
import pymysql
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
iSpec_location = '/home/virtual/shared/iSpec'
sys.path.insert(0, os.path.abspath(iSpec_location))
import ispec

# set up the locations of the various lines and mask files
strongLines = iSpec_location + '/input/regions/strong_lines/absorption_lines.txt'
telluricLines = iSpec_location + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
atomicMaskLines = iSpec_location + "/input/linelists/CCF/Narval.Sun.370_1048nm/mask.lst"
# "/input/linelists/CCF/Atlas.Arcturus.372_926nm/mask.lst""
# "/input/linelists/CCF/Atlas.Sun.372_926nm/mask.lst"
# "/input/linelists/CCF/HARPS_SOPHIE.A0.350_1095nm/mask.lst"
# "/input/linelists/CCF/HARPS_SOPHIE.F0.360_698nm/mask.lst"
# "/input/linelists/CCF/HARPS_SOPHIE.G2.375_679nm/mask.lst"
# "/input/linelists/CCF/HARPS_SOPHIE.K0.378_679nm/mask.lst"
# "/input/linelists/CCF/HARPS_SOPHIE.K5.378_680nm/mask.lst"
# "/input/linelists/CCF/HARPS_SOPHIE.M5.400_687nm/mask.lst"
# "/input/linelists/CCF/Synthetic.Sun.350_1100nm/mask.lst"
# "/input/linelists/CCF/VALD.Sun.300_1100nm/mask.lst"

# set up observatory
olat = 28.+(45./60.)-(37./3600.)
olon = -17.-(52./60.)-(46./3600.)
elev = 2332.
OBSERVATORY = EarthLocation(lat=olat*u.deg,lon=olon*u.deg,height=elev*u.m)

# have a kwarg dictionary for each instrument
INSTRUMENT = {'FIES': {'RA': 'OBJRA', \
                       'DEC': 'OBJDEC', \
                       'EXPTIME': 'EXPTIME', \
                       'WV_LLIM': 425.0, \
                       'WV_ULIM': 700.0, \
                       'RESOLUTION': 67000},
              'IDS': {'RA': 'CAT-RA', \
                      'DEC': 'CAT-DEC', \
                      'EXPTIME': 'EXPTIME', \
                      'WV_LLIM': 615.0, \
                      'WV_ULIM': 675.0, \
                      'RESOLUTION': 9500} # guestimate
              }

# acceptable comments for the CCF
valid_comments = ['sp', 'dp', 'np', 'mp', 'bp']

def argParse():
    parser = ap.ArgumentParser()
    parser.add_argument('instrument', help='instrument used for RVs', choices=['FIES', 'IDS'])
    parser.add_argument('--norm', help='normalise the continuum?', action='store_true')
    parser.add_argument('--plot', help='plot the different steps of analysis?', action='store_true')
    parser.add_argument('--wave', help='comma separated wavelength range to analyse - lower,upper', default='425.0,700.0')
    return parser.parse_args()

def getHostIP():
    """
    Get the VB host IP for database connection
    """
    return os.popen('netstat -rn').readlines()[2].split()[1]

def getSwaspIds():
    """
    Get a list of targets that require RV analysis
    These can be found by looking for NULL RV params in the database
    """
    swasp_ids = []
    qry = """
        SELECT distinct(swasp_id)
        FROM eblm_ids
        WHERE barycentric_velocity_iSpec IS NULL
        AND swasp_id IS NOT NULL
        AND n_traces = 1
        """
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            swasp_ids.append(row[0])
    return swasp_ids

def getDateObs(spec):
    """
    Return a tuple of utmid in the format:

    (year, month, day, hour, minute, second)

    for ispec.calculate_barycentric_velocity_correction
    """
    qry = """
        SELECT utmiddle, n_traces
        FROM eblm_ids
        WHERE image_id='{}'
        """.format(spec)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            utmid = Time(row[0], \
                      format='datetime', \
                      scale='utc', \
                      location=OBSERVATORY)
            n_traces = row[1]
    utmid_d = utmid.datetime
    year = utmid_d.year
    month = utmid_d.month
    day = utmid_d.day
    hour = utmid_d.hour
    minute = utmid_d.minute
    second = utmid_d.second
    return utmid, (year, month, day, hour, minute, second), n_traces

def getRaDec(spec):
    """
    Return a tuple of RA and Dec in the format:

    (hours, minutes, seconds, degrees, arcmin, arcsec)

    for ispec.calculate_barycentric_velocity_correction
    """
    with fits.open(spec) as fitsfile:
        ra = fitsfile[0].header[INSTRUMENT[args.instrument]['RA']]
        dec = fitsfile[0].header[INSTRUMENT[args.instrument]['DEC']]
    c = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
    ra_h = c.ra.hms.h
    ra_m = c.ra.hms.m
    ra_s = round(c.ra.hms.s, 2)
    dec_d = c.dec.dms.d
    dec_m = c.dec.dms.m
    dec_s = round(c.dec.dms.s, 2)
    return c, (ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)

def normaliseContinuum(spec):
    """
    Based on example.py
    normalize_whole_spectrum_strategy1_ignoring_prefixed_strong_lines function
    """
    model = 'Splines'
    degree = 2
    nknots = None
    from_resolution = INSTRUMENT[args.instrument]['RESOLUTION']

    # continuum fit
    order = 'median+max'
    median_wave_range = 0.01
    max_wave_range = 1.0

    strong_lines = ispec.read_line_regions(strongLines)
    continuum_model = ispec.fit_continuum(spec, \
                                          from_resolution=from_resolution, \
                                          nknots=nknots, \
                                          degree=degree, \
                                          median_wave_range=median_wave_range, \
                                          max_wave_range=max_wave_range, \
                                          model=model, \
                                          order=order, \
                                          automatic_strong_line_detection=True, \
                                          strong_line_probability=0.5, \
                                          use_errors_for_fitting=True)
    # continuum normalisation
    spec_norm = ispec.normalize_spectrum(spec, \
                                         continuum_model, \
                                         consider_continuum_errors=False)
    return spec_norm

def getBarycentricVelocity(dateobs, coords):
    """
    Get the barycentric velocity

    Based on example.py
    """
    return ispec.calculate_barycentric_velocity_correction(dateobs, coords)

def correctBarycentricVelocity(spectrum, barycentric_velocity):
    """
    Correct the spectrum to the barycentre
    """
    return ispec.correct_velocity(spectrum, barycentric_velocity)

def cleanTelluricRegions(spec):
    """
    Remove regions affected by tellurics

    Based on exmaple.py clean_telluric_regions() function
    """
    # determine telluric velocity shift
    linelist_telluric = ispec.read_telluric_linelist(telluricLines, minimum_depth=0.0)
    models, ccf = ispec.cross_correlate_with_mask(spec, \
                                                  linelist_telluric, \
                                                  lower_velocity_limit=-100, \
                                                  upper_velocity_limit=100, \
                                                  velocity_step=0.25, \
                                                  mask_size=2.0, \
                                                  mask_depth=0.01, \
                                                  fourier=False, \
                                                  only_one_peak=True)

    bv = np.round(models[0].mu(), 2) # km/s
    bv_err = np.round(models[0].emu(), 2) # km/s

    # clean regions that may be affected by tellurics
    min_vel = -30.0
    max_vel = +30.0
    # Only the 25% of the deepest ones:
    #dfilter = linelist_telluric['depth'] > np.percentile(linelist_telluric['depth'], 75)
    tfilter = ispec.create_filter_for_regions_affected_by_tellurics(spec['waveobs'], \
                                                                    linelist_telluric, \
                                                                    min_velocity=-bv+min_vel, \
                                                                    max_velocity=-bv+max_vel)
    clean_spec = spec[~tfilter]
    return bv, bv_err, clean_spec

def measureRadialVelocityWithMask(spec):
    """
    Radial velocity measurement using atomic line list

    Based on example.py determine_radial_velocity_with_mask() function
    """
    ccf_mask = ispec.read_linelist_mask(atomicMaskLines)
    models, ccf = ispec.cross_correlate_with_mask(spec, \
                                                  ccf_mask, \
                                                  lower_velocity_limit=-200, \
                                                  upper_velocity_limit=200, \
                                                  velocity_step=0.25, \
                                                  mask_depth=0.01, \
                                                  fourier=False)

    # Number of models represent the number of components
    components = len(models)
    # First component:
    rv = np.round(models[0].mu(), 2) # km/s
    rv_err = np.round(models[0].emu(), 2) # km/s
    return rv, rv_err, components, models, ccf

def getLightTravelTimes(target, time2corr):
    """
    Get the light travel times to the helio and
    barycentres

    Taken from NITES Process.py
    """
    ltt_bary = time2corr.light_travel_time(target)
    ltt_helio = time2corr.light_travel_time(target, 'heliocentric')
    return ltt_bary, ltt_helio

def getTargetParams(swasp_id):
    """
    Grab target parameters from the database
    """
    qry="""
        SELECT epoch, period
        FROM eblm_parameters
        WHERE swasp_id = '{}'
        LIMIT 1
        """.format(swasp_id)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            epoch = float(row[0]) + 2450000
            period = float(row[1])
    return epoch, period


def logRVsToDb(image_id,
           ccf_height,
           ccf_fwhm,
           atomic_rv,
           atomic_rv_err,
           v_tell,
           v_tell_err,
           barycentric_velocity,
           comment):
    """
    Insert the values from iSpec to the database
    """
    qry = """
        UPDATE eblm_ids SET
        ccf_height = {},
        ccf_fwhm = {},
        atomic_velocity = {},
        atomic_velocity_err = {},
        telluric_velocity = {},
        telluric_velocity_err = {},
        barycentric_velocity_iSpec = {},
        comment = '{}'
        WHERE
        image_id = '{}'
        """.format(ccf_height,
                   ccf_fwhm,
                   atomic_rv,
                   atomic_rv_err,
                   v_tell,
                   v_tell_err,
                   barycentric_velocity,
                   comment,
                   image_id)
    with db.cursor() as cur:
        cur.execute(qry)
        db.commit()

def getCurrentStatus(swasp_id):
    """
    Get current_status from parameters table
    """
    qry = """
        SELECT
        current_status
        FROM eblm_parameters
        WHERE swasp_id='{}'
        """.format(swasp_id)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            return row[0]

def updateTargetStatus(swasp_id, new_status):
    """
    Update a target's current status in parameters table
    """
    qry = """
        UPDATE eblm_parameters
        SET current_status='{}'
        WHERE swasp_id='{}'
        """.format(new_status, swasp_id)
    with db.cursor() as cur:
        cur.execute(qry)
        db.commit()

if __name__ == '__main__':
    # parse command line
    args = argParse()
    iSpec_dir = '/media/sf_all_spectra/'
    host_ip = getHostIP()
    db = pymysql.connect(host=host_ip,
                         db='eblm',
                         user='James',
                         password='mysqlpassword')
    # find the objects with NULL RV params in the database
    # loop over them and populate the databse. This method will allow for
    # adding more spectra later and then reducing only the newest ones.
    # add an option for just picking a particular object from the command line
    swasp_ids = getSwaspIds()
    for swasp_id in swasp_ids:
        target_dir = '{}/{}'.format(iSpec_dir, swasp_id)
        os.chdir(target_dir)
        # get the epoch and period from db
        epoch, period = getTargetParams(swasp_id)
        # get list of spectra
        spectra = g.glob('*.fits')
        print('Found {0:d} spectra for {1:s}'.format(len(spectra), swasp_id))
        results={}
        bjds, phase, rvs, rvs_errs = [], [], [], []
        ccf_heights, ccf_fwhms, comments = [], [], []
        for spectrum in spectra:
            print('Grabbing times from database for {}'.format(spectrum))
            # get the mid time in iSpec format
            utmid, dateobs, n_traces = getDateObs(spectrum)
            if n_traces > 1:
                print('BLEND, SKIPPING...')
                continue
            # get the target coords in iSpec format
            c, coords = getRaDec(spectrum)
            # get the barycentric velocity
            barycentric_velocity = getBarycentricVelocity(dateobs, coords)
            # bjd and hjd corrections
            ltt_bary, ltt_helio = getLightTravelTimes(c, utmid)
            bjd = utmid.tdb + ltt_bary
            hjd = utmid.utc + ltt_helio
            # read in the spectrum
            print('Reading {0:s}'.format(spectrum))
            spec = ispec.read_spectrum(spectrum)
            # correct the spectrum to the barycentre
            spec = correctBarycentricVelocity(spec, barycentric_velocity)
            # set up the plots
            if args.plot:
                fig, ax = plt.subplots(3,1, figsize=(10,10))
                ax[0].set_title('Spectral Analysis - {}'.format(swasp_id))
            # chop out the wavelength range we want to work on
            # accesss the data using column names
            # e.g. s['waveobs' | 'flux' | 'err']
            wave_filter = ispec.create_wavelength_filter(spec, \
                                          wave_base=INSTRUMENT[args.instrument]['WV_LLIM'], \
                                          wave_top=INSTRUMENT[args.instrument]['WV_ULIM'])
            spec = spec[wave_filter]
            if args.plot:
                ax[0].plot(spec['waveobs'], spec['flux'], 'r-')
                ax[0].set_xlabel('Wavelength (nm)')
            # normalise the continuum
            if args.norm:
                spec = normaliseContinuum(spec)
                if args.plot:
                    ax[0].plot(spec['waveobs'], spec['flux'], 'k-')
            # clean tellurics from the normalised spectrum
            v_tell, v_tell_err, spec = cleanTelluricRegions(spec)
            if args.plot:
                ax[0].plot(spec['waveobs'], spec['flux'], 'b-')
            # measure the radial velocity using an atomic line list
            atomic_rv, atomic_rv_err, components, models, ccf = measureRadialVelocityWithMask(spec)
            ccf_height = min(ccf['y'])
            ccf_heights.append(ccf_height)
            ccf_fwhm = models[0].sig()
            ccf_fwhms.append(ccf_fwhm)
            if args.plot:
                ax[2].set_title('CCF')
                ax[2].plot(ccf['x'],ccf['y'],'r-')
                ax[2].set_xlabel('Velocity (km/s)')
            results[bjd.jd] = (dateobs, \
                               barycentric_velocity, \
                               v_tell, \
                               v_tell_err, \
                               atomic_rv, \
                               atomic_rv_err)
            rvs.append(atomic_rv)
            rvs_errs.append(atomic_rv_err)
            bjds.append(bjd.jd)
            # python does the -ve wrap around check so no need to
            # do it manually like in javascript
            phase.append(((bjd.jd-epoch)/period)%1)
            if args.plot:
                # make some plots of the steps
                ax[1].set_title('Radial Velocity')
                ax[1].set_xlim(0, 1)
                ax[1].errorbar(phase, rvs, yerr=rvs_errs, fmt='r.')
                ax[1].set_xlabel('Orbital Phase')
                fig.subplots_adjust(hspace=0.4)
                plt.show()
            # print per-spectrum summary
            print(swasp_id,
                  spectrum,
                  round(atomic_rv, 4),
                  round(atomic_rv_err, 4),
                  round(ccf_height, 4),
                  round(ccf_fwhm, 4))
            # get a comment about the spectrum
            comment_valid = False
            while not comment_valid:
                comment = raw_input('Enter comment on spectrum: ')
                if comment.lower() in valid_comments:
                    comment_valid = True
                    comments.append(comment.lower())
            # log all the info plus comment to the database
            logRVsToDb(spectrum,
                       ccf_height,
                       ccf_fwhm,
                       atomic_rv,
                       atomic_rv_err,
                       v_tell,
                       v_tell_err,
                       barycentric_velocity,
                       comment)
        # print results and update the target in the database
        current_status = getCurrentStatus(swasp_id)
        print("\nSWASP_ID: {}".format(swasp_id))
        print("Epoch: {} Period: {}".format(round(epoch, 6), round(period, 6)))
        print("Current Status: {}".format(current_status))
        print("BJD  PHASE  RV  CCF_HEIGHT  CCF_FWHM  COMMENTS")
        for i, j, k, l, m, n in zip(bjds, phase, rvs, ccf_heights, ccf_fwhms, comments):
            print(i, round(j, 4), round(k, 4), round(l, 4), round(m, 4), n)
        update_status = raw_input('Update Status? (y/n): ')
        if update_status == 'y':
            new_status = raw_input('New Status: ')
            updateTargetStatus(swasp_id, new_status)
