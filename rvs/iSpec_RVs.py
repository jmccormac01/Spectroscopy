"""
Code to measure radial velocites using iSpec
and atomic line lists. Corrects for barycentric
velocity and telluric features before doing the
cross correlation

Steps:
    1. Measure barycentric velocity
    1. Measure velocity wrt Tellurics
    1. Measure RVs with atomic line mask
    1. Sum these 4 velocities to get the final
       radial velocity
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
iSpec_location = '/home/virtual/iSpec'
sys.path.insert(0, os.path.abspath(iSpec_location))
import ispec

# pylint: disable = superfluous-parens
# pylint: disable = invalid-name
# pylint: disable = no-member
# pylint: disable = redefined-outer-name

# set up the locations of the various lines and mask files
strongLines = iSpec_location + '/input/regions/strong_lines/absorption_lines.txt'
line_lists_parent = iSpec_location + '/input/linelists/CCF/'
telluricLines = line_lists_parent + "Synth.Tellurics.500_1100nm/mask.lst"
atomicMaskLines = {'A0': line_lists_parent + 'HARPS_SOPHIE.A0.350_1095nm/mask.lst',
                   'F0': line_lists_parent + 'HARPS_SOPHIE.F0.360_698nm/mask.lst',
                   'G2': line_lists_parent + 'HARPS_SOPHIE.G2.375_679nm/mask.lst',
                   'K0': line_lists_parent + 'HARPS_SOPHIE.K0.378_679nm/mask.lst',
                   'K5': line_lists_parent + 'HARPS_SOPHIE.K5.378_680nm/mask.lst',
                   'M5': line_lists_parent + 'HARPS_SOPHIE.M5.400_687nm/mask.lst'}

# set up observatory
olat = 28.+(45./60.)-(37./3600.)
olon = -17.-(52./60.)-(46./3600.)
elev = 2332.
OBSERVATORY = EarthLocation(lat=olat*u.deg, lon=olon*u.deg, height=elev*u.m)

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
# sp = single peaked
# dp = double peaked
# np = no peak
# mp = multiple peaked
# bp = broad peaked
# ap = asymetrical peaked
valid_comments = ['sp', 'dp', 'np', 'mp', 'bp', 'ap']

def argParse():
    """
    Parse the command line arguments
    """
    parser = ap.ArgumentParser()
    parser.add_argument('instrument',
                        help='instrument used for RVs',
                        choices=['FIES', 'IDS'])
    parser.add_argument('--norm',
                        help='normalise the continuum?',
                        action='store_true')
    parser.add_argument('--plot',
                        help='plot the different steps of analysis?',
                        action='store_true')
    parser.add_argument('--wave',
                        help='comma separated wavelength range to analyse - lower,upper',
                        default='425.0,700.0')
    parser.add_argument('--swasp_id',
                        help='swasp_id of single object to measure')
    parser.add_argument('--no_comment',
                        help=('skip adding comments to individual spectra, '
                              'only to be used for debugging'),
                        action='store_true')
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
        FROM eblm_ids_newest
        WHERE barycentric_velocity_iSpec IS NULL
        AND swasp_id IS NOT NULL
        AND analyse = 1
        ORDER BY swasp_id
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
        FROM eblm_ids_newest
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
    c = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.degree), frame='icrs')
    ra_h = c.ra.hms.h
    ra_m = c.ra.hms.m
    ra_s = round(c.ra.hms.s, 2)
    dec_d = c.dec.dms.d
    dec_m = c.dec.dms.m
    dec_s = round(c.dec.dms.s, 2)
    return c, (ra_h, ra_m, ra_s, dec_d, dec_m, dec_s)

def getListOfSpectralTypes():
    """
    Return a list of all spectral types
    """
    # list of all spectral types
    SPEC_TYPES = ['A', 'F', 'G', 'K', 'M']
    spec_types_list = []
    for typ in SPEC_TYPES:
        for i in range(0, 10):
            spec_types_list.append('{}{}'.format(typ, i))
    return spec_types_list

def getAtomicLineMaskIndexes(atomicMaskLines):
    """
    Get the indexes of the available atomic line
    masks from the full list of spectral types
    """
    atomicMaskLinesIndexes = {}
    spec_types = getListOfSpectralTypes()
    for key in atomicMaskLines:
        index = spec_types.index(key)
        atomicMaskLinesIndexes[key] = index
    return spec_types, atomicMaskLinesIndexes

def getCcfMaskType(swasp_id, spec_types, maskIndexes):
    """
    Find the closest spectral type to use for the RVs
    """
    if swasp_id.startswith('HD'):
        qry = """
            SELECT spec_type
            FROM rv_standards
            WHERE object_id = '{}'
            """.format(swasp_id)
    else:
        qry = """
            SELECT paramfit_spec_type
            FROM eblm_parameters
            WHERE swasp_id = '{}'
            """.format(swasp_id)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            spec_type = row[0]
    index = spec_types.index(spec_type)
    diff = 1E6
    for mask in maskIndexes:
        if abs(maskIndexes[mask] - index) < diff:
            diff = abs(maskIndexes[mask]-index)
            match = mask
    print ('Spectral Type: {} - Mask Match: {}'.format(spec_type, match))
    return spec_type, match

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

#def correctBarycentricVelocity(spectrum, barycentric_velocity):
#    """
#    Correct the spectrum to the barycentre
#    """
#    return ispec.correct_velocity(spectrum, barycentric_velocity)

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
                                                  velocity_step=0.10, \
                                                  mask_size=2.0, \
                                                  mask_depth=0.01, \
                                                  fourier=False, \
                                                  only_one_peak=True)
    try:
        bv = np.round(models[0].mu(), 2) # km/s
        bv_err = np.round(models[0].emu(), 2) # km/s
    except IndexError:
        print '\n\n\nPROBLEM MEASURING TELLURICS, SKIPPING...\n\n\n'
        return 0.0, 0.0, spec
    # not sure why we load this again, but iSpec example does, so...
    linelist_telluric = ispec.read_telluric_linelist(telluricLines, minimum_depth=0.0)
    # clean regions that may be affected by tellurics
    min_vel = -30.0
    max_vel = +30.0
    # Only the 25% of the deepest ones:
    dfilter = linelist_telluric['depth'] > np.percentile(linelist_telluric['depth'], 75)
    tfilter = ispec.create_filter_for_regions_affected_by_tellurics(spec['waveobs'], \
                                                                    linelist_telluric[dfilter],\
                                                                    min_velocity=-bv+min_vel, \
                                                                    max_velocity=-bv+max_vel)
    clean_spec = spec[~tfilter]
    return bv, bv_err, clean_spec

def measureRadialVelocityWithMask(spec, mask_type):
    """
    Radial velocity measurement using atomic line list

    Based on example.py determine_radial_velocity_with_mask() function
    """
    ccf_mask = ispec.read_cross_correlation_mask(atomicMaskLines[mask_type])
    models, ccf = ispec.cross_correlate_with_mask(spec, \
                                                  ccf_mask, \
                                                  lower_velocity_limit=-200, \
                                                  upper_velocity_limit=200, \
                                                  velocity_step=0.10, \
                                                  mask_depth=0.01, \
                                                  fourier=False)

    # Number of models represent the number of components
    components = len(models)
    # First component:
    try:
        rv = np.round(models[0].mu(), 2) # km/s
        rv_err = np.round(models[0].emu(), 2) # km/s
    except IndexError:
        print '\n\n\nPROBLEM RV WITH MASK, SKIPPING...\n\n\n'
        return 0.0, 0.0, 0, None, None
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
    qry = """
        SELECT epoch, period, paramfit_spec_type
        FROM eblm_parameters
        WHERE swasp_id = '{}'
        LIMIT 1
        """.format(swasp_id)
    with db.cursor() as cur:
        cur.execute(qry)
        result = cur.fetchone()
        if result is not None:
            epoch = float(result[0]) + 2450000
            period = float(result[1])
            spectral_type = result[2]
            return epoch, period, spectral_type
        else:
            return None, None, None

def logRVsToDb(spectrum,
               mask,
               mask_ccf_height,
               mask_ccf_fwhm,
               mask_rv,
               mask_rv_err,
               v_tell,
               v_tell_err,
               barycentric_velocity,
               comment):
    """
    Insert the values from iSpec to the database
    """
    qry = """
        UPDATE eblm_ids_newest SET
        mask = '{}',
        mask_ccf_height = {},
        mask_ccf_fwhm = {},
        mask_velocity = {},
        mask_velocity_err = {},
        telluric_velocity = {},
        telluric_velocity_err = {},
        barycentric_velocity_iSpec = {},
        comment = '{}'
        WHERE
        image_id = '{}'
        """.format(mask,
                   mask_ccf_height,
                   mask_ccf_fwhm,
                   mask_rv,
                   mask_rv_err,
                   v_tell,
                   v_tell_err,
                   barycentric_velocity,
                   comment,
                   spectrum)
    print(qry)
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
                         user='jmcc',
                         password='mysqlpassword')
    # set up wavelength range for this instrument
    wave_base = INSTRUMENT[args.instrument]['WV_LLIM']
    wave_top = INSTRUMENT[args.instrument]['WV_ULIM']
    # get a list of spectral types and the atomic mask
    # line list indexes. this is used to ID the best mask
    spec_types_list, atomicMaskLinesIndexes = getAtomicLineMaskIndexes(atomicMaskLines)
    # find the objects with NULL RV params in the database
    # loop over them and populate the databse. This method will allow for
    # adding more spectra later and then reducing only the newest ones.
    # here we check for whether a particular object has been selected or not
    if not args.swasp_id:
        swasp_ids = getSwaspIds()
    else:
        swasp_ids = [args.swasp_id]

    for swasp_id in swasp_ids:
        # get the epoch and period from db
        if not swasp_id.startswith('HD'):
            epoch, period, spectral_type = getTargetParams(swasp_id)
            print('{} E={} P={} ST={}'.format(swasp_id, epoch, period, spectral_type))
            if epoch is None or period is None or spectral_type is None:
                print('Skipping {}...'.format(swasp_id))
                continue
        target_dir = '{}/{}'.format(iSpec_dir, swasp_id)
        os.chdir(target_dir)
        # get list of spectra
        spectra = g.glob('*.fits')
        print('Found {} spectra for {}'.format(len(spectra), swasp_id))
        # get the mask type for this object
        spec_type, mask_type = getCcfMaskType(swasp_id, spec_types_list, atomicMaskLinesIndexes)

        # arrays/dicts to hold the results
        results = {}
        bjds, phase, comments = [], [], []
        mask_rvs, mask_rvs_errs = [], []
        mask_ccf_heights, mask_ccf_fwhms = [], []
        v_tells, v_tell_errs = [], []
        barycentric_velocities = []
        abs_rvs = []

        # loop over the spectra
        for spec_id, spectrum in enumerate(spectra):
            # get the mid time in iSpec format
            utmid, dateobs, n_traces = getDateObs(spectrum)
            print('{} Queried times from database...'.format(spectrum))
            # get the target coords in iSpec format
            c, coords = getRaDec(spectrum)
            print('{} Gathered target coordinates...'.format(spectrum))
            # get the barycentric velocity
            barycentric_velocity = getBarycentricVelocity(dateobs, coords)
            barycentric_velocities.append(barycentric_velocity)
            print('{} Barycentric velocity {} km/s...'.format(spectrum, barycentric_velocity))
            # bjd and hjd corrections
            ltt_bary, ltt_helio = getLightTravelTimes(c, utmid)
            bjd = utmid.tdb + ltt_bary
            hjd = utmid.utc + ltt_helio
            print('{} Light travel times calculated...'.format(spectrum))

            # set up the plots
            if args.plot:
                plt.ion()
                fig = plt.figure(1, figsize=(10, 10))
                ax_spec = plt.subplot2grid((4, 4), (0, 0), colspan=4, rowspan=2)
                ax_spec.set_title('{} {}'.format(swasp_id, spectrum))
                ax_mask_rv = plt.subplot2grid((4, 4), (2, 0), colspan=2, rowspan=2)
                ax_mask_rv.set_title('RV wrt {} (w/ bary_corr)'.format(mask_type))
                ax_mask_ccf = plt.subplot2grid((4, 4), (2, 2), colspan=2, rowspan=2)
                ax_mask_ccf.set_title('CCF wrt {}'.format(mask_type))
                plt.subplots_adjust(hspace=1.0)
                # figure for each ccf
                fig_ccf, ax_ccf = plt.subplots(1, figsize=(10, 10))
                ax_ccf.set_title('CCF wrt {}'.format(mask_type))
                plt.show()

            # read in the spectrum
            spec = ispec.read_spectrum(spectrum)
            print('{} Loaded...'.format(spectrum))

            # chop out the wavelength range we want to work on
            # accesss the data using column names
            # e.g. s['waveobs' | 'flux' | 'err']
            wave_filter = ispec.create_wavelength_filter(spec, \
                                                         wave_base=wave_base, \
                                                         wave_top=wave_top)
            spec = spec[wave_filter]
            print('{} Wavelength range restricted to {}-{}'.format(spectrum, \
                                                                   wave_base, \
                                                                   wave_top))
            if args.plot:
                ax_spec.plot(spec['waveobs'], spec['flux'], 'r-')

            # normalise the continuum
            if args.norm:
                spec = normaliseContinuum(spec)
                print('{} Continuum normalised...'.format(spectrum))
                if args.plot:
                    ax_spec.plot(spec['waveobs'], spec['flux'], 'r-')

            # clean tellurics from the normalised spectrum
            v_tell, v_tell_err, spec = cleanTelluricRegions(spec)
            v_tells.append(v_tell)
            v_tell_errs.append(v_tell_err)
            print('{} Tellurics regions cleaned...'.format(spectrum))
            print('v_Tell = {} v_Tell_err = {}'.format(v_tell, v_tell_err))
            if args.plot:
                ax_spec.plot(spec['waveobs'], spec['flux'], 'b-')
                ax_spec.set_xlabel('Wavelength (nm)')

            # measure the radial velocity using a atomic mask line list
            mask_rv, mask_rv_err, mask_components, mask_models, \
            mask_ccf = measureRadialVelocityWithMask(spec, mask_type)
            print('{} Cross correlated with {} mask...'.format(spectrum, mask_type))
            if mask_components < 1:
                mask_ccf_height = 1.0
                mask_ccf_heights.append(mask_ccf_height)
                mask_ccf_fwhm = 50.0
                mask_ccf_fwhms.append(mask_ccf_fwhm)
            else:
                mask_ccf_height = min(mask_ccf['y'])
                mask_ccf_heights.append(mask_ccf_height)
                mask_ccf_fwhm = mask_models[0].sig()
                mask_ccf_fwhms.append(mask_ccf_fwhm)
                if args.plot:
                    ax_mask_ccf.plot(mask_ccf['x'], mask_ccf['y'], 'r-')
                    ax_mask_ccf.set_xlabel('Velocity (km/s)')
                    ax_ccf.plot(mask_ccf['x'], mask_ccf['y'], 'r-')
                    ax_ccf.set_xlabel('Velocity (km/s)')
            mask_rvs.append(mask_rv)
            mask_rvs_errs.append(mask_rv_err)
            abs_rvs.append(mask_rv + barycentric_velocity)
            # make a dictionary of the results
            results[bjd.jd] = (dateobs, \
                               barycentric_velocity, \
                               v_tell, \
                               v_tell_err, \
                               mask_rv, \
                               mask_rv_err)
            bjds.append(bjd.jd)

            if not swasp_id.startswith('HD'):
                # python does the -ve wrap around check so no need to
                # do it manually like in javascript
                phase.append(((bjd.jd-epoch)/period)%1)
                if args.plot:
                    # make some plots of the steps
                    ax_mask_rv.set_xlim(0, 1)
                    ax_mask_rv.errorbar(phase, abs_rvs, yerr=mask_rvs_errs, fmt='r.')
                    ax_mask_rv.set_xlabel('Orbital Phase')
                    plt.show()

            # print per-spectrum mask summary
            print('Mask: {}'.format(mask_type),
                  swasp_id,
                  spectrum,
                  round(mask_rv, 4),
                  round(mask_rv_err, 4),
                  round(mask_ccf_height, 4),
                  round(mask_ccf_fwhm, 4))

            if not args.no_comment:
                # get a comment about the spectrum
                comment_valid = False
                while not comment_valid:
                    comment = raw_input('Enter comment on spectrum: ')
                    if comment.lower() in valid_comments:
                        comment_valid = True
                        comments.append(comment.lower())
            else:
                comments.append('nc')

            # log all the info plus comment to the database
            logRVsToDb(spectrum,
                       mask_type,
                       mask_ccf_height,
                       mask_ccf_fwhm,
                       mask_rv,
                       mask_rv_err,
                       v_tell,
                       v_tell_err,
                       barycentric_velocity,
                       comment)
            # save the ccf figure
            fig_ccf.savefig('{}_ccf_{}.png'.format(spectrum.split('.fits')[0],
                                                   mask_type), dpi=400)
            plt.close(fig_ccf)
        # don't do this for standard stars
        if not swasp_id.startswith('HD'):
            # print results and update the target in the database
            current_status = getCurrentStatus(swasp_id)
            print("\nSWASP_ID: {} SpecType: {}  Mask: {}".format(swasp_id,
                                                                 spec_type,
                                                                 mask_type))
            print("Epoch: {} Period: {}".format(round(epoch, 6), round(period, 6)))
            print("Current Status: {}".format(current_status))
            print("BJD             PHASE   BCV    V_TELL   MASK  CCF  FWHM  RV   COMMENTS")
            for i, j, k, l, m, n, o, p in zip(bjds,
                                              phase,
                                              barycentric_velocities,
                                              v_tells,
                                              mask_rvs,
                                              mask_ccf_heights,
                                              mask_ccf_fwhms,
                                              comments):
                print(round(i, 5), round(j, 4), round(k, 4),
                      round(l, 4), round(m, 4), round(n, 4),
                      round(o, 4), round(m + k, 4), p)
            if not args.no_comment:
                update_status = raw_input('Update Status? (y/n): ')
                if update_status == 'y':
                    new_status = raw_input('New Status: ')
                    updateTargetStatus(swasp_id, new_status)
        # do this instead
        else:
            print('\nID: {}'.format(swasp_id))
            print("BJD  BCV  V_TELL  MASK  CCF  FWHM  RV   COMMENTS")
            for i, j, k, l, m, n, o in zip(bjds,
                                           barycentric_velocities,
                                           v_tells,
                                           mask_rvs,
                                           mask_ccf_heights,
                                           mask_ccf_fwhms,
                                           comments):
                print(round(i, 5), round(j, 4), round(k, 4),
                      round(l, 4), round(m, 4), round(n, 4),
                      round(j + l, 4), o)
        plt.cla()
