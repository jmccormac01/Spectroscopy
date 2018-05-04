"""
Script to pull out HARPS CCFs

This was written to extract the CCFs for
NGTS-4b so I can fit both peaks in the CCF
"""

# TODO: Add normalisation of the CCF basline

import os
import sys
import glob as g
import argparse as ap
from datetime import (
    datetime,
    timedelta
    )
import numpy as np
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
olat = -29.-(15./60.)-(39.5/3600.)
olon = -70.-(43./60.)-(54.1/3600.)
elev = 2400.
OBSERVATORY = EarthLocation(lat=olat*u.deg, lon=olon*u.deg, height=elev*u.m)

# have a kwarg dictionary for each instrument
INSTRUMENT = {'HARPS': {'RA': 'RA', \
                        'DEC': 'DEC', \
                        'DATEOBS': 'DATE-OBS', \
                        'EXPTIME': 'EXPTIME', \
                        'RESOLUTION': 115000}}

def argParse():
    """
    Parse the command line arguments
    """
    parser = ap.ArgumentParser()
    parser.add_argument('instrument',
                        help='instrument used for RVs',
                        choices=['HARPS'])
    parser.add_argument('mask_type',
                        help='mask to use for extracting ccfs',
                        choices=atomicMaskLines.keys())
    parser.add_argument('data_dir',
                        help='location of the fits files')
    parser.add_argument('--norm',
                        help='normalise the continuum?',
                        action='store_true')
    parser.add_argument('--plot',
                        help='plot the different steps of analysis?',
                        action='store_true')
    return parser.parse_args()

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

def getDateObs(spec):
    """
    Return a tuple of utmid in the format:

    (year, month, day, hour, minute, second)

    for ispec.calculate_barycentric_velocity_correction
    """
    with fits.open(spec) as fitsfile:
        dateobs = fitsfile[0].header[INSTRUMENT[args.instrument]['DATEOBS']]
        dateobs = datetime.strptime(dateobs, "%Y-%m-%dT%H:%M:%S.%f")
        # correct the utstart to mid
        exptime = float(fitsfile[0].header[INSTRUMENT[args.instrument]['EXPTIME']])
        dateobs = dateobs + timedelta(seconds=exptime/2.)
        utmid = Time(dateobs, \
                  format='datetime', \
                  scale='utc', \
                  location=OBSERVATORY)
    year = dateobs.year
    month = dateobs.month
    day = dateobs.day
    hour = dateobs.hour
    minute = dateobs.minute
    second = dateobs.second
    return utmid, (year, month, day, hour, minute, second)

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

def cleanTelluricRegions(spec):
    """
    Remove regions affected by tellurics

    Based on exmaple.py clean_telluric_regions() function
    """
    # determine telluric velocity shift
    linelist_telluric = ispec.read_telluric_linelist(telluricLines, minimum_depth=0.0)
    models, _ = ispec.cross_correlate_with_mask(spec, \
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

if __name__ == '__main__':
    # parse command line
    args = argParse()
    os.chdir(args.data_dir)
    # get list of spectra
    spectra = sorted(g.glob('*.fits'))
    print('Found {} spectra in {}'.format(len(spectra), args.data_dir))

    # arrays/dicts to hold the results
    results = {}
    bjds, phase, comments = [], [], []
    mask_rvs, mask_rvs_errs = [], []
    mask_ccf_heights, mask_ccf_fwhms = [], []
    v_tells, v_tell_errs = [], []
    barycentric_velocities = []
    abs_rvs = []

    for spec_id, spectrum in enumerate(spectra):
        ccf_filename = '{}.ccf'.format(spectrum.split('.f')[0])
        if os.path.exists(ccf_filename):
            print('CCF file {} found, skipping'.format(ccf_filename))
            continue
        # get the mid time in iSpec format
        utmid, dateobs = getDateObs(spectrum)
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

        # read in the spectrum
        spec = ispec.read_spectrum(spectrum)
        print('{} Loaded...'.format(spectrum))

        # normalise the continuum
        if args.norm:
            spec = normaliseContinuum(spec)
            print('{} Continuum normalised...'.format(spectrum))

        # clean tellurics from the normalised spectrum
        v_tell, v_tell_err, spec = cleanTelluricRegions(spec)
        v_tells.append(v_tell)
        v_tell_errs.append(v_tell_err)
        print('{} Tellurics regions cleaned...'.format(spectrum))
        print('v_Tell = {} v_Tell_err = {}'.format(v_tell, v_tell_err))

        # measure the radial velocity using a atomic mask line list
        mask_rv, mask_rv_err, mask_components, mask_models, \
        mask_ccf = measureRadialVelocityWithMask(spec, args.mask_type)
        print('{} Cross correlated with {} mask...'.format(spectrum, args.mask_type))
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
        np.savetxt(ccf_filename,
                   np.c_[mask_ccf['x'], mask_ccf['y']],
                   fmt='%.3f  %.3f',
                   header='Velocity(km/s)  Contrast')
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
