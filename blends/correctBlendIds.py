"""
New script to disentangle the blends

To do:
    Add check for done objects
"""
import sys
import os
import urllib
from collections import defaultdict
from contextlib import contextmanager
import argparse as ap
import numpy as np
import pymysql
from astropy.io import fits
from astroquery.skyview import SkyView
import astropy.units as u
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.ndimage import rotate
from ds9 import (
    setupDs9,
    ds9Display
    )

DS9_NAME = 'BLENDS_DS9'

def argParse():
    """
    Parse command line arguments
    """
    p = ap.ArgumentParser()
    p.add_argument('--swasp_id',
                   help='swasp_id of single object to check')
    p.add_argument('--ds9',
                   help='display images in ds9?',
                   action='store_true')
    return p.parse_args()

@contextmanager
def openDb():
    with pymysql.connect(host='localhost', db='eblm',
                         password='mysqlpassword') as cur:
        yield cur

@contextmanager
def cd(path):
    """
    Improved os.chdir()
    """
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)

def getPartialBlendMap():
    """
    Read in a list of partially mapped blends to skip over
    """
    partial_map = []
    f = open('blend_map_ids_2017back.txt').readlines()
    for line in f:
        if not line.startswith('#'):
            partial_map.append(line.split()[0])
    return set(partial_map)

def getBlendedSpectra(swasp_id=None):
    blends = defaultdict(list)
    qry = """
        SELECT
        image_id, swasp_id, n_traces, sky_pa
        FROM eblm_ids_newest
        WHERE n_traces > 1
        """
    if swasp_id is not None:
        qry = qry + " AND swasp_id='{}'".format(swasp_id)
    with openDb() as cur:
        cur.execute(qry)
        results = cur.fetchall()
    for row in results:
        blends[row[1]].append({'image_id': row[0],
                               'n_traces': int(row[2]),
                               'sky_pa': int(row[3])})
    return blends

def plotUnnormalisedSpectra(ax, image_id, n_traces):
    """
    Plot the spectra in their order in the file
    so we can disentangle their swasp_ids
    """
    colours = ['k-', 'r-', 'b-']
    ax.set_ylabel('Flux')
    ax.set_xlabel('Wavelength')
    ax.set_title(image_id)
    with fits.open(image_id) as hdu:
        for i in range(n_traces):
            unnorm_spec = hdu[0].data[1][i]
            ax.plot(unnorm_spec, colours[i])
    ax.legend(tuple(range(n_traces)))
    return ax

def checkCache(cache_id):
    """
    Look for a cached image of this object
    """
    try:
        data = fits.open(cache_id)[0].data
        print('Found cached image {}'.format(cache_id))
    except FileNotFoundError:
        data = None
        print('No cached image found for {}'.format(cache_id))
    return data

def correctPa(angle):
    """
    Correct the PA for instrument config
    """
    return angle - 90

def drawSlit(ax, data, slit_width=20):
    """
    Work out the slit position from the image data
    """
    y, x = data.shape
    slit_upper_left = (0, int(y/2)-slit_width/2)
    rect = patches.Rectangle(slit_upper_left,
                             x-1, slit_width,
                             linewidth=1,
                             edgecolor='r',
                             facecolor='none')
    ax.add_patch(rect)
    return ax

def getSwaspIdRaDec(object_id):
    """
    Grab coordinates for swasp object
    """
    if object_id.startswith('1SWASP'):
        qry = """
            SELECT ra_hms, dec_dms
            FROM eblm_parameters
            WHERE swasp_id=%s
            """
    elif object_id.startswith('HD'):
        qry = """
            SELECT ra_hms, dec_dms
            FROM rv_standards
            WHERE object_id=%s
            """
    else:
        print('Object ID not understood, skipping...')
        return None, None
    with openDb() as cur:
        cur.execute(qry, (object_id,))
        result = cur.fetchone()
    ra = result[0]
    dec = result[1]
    return ra, dec

def plotFindingChart(ax, object_id, angle):
    """
    Make a plot of the field with the slit overlaid
    """
    ra, dec = getSwaspIdRaDec(object_id)
    if ra is None or dec is None:
        sys.exit(1)
    position = "{} {}".format(ra, dec)
    coordinates = "J2000"
    surveys = ["SDSSr", "SDSSi", "DSS2 Red", "DSS2 Blue"]
    # dimensions in arcmin
    width = 3
    height = 3
    # check for a cached image
    cache_location = "/Users/jmcc/Dropbox/PythonScripts/Finders/cache/"
    cache_id = "{}/{}_{}_{}.fits".format(cache_location,
                                         object_id,
                                         height,
                                         width)
    data = checkCache(cache_id)
    if data is None:
        sv = SkyView()
        for survey in surveys:
            print(survey)
            try:
                path = sv.get_images(position=position,
                                     coordinates=coordinates,
                                     survey=survey,
                                     width=width*u.arcmin,
                                     height=height*u.arcmin)
                data = path[0][0].data
                hdr = path[0][0].header
                # write out the fits file into the cache
                fits.writeto(cache_id, data, header=hdr)
                break
            except urllib.error.HTTPError:
                print('No image in {} for {}'.format(survey, position))
                continue
        # if we get here there are no images!
        else:
            print('No images found for {}, exiting...'.format(object_id))
            sys.exit(1)
    # if we make it here plot the finder
    ndata = rotate(data, correctPa(angle))
    # cameras need to have the images fliped about the Y axis
    # to match up with DS9 x and y axes
    ndata =  np.flip(ndata, 1)
    ax.imshow(ndata, origin='lower')
    ax = drawSlit(ax, ndata)
    ax.set_title('SKY_PA = {}'.format(angle))
    return ax

if __name__ == "__main__":
    args = argParse()
    top_dir = "/Users/jmcc/Dropbox/data/int/ids/eblm"
    all_spec_dir = "{}/all_spectra".format(top_dir)
    if args.ds9:
        setupDs9(DS9_NAME)
    # get blends for specific object or all of the blends
    if args.swasp_id:
        blends = getBlendedSpectra(args.swasp_id)
    else:
        blends = getBlendedSpectra()
    # get a list of partially mapped blends, if any
    partial_blend_map = getPartialBlendMap()
    if len(blends) > 0:
        # now go into each swasp directory and pull up the
        # finding chart for that object and cycle through the blends spectra
        os.chdir(all_spec_dir)
        for i, blend in enumerate(blends):
            if blend in partial_blend_map:
                print('{} blends mapped already, skipping...'.format(blend))
                continue
            with cd(blend):
                n_blends = len(blends[blend])
                print('[{}:{}] There are {} blended spectra for {}'.format(i+1,
                                                                           len(blends),
                                                                           n_blends,
                                                                           blend))
                # loop over each blended spectrum
                for j, row in enumerate(blends[blend]):
                    image_id = row['image_id']
                    sky_pa = row['sky_pa']
                    n_traces = row['n_traces']
                    if args.ds9:
                        # display the DS9 frame for this spectrum
                        ds9Set(DS9_NAME, 'frame show {}'.format(j+1))
                        orig_spec = image_id.split('_')[2]
                        original_spectrum = "{}/{}/original/{}.fits*".format(top_dir,
                                                                             night_dir,
                                                                             orig_spec)
                        ds9Display(DS9_NAME, original_spectrum)
                    fig, ax = plt.subplots(1, 2, figsize=(15, 5))
                    fig.suptitle(blend)
                    ax[0] = plotUnnormalisedSpectra(ax[0], image_id, n_traces)
                    ax[1] = plotFindingChart(ax[1], blend, sky_pa)
                    png_file = "{}_{}_{}.png".format(image_id, sky_pa, n_traces)
                    plt.subplots_adjust(top=0.80)
                    plt.savefig(png_file, dpi=300)
                plt.show()
    else:
        print('There are no blends to sort out, exiting.')
        sys.exit(0)
