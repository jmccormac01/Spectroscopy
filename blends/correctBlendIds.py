"""
New script to disentangle the blends

To do:
    Add code to split the traces up and log to the db
"""
import sys
import os
import time
import urllib
import glob as g
from collections import defaultdict
from contextlib import contextmanager
import argparse as ap
import numpy as np
import pymysql
from astropy.io import fits
import astropy.units as u
from astroquery.skyview import SkyView
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.ndimage import rotate
from ds9 import (
    setupDs9,
    ds9Display,
    ds9Set
    )

# pylint: disable=invalid-name
# pylint: disable=superfluous-parens
# pylint: disable=redefined-outer-name

DS9_NAME = 'BLENDS_DS9'

# set up some exceptions to work cross Python2 and 3
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

def argParse():
    """
    Parse command line arguments
    """
    p = ap.ArgumentParser()
    p.add_argument('--filter_blends',
                   help='look through all the blends manually check them',
                   action='store_true')
    p.add_argument('--split_blends',
                   help='split the previously filtered blended multispecs',
                   action='store_true')
    p.add_argument('--split_apply',
                   help='apply the splits',
                   action='store_true')
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
        image_id, swasp_id, n_traces, sky_pa, night
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
                               'sky_pa': int(row[3]),
                               'night': row[4]})
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
    ndata = np.flip(ndata, 1)
    ax.imshow(ndata, origin='lower')
    ax = drawSlit(ax, ndata)
    ax.set_title('SKY_PA = {}'.format(angle))
    return ax

def parseBlendMapFile(blend_map):
    """
    This file maps the original swasp_id and image_id to the the new
    swasp_ids for each trace. This file is made manually in the filter_blends
    part of this script

    Note the original swasp_id can be wrong if the wrong object was observed
    After the splitting the swasp_id for each trace is the right ID.
    """
    image_list = defaultdict(list)
    trace_map = {}
    f = open(blend_map).readlines()
    for line in f:
        if not line.startswith('#'):
            swasp_id, image_id, t1, t2, t3 = line.split()
            if t1 == 'None':
                t1 = None
            if t2 == 'None':
                t2 = None
            if t3 == 'None':
                t3 = None
            image_list[swasp_id].append(image_id)
            trace_map[image_id] = {1: t1,
                                   2: t2,
                                   3: t3}
    return image_list, trace_map

def scopy(iraf, image_id, trace_id, new_image_prefix):
    """
    Use IRAF scopy to extract a single trace from the multispec file
    """
    iraf.scopy.setParam('w1', 'INDEF')
    iraf.scopy.setParam('w2', 'INDEF')
    iraf.scopy.setParam('apmodulus', '0')
    iraf.scopy.setParam('format', 'onedspec')
    iraf.scopy.setParam('renumber', 'no')
    iraf.scopy.setParam('offset', '0')
    iraf.scopy.setParam('clobber', 'no')
    iraf.scopy.setParam('merge', 'no')
    iraf.scopy.setParam('rebin', 'no')
    iraf.scopy(input=image_id, output=new_image_prefix, aperture=trace_id)

def getOldBlendTableRow(image_id):
    """
    Grab the information from the old blended row that can be reused in
    each of the new per-trace rows
    """
    qry = """
        SELECT
        bjd_mid, hjd_mid, jd_mid, utmiddle, sky_pa, night
        FROM eblm_ids_newest
        WHERE image_id = %s
        """
    with openDb() as cur:
        cur.execute(qry, (image_id,))
        results = cur.fetchone()
    old_row_info = {'bjd': results[0],
                    'hjd': results[1],
                    'jd': results[2],
                    'ut': results[3].strftime('%Y-%m-%d %H:%M:%S'),
                    'pa': results[4],
                    'night': results[5].strftime('%Y-%m-%d')}
    return old_row_info

def updateOldBlendTableRow(image_id):
    """
    Set the old row in the table analyse=-1 to show this was a blended spectrum
    and now we don't want to analyse this image. At a later step we will
    include a new row per trace
    """
    qry = """
        UPDATE eblm_ids_newest
        SET analyse=-1
        WHERE image_id=%s
        """
    with openDb() as cur:
        cur.execute(qry, (image_id,))

def addNewUnblendedTableRow(new_image_prefix, new_swasp_id, old_row_info):
    """
    Insert a new table row for the unblended spectra
    """
    qry = """
        INSERT INTO eblm_ids_newest
        (image_id, swasp_id, object_name, bjd_mid, hjd_mid, jd_mid,
        utmiddle, n_traces, sky_pa, analyse, night)
        VALUES
        (%s, %s, %s, %s, %s, %s, %s, 1, %s, 1, %s)
        """
    qry_args = (new_image_prefix,
                new_swasp_id,
                new_swasp_id,
                old_row_info['bjd'],
                old_row_info['hjd'],
                old_row_info['jd'],
                old_row_info['ut'],
                old_row_info['pa'],
                old_row_info['night'],)
    with openDb() as cur:
        cur.execute(qry, qry_args)

if __name__ == "__main__":
    args = argParse()
    loginCl_location = "/Users/jmcc/"
    top_dir = "/Users/jmcc/Dropbox/data/int/ids/eblm"
    code_dir = "/Users/jmcc/Dropbox/PythonScripts/Spectroscopy/blends"
    all_spec_dir = "{}/all_spectra".format(top_dir)
    all_blends_dir = "{}/all_blends".format(all_spec_dir)
    map_file = "{}/blend_map_ids_2017back.txt".format(code_dir)
    # manually go through all the blends and match swasp_ids to traces
    # when a map file is made, run the splitting section of this script
    # on that file to divide up the 1D spectra for measuring RVs
    if args.filter_blends:
        if args.ds9:
            setupDs9(DS9_NAME)
        # get blends for specific object or all of the blends
        if args.swasp_id:
            blends = getBlendedSpectra(args.swasp_id)
        else:
            blends = getBlendedSpectra()
        # get a list of partially mapped blends, if any
        partial_blend_map = getPartialBlendMap()
        # this section is for sorting the blends into their
        # correct swasp_ids. Make a list of matching traces to swasp_ids
        # see blend_map_ids_2017back for an example
        # once this file is complete for all blends, we split the multispec files
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
                        night = row['night'].strftime('%Y%m%d')
                        if args.ds9:
                            # display the DS9 frame for this spectrum
                            ds9Set(DS9_NAME, 'frame frameno {}'.format(j+1))
                            orig_spec = image_id.split('_')[2]
                            spec_dir = "{}/{}/original/".format(top_dir,
                                                                night)
                            spec_list = g.glob('{}/{}*'.format(spec_dir, orig_spec))
                            if len(spec_list) > 0:
                                ds9Display(DS9_NAME, spec_list[0])
                            else:
                                print('No spectrum found for {}/{}*'.format(spec_dir,
                                                                            orig_spec))
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

    # work on splitting the blends filtered in the block above
    if args.split_blends:
        # load IRAF
        here = os.getcwd()
        os.chdir(loginCl_location)
        from pyraf import iraf
        # load onedspec module
        iraf.onedspec(_doprint=0)
        os.chdir(here)
        time.sleep(2)
        # read in the map file
        image_list, trace_map = parseBlendMapFile(map_file)
        # loop over each object and each image
        os.chdir(all_spec_dir)
        for swasp_id in image_list:
            with cd(swasp_id):
                # loop over each blended spectrum
                for image_id in image_list[swasp_id]:
                    if not os.path.exists(image_id):
                        print('{} does not exist. Already processed?'.format(image_id))
                        continue
                    # get the old row info once per blended spectrum
                    old_row_info = getOldBlendTableRow(image_id)
                    print(old_row_info)
                    # new image prefix
                    new_image_prefix = "{}".format(image_id.split('.')[0])
                    # set a flag to only update the original table row once
                    original_row_updated = False
                    # split the spectrum into its traces
                    for trace_id in trace_map[image_id]:
                        new_swasp_id = trace_map[image_id][trace_id]
                        if new_swasp_id is not None:
                            new_target_dir = "{}/{}".format(all_spec_dir,
                                                            new_swasp_id)
                            print('Splitting {} trace {} to {}'.format(image_id,
                                                                       trace_id,
                                                                       new_swasp_id))
                            # split out the aperture we want
                            if args.split_apply:
                                scopy(iraf, image_id, trace_id, new_image_prefix)
                                # update the row for that original_image, set analyse = -1
                                if not original_row_updated:
                                    updateOldBlendTableRow(image_id)
                                    original_row_updated = True

                                # copy the new file to the right folder
                                split_spectrum = "{}.000{}.fits".format(new_image_prefix,
                                                                        trace_id)
                                # make a new swasp_id directory if it doesn't exist
                                if not os.path.exists('{}/{}'.format(all_spec_dir,
                                                                     new_swasp_id)):
                                    os.mkdir('{}/{}'.format(all_spec_dir,
                                                            new_swasp_id))
                                # if the new swasp_id is not the current one,
                                # copy it to the right place and then copy all
                                # splits to the blend common area
                                mvcp_str = "{}.*{}.fits".format(new_image_prefix,
                                                                trace_id)
                                if swasp_id != new_swasp_id:
                                    os.system('cp {} {}/{}'.format(split_spectrum,
                                                                   all_spec_dir,
                                                                   new_swasp_id))
                                    # move the split files to the common blended folder
                                    os.system('mv {} {}'.format(mvcp_str, all_blends_dir))
                                # if the new and old swasp_ids are the same,
                                # just copy the split file to the blend area,
                                # leave it where it is and copy the other splits
                                # out too
                                else:
                                    print('Leaving split file in current folder')
                                    templist = g.glob(mvcp_str)
                                    for tempimg in templist:
                                        if tempimg == split_spectrum:
                                            os.system('cp {} {}'.format(split_spectrum,
                                                                        all_blends_dir))
                                        else:
                                            os.system('mv {} {}'.format(tempimg,
                                                                        all_blends_dir))
                                # add rows to the table for each new trace file
                                addNewUnblendedTableRow(split_spectrum,
                                                        new_swasp_id,
                                                        old_row_info)

                    # move the original file into a common blended folder
                    # but only after all traces have been extracted
                    os.system('mv {} {}'.format(image_id, all_blends_dir))

