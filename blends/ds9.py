"""
Functions for handling the DS9 tasks in the NITES pipeline
"""
import os
import re
import time
import platform
from collections import defaultdict
import numpy as np
import pyregion

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

# pylint: disable = invalid-name
# pylint: disable = line-too-long
# pylint: disable = too-many-arguments
# pylint: disable = too-many-locals

def checkForDs9(ds9_name):
    """
    Find DS9 window with ds9_name if it exists

    Parameters
    ----------
    ds9_name : str
        ID of the DS9 window to check for

    Returns
    -------
    ds9 : bool
        Boolean showing if DS9 was found

    Raises
    ------
    None
    """
    ds9 = False
    if platform.system().lower() == 'darwin':
        processes = os.popen('ps aux | grep {}'.format(ds9_name)).readlines()
    else:
        processes = os.popen('ps aux --cols 1024 | grep {}'.format(ds9_name)).readlines()
    for process in processes:
        if "ds9" in process and " -title {}".format(ds9_name) in process:
            print('Hurray! Found DS9 window process')
            print('Waiting 20s to be sure ds9 is open...')
            time.sleep(20)
            ds9 = True
            break
    else:
        print('Boo! Where the hell is the DS9 window?')
    return ds9

def startDs9(ds9_name):
    """
    Start up a DS9 window. Doing it like this is much
    more robust than using pyds9.

    Parameters
    ----------
    ds9_name : str
        ID of the DS9 window to open

    Returns
    -------
    None

    Raises
    ------
    None
    """
    os.system('ds9 -title {} &'.format(ds9_name))

def ds9Display(ds9_name, image_id):
    """
    Display an image in DS9

    Parameters
    ----------
    ds9_name : str
        ID of the ds9 window to display an image in
    image_id : str
        Name of the FITS image to display

    Returns
    -------
    None

    Raises
    ------
    None
    """
    os.system('xpaset {} fits < {}'.format(ds9_name, image_id))

def ds9Set(ds9_name, command):
    """
    Set command for updating a ds9 window

    Parameters
    ----------
    ds9_name : str
        ID of the DS9 window to update
    command : str
        xpaset command for updating DS9 window
        see http://ds9.si.edu/doc/ref/xpa.html for more details
    Returns
    -------
    None

    Raises
    ------
    None
    """
    os.system('xpaset -p {} {}'.format(ds9_name, command))

def setupDs9(ds9_name, filename=None):
    """
    Configure the DS9 image viewing window

    Parameters
    ----------
    ds9_name : str
        ID of the DS9 window to set up
    filename : str
        Name of the file to display in DS9

    Returns
    -------
    None

    Raises
    ------
    None
    """
    print('Checking for DS9...')
    if not checkForDs9(ds9_name):
        print('No DS9 found, starting DS9...')
        startDs9(ds9_name)
        while not checkForDs9(ds9_name):
            print('Checking for DS9...')
            time.sleep(1)

    time.sleep(5)
    ds9Set(ds9_name, 'scale mode zscale')
    ds9Set(ds9_name, 'preserve scale yes')
    ds9Set(ds9_name, 'preserve pan yes')
    if filename:
        ds9Display(ds9_name, filename)
        ds9Set(ds9_name, 'zoom to fit')
        region_file = "{0:s}.reg".format(filename.split('.')[0])
        if os.path.exists(region_file):
            ds9Set(ds9_name, 'regions < {0:s}'.format(region_file))

