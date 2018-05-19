"""
Script to extract and plot the CCFs from HARPS DRS

This script also makes a master CCF for meaursing the width
when fitting SB2s etc
"""
import os
import glob as g
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.interpolate import griddata
from publications.plots.defaults import (
    general,
    one_column
    )

if __name__ == "__main__":
    indir = '/Users/jmcc/Dropbox/NGTS/FullInstrument/Planets/NGTS-4b/HARPS/'
    ccfs = sorted(g.glob('{}/*_ccf_*.fits'.format(indir)))
    # phases from NGTS-4 paper notebook
    phases = np.array([0.14212619, 0.2629511, 0.35470913,
                       0.71084067, 0.34840988, 0.63096492,
                       0.90576128, 0.2404317,  0.48055739,
                       0.74425731])
    # from harps reduction with K5 mask
    rvs = np.array([38.71133, 38.48704, 38.54511, 39.62744,
                    38.56517, 39.47586, 39.47193, 38.51532,
                    38.93692, 39.67866])
    systemic = np.average(rvs)
    temp = zip(phases, rvs, ccfs)
    temp = sorted(temp)
    phases, rvs, ccfs = zip(*temp)

    # set up the figure
    general()
    one_column()
    fig, ax = plt.subplots(2, sharex=True)

    # interp grid - used for stacking ccfs to make a master
    interp_x = np.arange(0, 80, 0.25)
    data_aligned = []

    for i, ccf in enumerate(ccfs):
        ccf_filename = "{}/{}.ccf".format(indir, ccf.split('/')[-1].split('.f')[0])
        h = fits.open(ccf)
        delt = h[0].header['CDELT1']
        val = h[0].header['CRVAL1']
        x = np.arange(val, val+len(h[0].data[0])*delt, delt)
        data = np.sum(h[0].data, axis=0)
        data_norm = data / np.average(data[:120])
        # save the ccf to disc if not already there
        np.savetxt(ccf_filename, np.c_[x, data_norm], fmt='%.5f  %.5f', header='velocity_kms  contrast')
        ax[0].plot(x, data_norm + i*0.15, 'k-', lw=0.2)
        ax[0].text(10, 1+i*0.15, "phase={:.3f}".format(phases[i]), fontsize=6)

        # work out the difference needed to stack ccfs into master ccf
        diff = systemic - rvs[i]
        x_aligned = x + diff
        data_aligned.append(np.interp(interp_x, x_aligned, data))

    # stack the ccfs into the master and save it
    data_aligned = np.array(data_aligned)
    stacked = np.vstack(data_aligned)
    stacked_ccf = np.sum(stacked, axis=0)
    stacked_ccf_norm = stacked_ccf / np.average(stacked_ccf[:120])
    master_ccf_filename = "{}/HARPS_Master_10spectra.mccf".format(indir)
    if not os.path.exists(master_ccf_filename):
        np.savetxt(master_ccf_filename,
                   np.c_[interp_x, stacked_ccf_norm],
                   fmt='%.5f  %.5f',
                   header='velocity_kms  contrast')


    #ax[1].plot(interp_x, data_aligned + i*0.15, 'r-', lw=0.2)
    ax[0].axvline(systemic, lw=0.2)
    ax[0].set_ylabel('Contrast')

    # stacked ccf
    ax[1].plot(interp_x, stacked_ccf_norm, 'r-', lw=0.2)
    ax[1].axvline(systemic, lw=0.2)
    ax[1].set_xlabel('Radial velocity (km/s)')
    ax[1].set_ylabel('Contrast')
    fig.subplots_adjust(left=0.15, right=0.97, top=0.96,
                                bottom=0.12, hspace=0)
    fig.savefig('/Users/jmcc/Dropbox/NGTS/FullInstrument/Planets/NGTS-4b/HARPS/all_HARPS_CCFs.png')
