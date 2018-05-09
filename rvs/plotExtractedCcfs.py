"""
Script to plot extracted CCFs from iSpec,
colour coded based on orbital phase

The CCFs from iSpec are not as good as those from
HARPS DRS, which is done order by order. Use the
plotHarpsCcfs script instead!!!!
"""
import glob as g
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from publications.plots.defaults import (
    general,
    one_column
    )

if __name__ == "__main__":
    ccfs = sorted(g.glob('*.ccf'))
    general()
    one_column()
    fig, ax = plt.subplots(1)
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
    temp = zip(phases, ccfs)
    temp = sorted(temp)
    phases, ccfs = zip(*temp)
    for i, ccf in enumerate(ccfs):
        x, y = np.loadtxt(ccf, usecols=[0, 1], unpack=True)
        n_range = np.where(((x <= 0) | (x >= 100)))
        coeffs = np.polyfit(x[n_range], y[n_range], 1)
        besty = np.polyval(coeffs, x)
        y = y / besty
        ax.plot(x, y+i*0.15, 'k-', linewidth=0.2)
        ax.text(10, 1+i*0.15, "phase={:.3f}".format(phases[i]), fontsize=6)
    ax.set_xlabel('Radial velocity (km/s)')
    ax.set_ylabel('Contrast')
    #ax.set_ylim(0.55, 1.1)
    ax.set_xlim(0, 80)
    ax.axvline(systemic, lw=0.2)
    fig.subplots_adjust(left=0.15, right=0.97, top=0.96,
                                bottom=0.12, hspace=0)
    fig.savefig('/Users/jmcc/Dropbox/NGTS/FullInstrument/Planets/NGTS-4b/HARPS/all_CCFs_zoom.png')
