"""
Script to accurately measure CCF bisector spans

Take in a CCF
Apply lower and upper contrast cuts
Split it in half about the mid-point
Step through in even steps in contrast
Find the nearest two points
Fit a line between them and solve for the veloicty at the contrast step
Store contrast step and velocity in dictionary
Repeat for the other side of the CCF
Loop over the contrasts/velocities and find the mid point velocity
Plot these velocities
[Maybe fit some sort of function to them too]
"""
import os
import glob as g
import numpy as np
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt

def bootstrap(x, y):
    """
    Resample with replacement a bunch of times
    to get errors on the bisector measurements
    and also the final bisector slope
    """
    print('Determining the error on the gradient via bootstrap')
    stack = np.vstack((x, y)).transpose()
    gradients, intersections = [], []
    for i in range(0, 1000):
        samples = stack[np.random.choice(stack.shape[0], stack.shape[0], replace=True), :]
        x_sample = samples[:, 0]
        y_sample = samples[:, 1]
        coeffs_sample = np.polyfit(x_sample, y_sample, 1)
        gradients.append(coeffs_sample[0])
        intersections.append(coeffs_sample[1])
    slope_mean = np.average(gradients)
    slope_error = np.std(np.array(gradients))
    intersection_mean = np.average(intersections)
    intersection_error = np.std(np.array(intersections))
    return slope_mean, slope_error, intersection_mean, intersection_error

# load a test CCF
ccf_dir = "/Users/jmcc/Dropbox/NGTS/FullInstrument/Planets/NGTS-4b/HARPS"
os.chdir(ccf_dir)
ccfs = sorted(g.glob('*.ccf'))

base_scale = 0.10
tip_scale = 0.95
# from harps reduction with K5 mask
rvs = np.array([38.71133, 38.48704, 38.54511, 39.62744,
                38.56517, 39.47586, 39.47193, 38.51532,
                38.93692, 39.67866])

# make the cuts to get only the wings between the contrast limits
bisector_slope = np.empty(len(ccfs))
bisector_slope_error = np.empty(len(ccfs))
for i, ccf in enumerate(ccfs):
    velocity, contrast = np.loadtxt(ccf, usecols=[0, 1], unpack=True)
    # estimate the range to measure from the contrast
    cmin = min(contrast)
    cheight = 1 - cmin
    cstep = cheight/100
    cbase = 1 - base_scale*cheight
    ctip = 1 - tip_scale*cheight

    # mask out the CCF within the fitting range
    n = np.where(((contrast >= ctip) & (contrast <= cbase) & (velocity > 25) & (velocity < 50)))[0]
    velocity_fit = np.copy(velocity[n])
    contrast_fit = np.copy(contrast[n])
    # now split the velocities/contrasts about the RV point to get the two wings
    nleft = np.where(velocity_fit < rvs[i])[0]
    nright = np.where(velocity_fit > rvs[i])[0]

    velocity_fit_left = np.copy(velocity_fit[nleft])
    velocity_fit_right = np.copy(velocity_fit[nright])
    contrast_fit_left = np.copy(contrast_fit[nleft])
    contrast_fit_right = np.copy(contrast_fit[nright])

    # range of contrasts to fit for the bisectors
    contrast_range = np.arange(ctip, cbase, cstep)
    contrast_range_left = contrast_range[::-1]
    contrast_range_right = np.copy(contrast_range)
    # resample the left wing to fit the bisectors
    velocity_fit_left_resampled = np.interp(contrast_range_left[::-1],
                                            contrast_fit_left[::-1],
                                            velocity_fit_left[::-1])
    # resample the right wing to fit the bisectors
    velocity_fit_right_resampled = np.interp(contrast_range_right,
                                             contrast_fit_right,
                                             velocity_fit_right)

    # plot the CCF
    # set up the plots
    fig, ax = plt.subplots(2, figsize=(10, 10))
    ax[0].set_xlabel('Radial velocity (km/s)')
    ax[0].set_ylabel('Contrast')
    ax[1].set_ylabel('Radial velocity (km/s)')
    ax[1].set_xlabel('Contrast')

    # plot the profile and highlight the wings
    ax[0].plot(velocity, contrast, 'k--', lw=1.0)
    ax[0].plot(velocity_fit_left, contrast_fit_left, 'r-', lw=1.0)
    ax[0].plot(velocity_fit_right, contrast_fit_right, 'b-', lw=1.0)
    ax[0].axhline(cbase, lw=0.5, color='black', ls='dashed')
    ax[0].axhline(ctip, lw=0.5, color='black', ls='dashed')

    # plot the profile intersection points
    ax[0].plot(velocity_fit_left_resampled, contrast_range_left[::-1], 'ko', ms=0.5)
    ax[0].plot(velocity_fit_right_resampled, contrast_range_right, 'ko', ms=0.5)

    # plot the profile intersection lines
    bisector_velocity = np.empty(len(velocity_fit_left_resampled))
    for j in range(0, len(velocity_fit_left_resampled)):
        x = [velocity_fit_left_resampled[j], velocity_fit_right_resampled[j]]
        y = [contrast_range_left[::-1][j], contrast_range_right[j]]
        bisector_velocity[j] = x[0] + ((x[1] - x[0])/2.)
        ax[0].plot(x, y, 'k-', lw=0.5)

    # plot the mid-points of the intersection points
    ax[0].plot(bisector_velocity, contrast_range, 'ko', ms=1)

    # fit the mid-points of the insersections and work out the slope
    slope_mean, slope_error, intersection_mean, intersection_error  = bootstrap(contrast_range,
                                                                                bisector_velocity)
    besty = np.polyval(np.array([slope_mean, intersection_mean]), contrast_range)
    ax[0].plot(besty, contrast_range, 'g-', lw=1.5)
    bisector_slope[i] = slope_mean
    bisector_slope_error[i] = slope_error

    # plot the fit
    ax[1].plot(contrast_range, bisector_velocity, 'ko', ms=2)
    ax[1].plot(contrast_range, besty, 'g-', lw=1.5)
    ax[1].legend(('Bisector midpoints', 'y={:.3f} ({:.3f})x + {:.3f} ({:.3f})'.format(slope_mean,
                                                                                      slope_error,
                                                                                      intersection_mean,
                                                                                      intersection_error)), )
    fig.subplots_adjust(bottom=0.08, top=0.95, left=0.10, right=0.95)
    fig.savefig('{}/{}_bisectors.png'.format(ccf_dir, ccf.split('.ccf')[0]), dpi=300)

# now plot the bisectors and fit them
bisector_trend_mean, bisector_trend_error, bisector_intersection, bisector_intersection_error = bootstrap(rvs, bisector_slope)
besty2 = np.polyval(np.array([bisector_trend_mean, bisector_intersection]), rvs)
fig2, ax2 = plt.subplots(1, figsize=(10, 10))
ax2.errorbar(rvs, bisector_slope, yerr=bisector_slope_error, fmt='ko')
ax2.plot(rvs, besty2, 'k-')
ax2.set_xlabel('Radial velocity (km/s)')
ax2.set_ylabel('Bisector slope')
ax2.legend(('Bisector slopes', 'y={:.3f} ({:.3f})x + {:.3f} ({:.3f})'.format(bisector_trend_mean,
                                                                             bisector_trend_error,
                                                                             bisector_intersection,
                                                                             bisector_intersection_error),), )
fig2.subplots_adjust(bottom=0.08, top=0.95, left=0.10, right=0.95)
fig2.savefig('{}/BISECTORS_JMCC.png'.format(ccf_dir), dpi=300)
plt.show()
