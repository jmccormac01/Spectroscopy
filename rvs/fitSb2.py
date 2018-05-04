"""
Code to fit double gaussians to SB2 CCFs

Need to generalise this code for priors and the inputs
"""
import sys
import argparse as ap
import pickle
from collections import OrderedDict
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np
import emcee
import corner
from astropy.modeling.models import Voigt1D

# pylint: disable=invalid-name
# pylint: disable=superfluous-parens
# pylint: disable=redefined-outer-name

def argParse():
    """
    Parse the command line arguments
    """
    p = ap.ArgumentParser()
    p.add_argument('--plot_ccf_only',
                   help='used to plot CCF only (to estimate starting position',
                   action='store_true')
    p.add_argument('--norm',
                   help='normalise the flat part of the CCF to 1?',
                   action='store_true')
    return p.parse_args()

def readInitialAndPriors(infile):
    """
    Read in the initial guess parameters
    and their associated uniform priors
    """
    init_priors = OrderedDict()
    f = open(infile).readlines()
    for line in f:
        print(line.split())
        if line.split()[0] == 'norm_range':
            name, lower, upper = line.split()
            init_priors[name] = (float(lower), float(upper))
        else:
            name, initial, sigma_walkers, prior_l, prior_h = line.split()
            init_priors[name] = OrderedDict()
            init_priors[name]['initial'] = float(initial)
            init_priors[name]['sigma_walkers'] = float(sigma_walkers)
            init_priors[name]['prior_l'] = float(prior_l)
            init_priors[name]['prior_h'] = float(prior_h)
        if 'norm_range' not in init_priors.keys():
            init_priors['norm_range'] = None
    return init_priors

def lnprior(theta, init_priors, voigt=False):
    """
    """
    # put very wide liberal uniform priors on all variables
    if voigt:
        m1, s1g, s1l, a1, m2, s2g, s2l, a2 = theta
        if init_priors['m1']['prior_l'] < m1 < init_priors['m1']['prior_h'] and \
            init_priors['s1g']['prior_l'] < s1g < init_priors['s1g']['prior_h'] and \
            init_priors['s1l']['prior_l'] < s1l < init_priors['s1l']['prior_h'] and \
            init_priors['a1']['prior_l'] < a1 < init_priors['a1']['prior_h'] and \
            init_priors['m2']['prior_l'] < m2 < init_priors['m2']['prior_h'] and \
            init_priors['s2g']['prior_l'] < s2g < init_priors['s2g']['prior_h'] and \
            init_priors['s2l']['prior_l'] < s2l < init_priors['s2l']['prior_h'] and \
            init_priors['a2']['prior_l'] < a2 < init_priors['a2']['prior_h']:
            return 0.0
        else:
            return -np.inf
    else:
        m1, s1, a1, m2, s2, a2 = theta
        if init_priors['m1']['prior_l'] < m1 < init_priors['m1']['prior_h'] and \
            init_priors['s1']['prior_l'] < s1 < init_priors['s1']['prior_h'] and \
            init_priors['a1']['prior_l'] < a1 < init_priors['a1']['prior_h'] and \
            init_priors['m2']['prior_l'] < m2 < init_priors['m2']['prior_h'] and \
            init_priors['s2']['prior_l'] < s2 < init_priors['s2']['prior_h'] and \
            init_priors['a2']['prior_l'] < a2 < init_priors['a2']['prior_h']:
            return 0.0
        else:
            return -np.inf

def lnlike(theta, x, y, yerr, voigt=False):
    """
    """
    # voigt or gaussian?
    if voigt:
        m1, s1g, s1l, a1, m2, s2g, s2l, a2 = theta
        model = doubleVoigt(x, m1, s1g, s1l, a1, m2, s2g, s2l, a2)
    else:
        m1, s1, a1, m2, s2, a2 = theta
        model = doubleGaussian(x, m1, s1, a1, m2, s2, a2)

    if True in np.isnan(model) or np.min(model) <= 0:
        lnlike = -np.inf
    else:
        inv_sigma2 = 1.0/(yerr**2)
        eq_p1 = (y-model)**2*inv_sigma2 - np.log(inv_sigma2)
        lnlike = -0.5*(np.sum(eq_p1) - np.log(len(y) + 1))
    return lnlike

def lnprob(theta, init_priors, x, y, yerr, voigt=False):
    """
    """
    lp = lnprior(theta, init_priors, voigt=voigt)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr, voigt=voigt)

def doubleGaussian(x, m1, s1, a1, m2, s2, a2):
    """
    Model for a double Gaussian profile
    """
    # primary peak
    g1 = np.exp(-0.5*((x-m1)/s1)**2)
    # secondary peak
    g2 = np.exp(-0.5*((x-m2)/s2)**2)
    # total model
    mod1 = 1 - a1 * g1
    mod2 = 1 - a2 * g2
    modt = mod1 + mod2 - 1
    return modt

def doubleVoigt(x, m1, s1g, s1l, a1, m2, s2g, s2l, a2):
    """
    Model for a double Voigt profile
    """
    v1 = Voigt1D(x_0=m1, amplitude_L=a1, fwhm_L=s1l, fwhm_G=s1g)
    v2 = Voigt1D(x_0=m2, amplitude_L=a2, fwhm_L=s2l, fwhm_G=s2g)
    mod1 = 1 - v1(x)
    mod2 = 1 - v2(x)
    modt = mod1 + mod2 - 1
    return modt

    # plot the final model and output the results
if __name__ == "__main__":
    args = argParse()
    # read in the CCF to fit
    data_dir = "/Users/jmcc/Dropbox/NGTS/FullInstrument/Planets/NGTS-4b/HARPS"
    data_file_base = "HARPS.2017-12-23T03_53_23.073_s1d_A"
    data_file = "{}/{}.ccf".format(data_dir, data_file_base)
    init_priors_file = "{}/{}.sb2".format(data_dir, data_file_base)
    # used for CCF text files
    x, y = np.loadtxt(data_file, usecols=[0, 1], unpack=True)
    yerr = np.array([0.005]*len(y))
    # get initial parameters and priors
    init_priors = readInitialAndPriors(init_priors_file)

    if args.plot_ccf_only:
        fig, ax = plt.subplots(1, figsize=(10, 10))
        ax.plot(x, y, 'k-')
        ax.set_ylabel('Contrast')
        ax.set_xlabel('Radial Velocity (km/s)')
        if not args.norm:
            plt.show()
            sys.exit()

    # normalise the flat part of CCF
    if args.norm:
        if init_priors['norm_range']:
            n_range = np.where(((x <= init_priors['norm_range'][0]) | (x >= init_priors['norm_range'][1])))
            coeffs = np.polyfit(x[n_range], y[n_range], 1)
            besty = np.polyval(coeffs, x)
            y = y / besty
            if args.plot_ccf_only:
                ax.plot(x, y, 'r-')
                plt.show()
                sys.exit()
        else:
            print('Normalising requested but no norm_range specified')

    # used for CERES CCFs
    #data = pickle.load(open(data_file, 'r'))
    #x = data['rvels']
    #y = data['rxc_av']
    #yerr = np.array([0.01]*len(y))

    # determine if using gaussian or voigt profile?
    if 's1g' in init_priors.keys() and 's1l' in init_priors.keys():
        print('Using Voigt profile')
        FIT_TYPE = 'voigt'
        # set up the MCMC arrays
        initial = [init_priors['m1']['initial'],
                   init_priors['s1g']['initial'],
                   init_priors['s1l']['initial'],
                   init_priors['a1']['initial'],
                   init_priors['m2']['initial'],
                   init_priors['s2g']['initial'],
                   init_priors['s2l']['initial'],
                   init_priors['a2']['initial']]
        parameters = [key for key in init_priors if key != 'norm_range']
        weights = [init_priors['m1']['sigma_walkers'],
                   init_priors['s1g']['sigma_walkers'],
                   init_priors['s1l']['sigma_walkers'],
                   init_priors['a1']['sigma_walkers'],
                   init_priors['m2']['sigma_walkers'],
                   init_priors['s2g']['sigma_walkers'],
                   init_priors['s2l']['sigma_walkers'],
                   init_priors['a2']['sigma_walkers']]
    else:
        print('Using Gaussian profile')
        FIT_TYPE = 'gaussian'
        # set up the MCMC arrays
        initial = [init_priors['m1']['initial'],
                   init_priors['s1']['initial'],
                   init_priors['a1']['initial'],
                   init_priors['m2']['initial'],
                   init_priors['s2']['initial'],
                   init_priors['a2']['initial']]
        parameters = [key for key in init_priors if key != 'norm_range']
        weights = [init_priors['m1']['sigma_walkers'],
                   init_priors['s1']['sigma_walkers'],
                   init_priors['a1']['sigma_walkers'],
                   init_priors['m2']['sigma_walkers'],
                   init_priors['s2']['sigma_walkers'],
                   init_priors['a2']['sigma_walkers']]
    # set up sampler
    ndim = len(initial)
    nwalkers = 4*ndim*2
    nsteps = 5000
    pos = [initial + weights*np.random.randn(ndim) for i in range(nwalkers)]
    # sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                                    args=(init_priors, x, y, yerr, FIT_TYPE),
                                    threads=6)
    print("Running MCMC...")
    # run sampling with progress status
    #sampler.run_mcmc(pos, nsteps, rstate0=np.random.get_state())
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps,
                                              rstate0=np.random.get_state())):
        if (i+1) % 100 == 0:
            print("{0:5.1%}".format(float(i) / nsteps))
    print("Done.")
    # plot and save the times series of each parameter
    for i, (initial_param, label) in enumerate(zip(initial, parameters)):
        fig, ax = plt.subplots(1, figsize=(10, 10))
        ax.plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
        ax.axhline(initial_param, color="#888888", lw=2)
        ax.set_ylabel(label)
        ax.set_xlabel('step number')
        fig.savefig('{}/chain_{}steps_{}walkers_{}.png'.format(data_dir,
                                                               nsteps,
                                                               nwalkers,
                                                               label))
    burnin = int(raw_input('Enter burnin period: '))
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))

    if FIT_TYPE == 'voigt':
        # extract the results from post-burnin
        m1 = np.median(samples[:, 0])
        s1g = np.median(samples[:, 1])
        s1l = np.median(samples[:, 2])
        a1 = np.median(samples[:, 3])
        m2 = np.median(samples[:, 4])
        s2g = np.median(samples[:, 5])
        s2l = np.median(samples[:, 6])
        a2 = np.median(samples[:, 7])
        em1 = np.std(samples[:, 0])
        es1g = np.std(samples[:, 1])
        es1l = np.std(samples[:, 2])
        ea1 = np.std(samples[:, 3])
        em2 = np.std(samples[:, 4])
        es2g = np.std(samples[:, 5])
        es2l = np.std(samples[:, 6])
        ea2 = np.std(samples[:, 7])

        # print the results
        print('m1 = {} +/- {}'.format(m1, em1))
        print('s1g = {} +/- {}'.format(s1g, es1g))
        print('s1l = {} +/- {}'.format(s1l, es1l))
        print('a1 = {} +/- {}'.format(a1, ea1))
        print('m2 = {} +/- {}'.format(m2, em2))
        print('s2g = {} +/- {}'.format(s2g, es2g))
        print('s2l = {} +/- {}'.format(s2l, es2l))
        print('a2 = {} +/- {}'.format(a2, ea2))

        # make a corner plot
        fig = corner.corner(samples,
                            labels=["$m1$",
                                    "$s1g$",
                                    "$s1l$",
                                    "$a1$",
                                    "$m2$",
                                    "$s2g$",
                                    "$s2l$",
                                    "$a2$"],
                            truths=initial,
                            plot_contours=False)
        fig.savefig('{}/corner_{}steps_{}walkers.png'.format(data_dir, nsteps, nwalkers))
    else:
        # extract the results from post-burnin
        m1 = np.median(samples[:, 0])
        s1 = np.median(samples[:, 1])
        a1 = np.median(samples[:, 2])
        m2 = np.median(samples[:, 3])
        s2 = np.median(samples[:, 4])
        a2 = np.median(samples[:, 5])
        em1 = np.std(samples[:, 0])
        es1 = np.std(samples[:, 1])
        ea1 = np.std(samples[:, 2])
        em2 = np.std(samples[:, 3])
        es2 = np.std(samples[:, 4])
        ea2 = np.std(samples[:, 5])

        # print the results
        print('m1 = {} +/- {}'.format(m1, em1))
        print('s1 = {} +/- {}'.format(s1, es1))
        print('a1 = {} +/- {}'.format(a1, ea1))
        print('m2 = {} +/- {}'.format(m2, em2))
        print('s2 = {} +/- {}'.format(s2, es2))
        print('a2 = {} +/- {}'.format(a2, ea2))

        # make a corner plot
        fig = corner.corner(samples,
                            labels=["$m1$",
                                    "$s1$",
                                    "$a1$",
                                    "$m2$",
                                    "$s2$",
                                    "$a2$"],
                            truths=initial,
                            plot_contours=False)
        fig.savefig('{}/corner_{}steps_{}walkers.png'.format(data_dir, nsteps, nwalkers))

    # plot the final model and output the results
    fig2, ax2 = plt.subplots(2, figsize=(10, 10))
    # final model
    if FIT_TYPE == 'voigt':
        final_model = doubleVoigt(x, m1, s1g, s1l, a1, m2, s2g, s2l, a2)
        # plot the final model
        ax2[0].plot(x, y, 'k-')
        ax2[0].plot(x, final_model, 'r-')
        ax2[0].legend(('CCF', 'SB2tool model'), loc=2)
        ax2[0].set_ylabel('CCF')
        ax2[0].set_xlabel('Radial velocity (km/s)')
        ax2[1].plot(x, y, 'k-')
        ax2[1].plot(x, final_model, 'r-')
        ax2[1].set_xlim(m1-2*s1g, m1+2*s1g)
        ax2[1].set_ylim(1-a1-0.05, 1-a1+0.15)
        ax2[1].legend(('CCF', 'SB2tool model'), loc=2)
        ax2[0].set_ylabel('CCF')
        ax2[1].set_xlabel('Radial velocity (km/s)')
        fig2.savefig('{}/fitted_double_voigt.png'.format(data_dir))
        # save the results to file
        with open('{}/fitted_double_voigt.log'.format(data_dir), 'w') as of:
            of.write('m1 = {} +/- {}\n'.format(m1, em1))
            of.write('s1g = {} +/- {}\n'.format(s1g, es1g))
            of.write('s1l = {} +/- {}\n'.format(s1l, es1l))
            of.write('a1 = {} +/- {}\n'.format(a1, ea1))
            of.write('m2 = {} +/- {}\n'.format(m2, em2))
            of.write('s2g = {} +/- {}\n'.format(s2g, es2g))
            of.write('s2l = {} +/- {}\n'.format(s2l, es2l))
            of.write('a2 = {} +/- {}\n'.format(a2, ea2))
    else:
        final_model = doubleGaussian(x, m1, s1, a1, m2, s2, a2)
        # plot the final model
        ax2[0].plot(x, y, 'k-')
        ax2[0].plot(x, final_model, 'r-')
        ax2[0].legend(('CCF', 'SB2tool model'), loc=2)
        ax2[0].set_ylabel('CCF')
        ax2[0].set_xlabel('Radial velocity (km/s)')
        ax2[1].plot(x, y, 'k-')
        ax2[1].plot(x, final_model, 'r-')
        ax2[1].set_xlim(m1-2*s1, m1+2*s1)
        ax2[1].set_ylim(1-a1-0.05, 1-a1+0.15)
        ax2[1].legend(('CCF', 'SB2tool model'), loc=2)
        ax2[0].set_ylabel('CCF')
        ax2[1].set_xlabel('Radial velocity (km/s)')
        fig2.savefig('{}/fitted_double_gaussian.png'.format(data_dir))
        # save the results to file
        with open('{}/fitted_double_gaussian.log'.format(data_dir), 'w') as of:
            of.write('m1 = {} +/- {}\n'.format(m1, em1))
            of.write('s1 = {} +/- {}\n'.format(s1, es1))
            of.write('a1 = {} +/- {}\n'.format(a1, ea1))
            of.write('m2 = {} +/- {}\n'.format(m2, em2))
            of.write('s2 = {} +/- {}\n'.format(s2, es2))
            of.write('a2 = {} +/- {}\n'.format(a2, ea2))
    plt.show()

