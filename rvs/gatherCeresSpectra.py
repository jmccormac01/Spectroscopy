"""
Script to gather up the NOT/FIES spectra reduced by CERES
from CAFE/FIES into per object folders

This is useful to check the CCFs etc
"""
import os
import argparse as ap
import glob as g

# pylint: disable=invalid-name
# pylint: disable=superfluous-parens

def argParse():
    """
    Parse the command line
    """
    p = ap.ArgumentParser()
    p.add_argument('instrument',
                   help='name of instrument to gather spectra from',
                   choices=['fies', 'cafe'])
    return p.parse_args()

if __name__ == "__main__":
    args = argParse()
    top_dir = '/Users/jmcc/Dropbox/data/{}'.format(args.instrument)
    all_spec_dir = '{}/all_spectra'.format(top_dir)

    # get a list of reduced nights
    os.chdir(top_dir)
    dir_list = g.glob('*_red*')

    for directory in dir_list:
        os.chdir('{}/proc'.format(directory))
        f = open('results.txt').readlines()
        swasp_ids = []
        for line in f:
            s = line.split()
            swasp_id = s[0]
            if swasp_id not in swasp_ids:
                swasp_ids.append(swasp_id)
                if not os.path.exists('{}/{}'.format(all_spec_dir, swasp_id)):
                    os.mkdir('{}/{}'.format(all_spec_dir, swasp_id))
            pdf = s[-1]
            # copy the pdfs
            print('{} {}'.format(swasp_id, pdf))
            if not os.path.exists('{}/{}/{}'.format(all_spec_dir, swasp_id, pdf)):
                comm = 'cp {} {}/{}'.format(pdf, all_spec_dir, swasp_id)
                print(comm)
                os.system(comm)
            else:
                print('{}/{}/{} exists...'.format(all_spec_dir, swasp_id, pdf))
        # copy the fits spectra
        for swasp_id in swasp_ids:
            spectra = g.glob('{}*.fits'.format(swasp_id))
            for spectrum in spectra:
                if not os.path.exists('{}/{}/{}'.format(all_spec_dir, swasp_id, spectrum)):
                    comm2 = 'cp {} {}/{}'.format(spectrum, all_spec_dir, swasp_id)
                    print(comm2)
                    os.system(comm2)
                else:
                    print('{}/{}/{} exists'.format(all_spec_dir, swasp_id, spectrum))
        os.chdir('../../')
