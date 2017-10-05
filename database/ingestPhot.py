"""
Script to ingest phot for interesting objects

Process:
    Check there is no ingested phot already
    If so, stop (remove later manuall?)
    If not, ingest it
"""

import glob as g
import argparse as ap
import pymysql
import numpy as np

# pylint: disable = superfluous-parens
# pylint: disable = invalid-name
# pylint: disable = redefined-outer-name

db = pymysql.connect(host='localhost',
                     db='eblm',
                     password='mysqlpassword')

def argParse():
    """
    Parse the command line argments
    """
    description = '''Script for ingesting photometry into database
                    Select the instrument e.g. WASP | NITES
                    Data should be in diff_mags for ingesting
                    You can use the options
                    swasp_id = 'find' and
                    infile = 'find'
                    to ingest many light curves of the same instrument set up
                    '''
    filters = ["U", "B", "V", "R", "I", "clear", "u", "g", "r", "i", "z"]
    parser = ap.ArgumentParser(description=description)
    parser.add_argument('swasp_id', help='swasp_id')
    parser.add_argument('infile', help='path/to/phot/file')
    parser.add_argument('instrument', choices=['wasp', 'nites'], help='instrument')
    parser.add_argument('filter', choices=filters, help='filter used for observations')
    parser.add_argument('--commit', help='filter used for observations', action='store_true')
    return parser.parse_args()


def checkForPhot(swasp_id, instrument):
    """
    Check if there is phot so we don't duplicate
    """
    print('Checking for phot for {} {}'.format(swasp_id, instrument))
    qry = """
        SELECT count(*)
        FROM eblm_phot
        WHERE swasp_id='{}' AND
        instrument='{}'
        """.format(swasp_id, instrument)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            if int(row[0]) < 1:
                print('No phot found...')
                return False
            else:
                print('PHOT FOUND!')
                return True

def ingest(t, m, e, swasp_id, instrument, filt):
    """
    Ingest the new photometry
    """
    print('Ingesting phot for {} {}'.format(swasp_id, instrument))
    if args.commit:
        with db.cursor() as cur:
            for i, j, k in zip(t, m, e):
                qry = """
                    INSERT INTO eblm_phot
                    (swasp_id, hjd_mid, mag, mag_err, instrument, filter)
                    VALUES
                    ('{}', {}, {}, {}, '{}', '{}')
                    """.format(swasp_id, i, j, k, instrument, filt)
                cur.execute(qry)
                db.commit()

if __name__ == "__main__":
    args = argParse()
    # ingest many of the same instrument setup
    if args.infile == 'find' and args.swasp_id == 'find':
        lis = g.glob('1SWASP*')
        for i in lis:
            swasp_id = i.split('_')[0]
            print swasp_id
            # check for phot already ingested
            res = checkForPhot(swasp_id, args.instrument)
            if not res:
                t, m, e = np.loadtxt(i, usecols=[0, 1, 2], unpack=True)
                ingest(t, m, e, swasp_id, args.instrument, args.filter)
                print("Ingested {}".format(i))
            else:
                print("{} {} PHOT ALREADY EXISTS, CHECK! Skipping...".format(swasp_id, args.instrument))
    else:
        # check for phot already ingested
        res = checkForPhot(args.swasp_id, args.instrument)
        if not res:
            t, m, e = np.loadtxt(args.infile, usecols=[0, 1, 2], unpack=True)
            ingest(t, m, e, args.swasp_id, args.instrument, args.filter)
            print("Ingested {}".format(args.infile))
        else:
            print("{} {} PHOT ALREADY EXISTS, CHECK!!".format(args.swasp_id, args.instrument))

