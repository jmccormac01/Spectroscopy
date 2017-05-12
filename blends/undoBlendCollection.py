"""
Run this script to undo the file changes

Also do the following:
    Reset analyse = 1 where analyse = -1 in db
    Remove rows that have image_ids like:
        .0001.fits, .0002.fits and .0003.fits
"""

import os
from contextlib import contextmanager
import glob as g
import pymysql

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
if __name__ == "__main__":
    print('Returning multispec files to original folders')
    os.chdir('/Users/jmcc/Dropbox/data/int/ids/eblm/all_spectra/all_blends')
    t = g.glob('*_tn.ms.fits')
    for i in t:
        qry = """
            SELECT swasp_id
            FROM eblm_ids_newest
            WHERE image_id = %s
            """
        with openDb() as cur:
            cur.execute(qry, (i, ))
            res = cur.fetchone()
            if res is None:
                print('No match for {}, exiting!'.format(i))
            else:
                print('{} --> {}'.format(i, res[0]))
                os.system('mv {} ../{}'.format(i, res[0]))

    # now look for *.000X.fits files in any folder and remove them
    print('Looking for lurking split files, get em!')
    os.chdir('/Users/jmcc/Dropbox/data/int/ids/eblm/all_spectra')
    t = g.glob('1SWASP*')
    for i in t:
        with cd(i):
            t2 = g.glob('*.000*.fits')
            if len(t2) > 0:
                print(i)
