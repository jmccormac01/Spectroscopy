"""
Here were check the database for all spectra for 
each object and check they are all in one place 
for analysis. If they are do nothing, if they are not
then copy them to the iSpec folder. 
"""

import os
import argparse as ap
from datetime import timedelta
import pymysql

# set up the per object master directory
top_dir = '/Users/James/Dropbox/data/int/ids/eblm'
iSpec_dir = '{}/all_spectra'.format(top_dir)
if not os.path.exists(iSpec_dir):
    os.mkdir(iSpec_dir)

# connect to the database
db = pymysql.connect(host='localhost',
                     db='eblm',
                     password='mysqlpassword')

def argParse():
    """
    Parse the command-line arguments
    """
    description = "A script to gather up spectra on per object basis"
    parser = ap.ArgumentParser(description=description)
    parser.add_argument('--copy',
                        help='copy the spectra',
                        action='store_true')
    return parser.parse_args()

def getUniqueTargetList():
    """
    Gather up a list of unique eblm targets in database
    """
    qry = """
        SELECT distinct(swasp_id)
        FROM eblm_ids_final
        WHERE swasp_id IS NOT NULL
        """
    swasp_ids = []
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            if row[0] != 'None':
                swasp_ids.append(row[0])
    return swasp_ids

def gatherAllTargetSpectra(swasp_id):
    """
    Gather all the spectra for a particular target
    in one place for analysis with iSpec.

    If the spectrum exists in the target location
    do not copy it again
    """
    # keep a count of the number of spectra to copy
    total = 0
    # check for the directory that is to contain the spectra
    target_dir = '{}/{}'.format(iSpec_dir, swasp_id)
    if not os.path.exists(target_dir) and args.copy:
        os.mkdir(target_dir)
    # qry all the spectra in the database for this object
    qry = """
        SELECT image_id,night
        FROM eblm_ids_final
        WHERE swasp_id = '{}'
        """.format(swasp_id)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            # look in the destination directory for each spectrum
            # if it exists continue, if not copy it to the destination
            image_id = row[0]
            night = row[1].strftime('%Y%m%d')
            spectrum_source = '{}/{}/reduced/spectra_r/{}'.format(top_dir, night, image_id)
            spectrum_dest = '{}/{}/{}'.format(iSpec_dir, swasp_id, image_id)
            if not os.path.exists(spectrum_dest):
                if os.path.exists(spectrum_source):
                    print('Copying {} to {}'.format(spectrum_source, spectrum_dest))
                    if args.copy:
                        os.system('cp {} {}/{}/'.format(spectrum_source, iSpec_dir, swasp_id))
                    total += 1
                else:
                    print('{} DOES NOT EXIST!'.format(spectrum_source))
            else:
                print('{} EXISTS IN COMBINED FOLDER'.format(spectrum_dest))
    print('\n\n')
    return total

if __name__ == '__main__':
    args = argParse()
    swasp_ids = getUniqueTargetList()
    total = 0
    for swasp_id in sorted(swasp_ids):
        subtotal = gatherAllTargetSpectra(swasp_id)
        total += subtotal
    print('{} spectra to be copied'.format(total))
