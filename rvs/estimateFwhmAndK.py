"""
Script to look through the SB1s and estimate K and FWHM_avg from available data

Run this script ocassionally to update the CCF average FWHM and K in the parameters table
"""
import argparse as ap
from collections import OrderedDict
import numpy as np
import pymysql

# pylint: disable=invalid-name
# pylint: disable=redefined-outer-name

table_map = {'fies': 'eblm_fies',
             'cafe': 'eblm_cafe',
             'ids': 'eblm_ids_newest'}

def argParse():
    """
    Parse the command line args
    """
    p = ap.ArgumentParser()
    p.add_argument('--swasp_id',
                   help='swasp_id to update')
    return p.parse_args()

def getListOfSb1s():
    """
    Get a list of SB1s to analyse
    """
    swasp_ids = []
    qry = """
        SELECT swasp_id
        FROM eblm_parameters
        WHERE current_status LIKE 'SB1%'
        OR current_status = 'CONTINUE'
        ORDER BY swasp_id ASC
        """
    with pymysql.connect(host='localhost', db='eblm',
                         user='jmcc', password='mysqlpassword') as cur:
        cur.execute(qry)
        results = cur.fetchall()
        for result in results:
            swasp_ids.append(result[0])
    return swasp_ids

def getFwhms(swasp_ids):
    """
    Make a dictionary of FWHMS forr each object from each spectrograph
    """
    fwhms = OrderedDict()
    with pymysql.connect(host='localhost', db='eblm',
                         user='jmcc', password='mysqlpassword') as cur:
        for swasp_id in swasp_ids:
            fwhms[swasp_id] = {}
            for instrument in table_map:
                print('{}: {}'.format(swasp_id, instrument))
                qry = """
                    SELECT mask_ccf_fwhm
                    FROM {}
                    WHERE swasp_id = '{}'
                    AND analyse = 1
                    AND mask_ccf_fwhm IS NOT NULL
                    """.format(table_map[instrument], swasp_id)
                cur.execute(qry)
                results = cur.fetchall()
                if len(results) > 0:
                    temp = np.empty(len(results))
                    for i, row in enumerate(results):
                        temp[i] = row[0]
                    K = np.average(temp)
                    fwhms[swasp_id][instrument] = round(K, 3)
                else:
                    fwhms[swasp_id][instrument] = None
    return fwhms

def getPhasedRvData(swasp_id, instrument):
    """
    Grab all RVs for a specific object from a specific instrument
    """
    if instrument == 'ids':
        qry = """
            SELECT bjd_mid, mask_velocity,
            barycentric_velocity_iSpec
            FROM {}
            WHERE swasp_id = '{}'
            AND analyse = 1
            """.format(table_map[instrument], swasp_id)
    else:
        qry = """
            SELECT bjd_mid, mask_velocity
            FROM {}
            WHERE swasp_id = '{}'
            AND analyse = 1
            """.format(table_map[instrument], swasp_id)
    qry2 = """
        SELECT epoch, period
        FROM eblm_parameters
        WHERE swasp_id = '{}'
        """.format(swasp_id)
    with pymysql.connect(host='localhost', db='eblm',
                         user='jmcc', password='mysqlpassword') as cur:
        cur.execute(qry)
        results = cur.fetchall()
        if len(results) > 0:
            bjd = np.empty(len(results))
            rv = np.empty(len(results))
            for i, result in enumerate(results):
                bjd[i] = float(result[0])
                if instrument == 'ids':
                    rv[i] = float(result[1]) + float(result[2])
                else:
                    rv[i] = float(result[1])
            # now phase the data
            cur.execute(qry2)
            results2 = cur.fetchone()
            epoch = float(results2[0])
            period = float(results2[1])
            phase = (bjd - (epoch + 2450000))/period%1
            return phase, rv
        else:
            return None, None

def calculateSemiAmplitude(phase, rv):
    """
    Take the phased data and work out K if possible
    """
    q1, q2 = [], []
    for i, p in enumerate(phase):
        if p >= 0.15 and p <= 0.35:
            q1.append(rv[i])
        if p >= 0.65 and p <= 0.85:
            q2.append(rv[i])
    if len(q1) > 0 and len(q2) > 0:
        return round(abs(np.average(q2) - np.average(q1))/2.0, 3)
    else:
        return None

def getSemiAmplitudes(swasp_ids):
    """
    Work out K for any spectra possible
    """
    Ks = OrderedDict()
    for swasp_id in swasp_ids:
        Ks[swasp_id] = {}
        for instrument in table_map:
            phase, rv = getPhasedRvData(swasp_id, instrument)
            if phase is not None and rv is not None:
                K = calculateSemiAmplitude(phase, rv)
            else:
                K = None
            Ks[swasp_id][instrument] = K
    return Ks

def updateDatabaseRow(swasp_id, column, value):
    """
    Add a row to the database. If it exists, prompt to clobber
    """
    # check it is there first
    qry = """
        SELECT {}
        FROM eblm_parameters
        WHERE swasp_id = '{}'
        """.format(column, swasp_id)
    with pymysql.connect(host='localhost', db='eblm',
                         user='jmcc', password='mysqlpassword') as cur:
        cur.execute(qry)
        result = cur.fetchone()
        if result[0] is None:
            go = True
        else:
            yn = raw_input('Clobber value {}? (y/n)'.format(result[0]))
            if yn.lower() == 'y':
                go = True
            else:
                go = False
        if go:
            qry2 = """
                UPDATE eblm_parameters
                SET {}={}
                WHERE swasp_id='{}'
                LIMIT 1
                """.format(column, value, swasp_id)
            cur.execute(qry2)

if __name__ == "__main__":
    args = argParse()
    if args.swasp_id:
        swasp_ids = [args.swasp_id]
    else:
        swasp_ids = getListOfSb1s()
    fwhms = getFwhms(swasp_ids)
    ks = getSemiAmplitudes(swasp_ids)

    # work out which to ingest - fies > cafe > ids
    ingest_fwhms = {}
    for swasp_id in fwhms:
        fwhm_to_log = None
        if fwhms[swasp_id]['ids'] is not None:
            fwhm_to_log = fwhms[swasp_id]['ids']
        if fwhms[swasp_id]['cafe'] is not None:
            fwhm_to_log = fwhms[swasp_id]['cafe']
        if fwhms[swasp_id]['fies'] is not None:
            fwhm_to_log = fwhms[swasp_id]['fies']
        ingest_fwhms[swasp_id] = fwhm_to_log

    ingest_ks = {}
    for swasp_id in ks:
        k_to_log = None
        if ks[swasp_id]['ids'] is not None:
            k_to_log = ks[swasp_id]['ids']
        if ks[swasp_id]['cafe'] is not None:
            k_to_log = ks[swasp_id]['cafe']
        if ks[swasp_id]['fies'] is not None:
            k_to_log = ks[swasp_id]['fies']
        ingest_ks[swasp_id] = k_to_log

    # now do the ingesting
    for swasp_id in ingest_fwhms:
        if ingest_fwhms[swasp_id]:
            print('Updating {} ccf_fwhm_avg --> {}'.format(swasp_id, ingest_fwhms[swasp_id]))
            updateDatabaseRow(swasp_id, 'ccf_fwhm_avg', ingest_fwhms[swasp_id])
        else:
            print('Skipping {} ccf_fwhm_avg update'.format(swasp_id))

    for swasp_id in ingest_ks:
        if ingest_ks[swasp_id]:
            print('Updating {} semi_amplitude --> {}'.format(swasp_id, ingest_ks[swasp_id]))
            updateDatabaseRow(swasp_id, 'semi_amplitude', ingest_ks[swasp_id])
        else:
            print('Skipping {} semi_amplitude update'.format(swasp_id))
