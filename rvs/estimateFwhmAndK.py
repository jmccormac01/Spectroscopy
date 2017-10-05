"""
Script to look through the SB1s and estimate K and FWHM_avg from available data

Run this script ocassionally to update the CCF average FWHM and K in the parameters table
"""
import argparse as ap
import pymysql

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
    fwhms = {}
    with pymysql.connect(host='localhost', db='eblm',
                         user='jmcc', password='mysqlpassword') as cur:
        for swasp_id in swasp_ids:
            fwhms[swasp_id] = {'fies': None, 'cafe': None, 'ids': None}
            for instrument in table_map:
                qry = """
                    SELECT mask_ccf_fwhm
                    FROM {}
                    WHERE current_status LIKE 'SB1%'
                    OR current_status = 'CONTINUE'
                    """.format(table_map[instrument])
                cur.execute(qry)
                results = cur.fetchall()
                if results is not None:
                    temp = np.array(len(results))
                    for i, row in enumerate(results):
                        temp[i] = row[0]
                    K = np.average(temp)
                    fwhms[swasp_id][instrument] = K
    return fwhms

if __name__ == "__main__":
    args = argParse()
    if arg.swasp_id:
        swasp_ids = [args.swasp_id]
    else:
        swasp_ids = getListOfSb1s()
    fwhms = getFwhms(swasp_ids)

