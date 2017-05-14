"""
Script to check for all spectra being in the all_spectra folder
"""
import os
import glob as g
from collections import defaultdict
import pymysql

if __name__ == "__main__":
    swasp_ids = defaultdict(list)
    object_names = defaultdict(list)
    failed = defaultdict(list)
    spec_dir = "/Users/jmcc/Dropbox/data/int/ids/eblm/all_spectra"
    qry = """
        SELECT image_id, swasp_id, object_name
        FROM eblm_ids_newest
        WHERE n_traces = 1
        """
    with pymysql.connect(host='localhost', db='eblm',
                         password='mysqlpassword') as cur:
        cur.execute(qry)
        results = cur.fetchall()
    for row in results:
        if row[1] is None:
            object_names[row[2]].append(row[0])
        else:
            swasp_ids[row[1]].append(row[0])
    for swasp_id in swasp_ids:
        if swasp_id != 'None':
            os.chdir('{}/{}'.format(spec_dir, swasp_id))
            templist = g.glob('*.fits')
            if len(templist) != len(swasp_ids[swasp_id]):
                failed[swasp_id].append('Mismatch!')
            for image in swasp_ids[swasp_id]:
                if image not in templist:
                    failed[swasp_id].append(image)
    for f in failed:
        print("{}:".format(f))
        print(failed[f])
