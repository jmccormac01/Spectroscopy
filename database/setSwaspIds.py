"""
Script to set the swasp_id for reduced IDS spectra
"""
import argparse as ap
from datetime import timedelta
from collections import defaultdict
import pymysql

db = pymysql.connect(host='localhost',
                     db='eblm',
                     password='mysqlpassword')

def argParse():
    """
    Parse the command line arguments
    """
    parser = ap.ArgumentParser()
    parser.add_argument('--update',
                        help='update the database table?',
                        action='store_true')
    return parser.parse_args()

def getObjectsImageIds():
    """
    Get a list of images that have no swasp_id
    """
    objects = defaultdict(list)
    qry = """
        SELECT image_id, object_name, night
        FROM eblm_ids_newest
        WHERE swasp_id IS NULL
        """
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            objects[row[1]].append([row[0], row[2]])
    return objects

def getMatchingSwaspIds(object_name, objects):
    """
    Get a list of matching swasp_ids from 
    the database
    """
    swasp_ids = []
    qry = """
        SELECT swasp_id
        FROM eblm_parameters
        WHERE swasp_id LIKE '%{}%'
        """.format(object_name)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            swasp_ids.append(row[0])
    # make the checks for the number of matches
    if len(swasp_ids) == 0:
        print('No swasp_id found for {}'.format(object_name))
        for image_id in objects[object_name]:
            print(image_id)
        print('Find the object and add it to the database then rerun the script')
    elif len(swasp_ids) > 1:
        print('Multiple objects found for {}'.format(object_name))
        print('Pick one from those below:')
        for i, swasp_id in enumerate(swasp_ids):
            print('\t{}\t{}'.format(i, swasp_id))
        for image_id in objects[object_name]:
            print(image_id)
        choice = 1E6
        while choice not in range(len(swasp_ids)):
            choice = int(input('Choice: '))
        print('{} matched to {}'.format(object_name, swasp_ids[choice]))
        return swasp_ids[choice]
    else:
        return swasp_ids[0]

def updateDatabaseSwaspIds(object_name, swasp_id):
    """
    Take the paired object_names and swasp_ids and
    update the database now
    """
    qry = """
        UPDATE eblm_ids_newest
        SET swasp_id = '{}'
        WHERE object_name = '{}'
        """.format(swasp_id, object_name)
    with db.cursor() as cur:
        cur.execute(qry)
        db.commit()

if __name__ == "__main__":
    args = argParse()
    objects = getObjectsImageIds()
    swasp_ids = {}
    for object_name in sorted(objects):
        if not object_name.startswith('HD'):
            swasp_id = getMatchingSwaspIds(object_name, objects)
            swasp_ids[object_name] = swasp_id
        else:
            print('Skipping {}...'.format(object_name))
    # now update the database setting the correct swasp_id
    for object_name in sorted(swasp_ids):
        print(object_name, swasp_ids[object_name])
        if args.update:
            updateDatabaseSwaspIds(object_name, swasp_ids[object_name])

