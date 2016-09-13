"""
Script to populate the RA/DEC of all the targets in
the database based on their swasp_id

Set both the hms/dms and degrees formats using SkyCoord
"""
import pymysql
from astropy.coordinates import SkyCoord
import astropy.units as u

db = pymysql.connect(host='localhost',
                     db='eblm',
                     password='mysqlpassword')

def fetchSwaspIds():
    """
    Grab all swasp_ids from database
    """
    swasp_ids = []
    qry="""
        SELECT swasp_id
        FROM eblm_parameters
        WHERE ra_hms IS NULL
        """
    print(qry)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            swasp_ids.append(row[0])
    return swasp_ids

def parseSwaspId(obj):
    """
    Parse the swasp_id and return as SkyCoord object
    """
    ra = "{0:s}:{1:s}:{2:s}".format(obj[7:9],obj[9:11],obj[11:16])
    dec = "{0:s}:{1:s}:{2:s}".format(obj[16:19],obj[19:21],obj[21:25])
    c = SkyCoord('{} {}'.format(ra, dec), unit=(u.hourangle, u.degree), frame='icrs')
    return c, ra, dec

def updateRaDec(obj, c, ra, dec):
    """
    Take the swasp_id and SkyCoord object
    and update the database table
    """
    qry = """
        UPDATE eblm_parameters
        SET ra_hms = '{}',
        dec_dms = '{}',
        ra_deg = {},
        dec_deg = {}
        WHERE swasp_id = '{}'
        """.format(ra,
                   dec,
                   round(c.ra.deg, 6),
                   round(c.dec.deg, 6),
                   obj)
    print(qry)
    with db.cursor() as cur:
        cur.execute(qry)
        db.commit()

if __name__ == '__main__':
    swasp_ids = fetchSwaspIds()
    for swasp_id in swasp_ids:
        c, ra, dec = parseSwaspId(swasp_id)
        updateRaDec(swasp_id, c, ra, dec)

