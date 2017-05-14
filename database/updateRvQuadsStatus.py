"""
Code to update the database with what quads have been observed
"""
import pymysql
import argparse as ap

# connect to the database
db = pymysql.connect(host='localhost', db='eblm', password='mysqlpassword')

def argParse():
    """
    Parse the command line arguments
    """
    parser = ap.ArgumentParser()
    parser.add_argument('--v',
                        help='increased verbosity',
                        action='store_true')
    parser.add_argument('--commit',
                        help='commit changes to the database',
                        action='store_true')
    return parser.parse_args()

def checkQuads(swasp_id, E, P):
    qry2 = """
        SELECT
        image_id, hjd_mid
        FROM eblm_ids_newest
        WHERE swasp_id='{}'
        AND analyse = 1
        """.format(swasp_id)
    if args.v:
        print(qry2)
    q1_list, q2_list = [], []
    with db.cursor() as cur2:
        cur2.execute(qry2)
        for row2 in cur2:
            image_id = row2[0]
            hjd_mid = float(row2[1])
            phase = ((hjd_mid-E)/P)%1
            if args.v:
                print(image_id, hjd_mid, phase)
            if phase >= 0.2 and phase <= 0.3:
                q1_list.append(image_id)
                if args.v:
                    print('Q1 HIT!')
            elif phase >= 0.7 and phase <= 0.8:
                q2_list.append(image_id)
                if args.v:
                    print('Q2 HIT! ')
            else:
                if args.v:
                    print('Q1/2 MISS...')
                pass
        # set the 0's again in case the period changes
        # and the quads are no longer hit
        # Q1
        if len(q1_list) > 0:
            updateQuad('q1', swasp_id, 1)
        else:
            updateQuad('q1', swasp_id, 0)
        # Q2
        if len(q2_list) > 0:
            updateQuad('q2', swasp_id, 1)
        else:
            updateQuad('q2', swasp_id, 0)

def updateQuad(quad, swasp_id, status):
    qry = """
        UPDATE eblm_parameters
        SET {} = {}
        WHERE swasp_id='{}'
        """.format(quad, status, swasp_id)
    if args.v:
        print(qry)
    if args.commit:
        with db.cursor() as cur:
            cur.execute(qry)
            db.commit()

if __name__ == "__main__":
    args = argParse()
    # loop over targets and set the quads toggles
    qry = """
          SELECT
          swasp_id, epoch, period, q1,q2
          FROM
          eblm_parameters
          """
    if args.v:
        print(qry)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            swasp_id = row[0]
            E = 2450000 + float(row[1])
            P = float(row[2])
            Q1 = int(row[3]) if row[3] else None
            Q2 = int(row[4]) if row[4] else None
            if args.v:
                print(swasp_id, E, P, Q1, Q2)
            # check the Quads
            checkQuads(swasp_id, E, P)

