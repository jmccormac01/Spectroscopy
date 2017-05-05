# Here we take updated tables of WASP EBLMs and 
# compare it to what is in the db already, then update the vlaues
# this should be done before each observation run to be sure
# that we have all the most recent eblms and their values are up to date
# take periods and epochs from transvis, mcmc then orion in order of 
# reliability
import sys
import pymysql

db = pymysql.connect(host='localhost', db='eblm')

eblm_file = '/Users/jmcc/Dropbox/Observing/INTObs/EBLM/201608_10/eblms.csv'
f = open(eblm_file, 'r').readlines()
# skip the title row
for i in range(1, len(f)):
    print('[{}:{}]'.format(i,len(f)))
    ingest = False
    vals = f[i].split()
    swasp_id = '{0:s}{1:s}'.format(vals[0], vals[1])

    # get what is in the eblm db already
    with db.cursor() as cur:
        qry = "SELECT epoch, period, current_status, last_update FROM eblm_parameters WHERE swasp_id='{0:s}'".format(swasp_id)
        qry_len = cur.execute(qry)
        if qry_len > 0:
            for row in cur:
                epoch = row[0]
                period = row[1]
                status = row[2]
                update = row[3]
            print('\n{} DB values {} {} {} {}'.format(swasp_id, epoch, period, status, update))
        else:
            print('\n{0:s} not in database, ingest!'.format(swasp_id))
            # ingest function here
            ingest = True

    # get the best period, updated, mcmc then orion
    if vals[6] != 'NULL' and vals[7] != 'NULL':
        b_epoch = vals[6]
        b_period = vals[7]
    elif vals[4] != 'NULL' and vals[5] != 'NULL':
        b_epoch = vals[4]
        b_period = vals[5]
    else:
        b_epoch = vals[2]
        b_period = vals[3]

    print('{} has updated ephem {} {}'.format(swasp_id, vals[6], vals[7]))
    print('{} has mcmc ephem {} {}'.format(swasp_id, vals[4], vals[5]))
    print('{} has orion ephem {} {}'.format(swasp_id, vals[2], vals[3]))
    print('{} BEST VALUES {} {}'.format(swasp_id, b_epoch, b_period))

    # ingest new objects here
    if ingest:
        qry = """
              INSERT INTO eblm_parameters
              (swasp_id, epoch, period, current_status, Vmag)
              VALUES
              ('{}', {}, {}, 'OBSERVE', {})
              """.format(swasp_id, b_epoch, b_period, vals[8])
        print(qry)
        with db.cursor() as cur:
            cur.execute(qry)
            db.commit()

    # or update or not current objects
    else:
        if str(b_epoch) != str(epoch) or str(b_period) != str(period):
            print('{} Different values Hunter/DB, check: '.format(swasp_id))
            choice = input('u = update, x=quit, k=keep: ')
            if choice == 'u':
                qry = """
                      UPDATE eblm_parameters
                      SET epoch = {}, period = {}
                      WHERE swasp_id = '{}'
                      """.format(b_epoch, b_period, swasp_id)
                print(qry)
                with db.cursor() as cur:
                    cur.execute(qry)
                    db.commit()
            elif choice == 'x':
                sys.exit(1)
            else:
                continue
        else:
            print('{} DB values match best WASP values'.format(swasp_id))
            continue

