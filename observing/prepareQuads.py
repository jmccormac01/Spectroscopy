"""
Prepare quads for CAFE/IDS obs

TODO:

"""
import math
import sys
import argparse as ap
from collections import defaultdict
import pymysql
import ephem
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from numpy.polynomial import Polynomial as P

# connect to database
db = pymysql.connect(host='localhost', db='eblm', password='mysqlpassword')

class Range(object):
    """
    Range class for use with argparse
    """
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
    def __repr__(self):
        return '{0}-{1}'.format(self.start, self.end)

def argParse():
    description = "Code to calculate QUADS for CAFE/IDS objects"
    status_choices = ["CONTINUE-CAFE",
                      "CONTINUE-STABILIZED",
                      "CONTINUE",
                      "OBSERVE",
                      "PHASE_COVERAGE",
                      "1SWASP"]
    instruments = ['CAFE', 'IDS']
    parser = ap.ArgumentParser(description=description)
    parser.add_argument('instrument',
                        help='Instrument to plan observations for',
                        choices=instruments)
    parser.add_argument('status_flag',
                        choices = status_choices,
                        help='Status flag to ID objects in DB')
    parser.add_argument('night1', help='YYYY-MM-DD of night1 of run')
    parser.add_argument('night2', help='YYYY-MM-DD of night2 of run')
    parser.add_argument('--phase1',
                        help='Orbital phase 1 to schedule',
                        type=float, default=0.25,
                        choices=[Range(0.0, 1.0)] )
    parser.add_argument('--phase2',
                        help='Orbital phase 2 to schedule',
                        type=float, default=0.75,
                        choices=[Range(0.0, 1.0)])
    parser.add_argument('--snr', help='Target SNR for obs', type=int, default=30)
    parser.add_argument('--giants', help='search giants also?', action='store_true')
    return parser.parse_args()

# parse the command line
args = argParse()
# observatory set up
telescope = ephem.Observer()
if args.instrument == 'CAFE':
    telescope.lat = str(37. + (13./60.) + (25./3600.))
    telescope.lon = str(-2. - (32./60.) - (46./3600.))
    telescope.elev = 2168
    telescope.horizon = '-18'
else:
    telescope.lon=str(-17.-(52./60.))
    telescope.lat=str(28.+(40./60.))
    telescope.elev=2326
    telescope.horizon = '-12'
telescope.pressure = 750
# some observing constraints
MOON_ANGLE_LIMIT= 30
ELEVATION_LIMIT = 40
DECLINATION_LIMIT = -20

def Deg(l):
    if float(l[0]) >= 0:
        return (float(l[0]) + float(l[1])/60. + float(l[2])/60./60.)
    else:
        return (float(l[0]) - (float(l[1])/60.) - (float(l[2])/3600.))

def getCafeExptime(swasp_id, vmag):
    V=np.arange(8,14.0,0.5)
    t=np.array([30,120,300,900,2700])

    snr=defaultdict(list)
    snr[8.]   = np.array([23.381,46.762,73.937,128.063,221.812]) 
    snr[8.5]  = np.array([18.572,37.144,58.731,101.724,176.192])
    snr[9.]   = np.array([15.306,30.612,48.402,83.835,145.206])
    snr[9.5]  = np.array([12.158,24.316,38.447,66.592,115.341])
    snr[10.]  = np.array([9.657,19.315,30.540,52.896,91.619])
    snr[10.5] = np.array([7.671,15.342,24.258,42.017,72.776])
    snr[11.]  = np.array([6.093,12.817,19.269,33.375,57.807])
    snr[11.5] = np.array([4.840,9.680,15.306,26.511,45.918])
    snr[12.]  = np.array([3.844,7.689,12.158,21.058,36.474])
    snr[12.5] = np.array([3.05,6.108,9.657,16.727,28.972])
    snr[13.]  = np.array([2.426,4.852,7.671,13.287,23.014])
    snr[13.5] = np.array([1.923,3.854,6.093,10.554,18.280])

    coeffs_store=defaultdict(list)
    for i in sorted(snr.keys()):
        coeffs=np.polyfit(t,snr[i],2)
        coeffs_store[i]=coeffs
        tn=np.arange(30,3060,60)
        besty=np.polyval(coeffs,tn)
    diff=V-vmag
    n=np.where(abs(diff)==min(abs(diff)))[0][0]
    p=P.fit(t,snr[V[n]],2)
    t1,t2=(p-args.snr).roots()
    if t1<t2:
        predicted=t1.real
    else:
        predicted=t2.real
    print('Predicted exptime {0:.2f}'.format(predicted))
    # round to the nearest 60
    predicted = int(math.ceil(predicted/60))*60
    print('Rounding up to nearest 60s {0:d}'.format(predicted))
    # check for texp>2700, scale to right number of
    # spectra to get required SNR
    if predicted > 2700:
        print('Predicted time > 2700s, looking for good combo of spectra')
        snr_max = snr[V[n]][-1]
        print('SNR_max @ 2700 = {0:.2f}'.format(snr_max))
        n_spectra = (args.snr/snr[V[n]][-1])**2
        print('This needs {0:.2f} spectra @ 2700'. format(n_spectra))
        n_spectra = int(math.ceil(n_spectra))
        print('Rounding up to next integer n_spectra = {0:d}'.format(n_spectra))
        # now work out the best exptime to use to combine 
        # n_spectra spectra to get the desired args.snr
        target_single_spectra_snr = args.snr/math.sqrt(n_spectra)
        print('Working out the target SNR per spectra...')
        print('To achive SNR_total = {0:d} we need {1:d} spectra of SNR_i = {2:.2f}'.format(args.snr, n_spectra, target_single_spectra_snr))
        snr_range = np.polyval(coeffs_store[V[n]], tn)
        diff_ss_snr = abs(snr_range - target_single_spectra_snr)
        loc_ss_snr = np.where(diff_ss_snr == min(diff_ss_snr))[0][0]
        predicted = tn[loc_ss_snr]
        print('The best exposure time for this combination of spectra is {0:d} x {1:d}s'.format(n_spectra, predicted))
    elif predicted < 0:
        predicted = 60
        n_spectra=1
    else:
        n_spectra=1
    print("{0:s} {1:.2f} {2:d} x {3:.2f}s".format(swasp_id,vmag,n_spectra,predicted))
    return n_spectra, predicted

def getIdsExptime(mag):
    """
    Below are exposure times to give SNR ~50-60
    with IDS assuming bright time, airmass 1.2,
    slit width 1.4 and seeing 1.2'' (H1800V +
    6500A)
    """
    mags = np.arange(8.0, 14.0, 0.5)
    exptimes = np.array([8, 12, 20, 30, 50, 80, 120,
                         180, 300, 500, 700, 1200])
    diff = abs(mags - mag)
    n = np.where(diff == min(diff))[0][-1]
    return exptimes[n]

def getObjects(status, giants):
    objects = {}
    qry = """
        SELECT
        swasp_id, period, epoch, Vmag, paramfit_spec_type, q1, q2
        FROM
        eblm_parameters
        WHERE
        current_status = '{0:s}'
        """.format(status)
    if not giants:
        qry = qry + " AND (giant_flag IS NULL OR giant_flag != 'giant')"
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            objects[row[0]] = (row[1], row[2]+2450000, row[3], row[4], row[5], row[6])
    return objects

def getSwaspObject(swasp_id):
    objects = {}
    qry = """
        SELECT
        swasp_id, period, epoch, Vmag, paramfit_spec_type, q1, q2
        FROM
        eblm_parameters
        WHERE
        swasp_id = '{0:s}'
        """.format(swasp_id)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            objects[row[0]] = (row[1], row[2]+2450000, row[3], row[4], row[5], row[6])
    return objects

def getObjectsForPhaseCoverage(giants):
    """
    For these objects we don't care about the time

    Just output them all with their exptimes and mags
    """
    qry = """
        SELECT swasp_id, Vmag, period, paramfit_spec_type
        FROM eblm_parameters
        WHERE current_status = 'PHASE_COVERAGE_IDS'
        """
    if not giants:
        qry = qry + " AND (giant_flag IS NULL OR giant_flag != 'giant')"
    catalogue = []
    with db.cursor() as cur:
        cur.execute(qry)
        print('swasp_id', 'V mag', 'period', 'spec_type', 'exptime','n_spectra')
        for row in cur:
            swasp_id = row[0]
            period = round(float(row[2]), 5)
            spec_type = row[3]
            mag = round(float(row[1]),2)
            exptime = getIdsExptime(mag)
            qry2 = """
                SELECT count(*)
                FROM eblm_ids_newest
                WHERE swasp_id='{}'
                """.format(swasp_id)
            with db.cursor() as cur2:
                cur2.execute(qry2)
                for row2 in cur2:
                    count = row2[0]
            print(swasp_id, mag, period, spec_type, exptime, count)
            catalogue.append(swasp_id)
    print('Catalogue:')
    for c in catalogue:
        printCatalogueLine(c)

def printCatalogueLine(line):
    """
    Print ING formatted catalogue line
    """

    name = '{}'.format(line.split('1SWASP')[1])
    ra1 = name[1:3]
    ra2 = name[3:5]
    ra3 = name[5:10]
    dec1 = name[10:13]
    dec2 = name[13:15]
    dec3 = name[15:]
    print("{} {} {} {} {} {} {} J2000".format(name, ra1, ra2, ra3, dec1, dec2, dec3))

if __name__ == '__main__':
    if len(args.night1.split('-')) != 3 or len(args.night2.split('-')) != 3:
        print('Dates must have YYYY-MM-DD format')
        sys.exit(1)
    if args.status_flag == 'PHASE_COVERAGE':
        getObjectsForPhaseCoverage(args.giants)
        sys.exit(0)
    JD1 = (Time(args.night1, format='isot', scale='utc', in_subfmt='date') + 0.5*u.day).jd
    JD2 = (Time(args.night2, format='isot', scale='utc', in_subfmt='date') + 1.5*u.day).jd
    # dictionary to hold the plan
    obs = {}
    # list to hold the catalogue for TCS
    catalogue = []
    # list to hold quads
    Q1_sched, Q2_sched = [],[]
    # get the object(s) from the database
    if args.status_flag == '1SWASP':
        in_swasp_id = input('Enter the swasp_id: ')
        objects = getSwaspObject(in_swasp_id)
    else:
        objects = getObjects(args.status_flag, args.giants)
    for obj in objects:
        Epoch = objects[obj][1]
        Period = objects[obj][0]
        Vmag = objects[obj][2]
        # get expected exposure time for args.snr
        if args.instrument == 'CAFE':
            n_spectra, exptime = getCafeExptime(obj, Vmag)
        # for IDS we just use the 50-60 SNR values from always
        else:
            n_spectra = 1
            exptime = getIdsExptime(Vmag)
        # generate an ephem star
        ra = "{0:s}:{1:s}:{2:s}".format(obj[7:9],obj[9:11],obj[11:16])
        dec = "{0:s}:{1:s}:{2:s}".format(obj[16:19],obj[19:21],obj[21:25])
        star=ephem.FixedBody()
        star._ra=ephem.hours(ra)
        star._dec=ephem.degrees(dec)
        print(obj, ra, dec,  objects[obj][0], objects[obj][1])

        # using astropy SkyCoord to cut on declination
        skycoord = SkyCoord(ra, dec, unit=(u.hourangle, u.degree), frame='icrs')
        if skycoord.dec.deg < DECLINATION_LIMIT:
            print('Target {} too low, skipping...'.format(obj))
            continue

        # generate list of quads
        nepoch = 30000
        ec = np.empty(nepoch)
        for j in range(0, nepoch):
            ec[j] = Epoch + j*Period
        diff = ec - JD1
        if min(abs(diff)) > Period:
            print('Add more epochs, quitting...')
            sys.exit()
        n=np.where(np.abs(diff) == np.min(np.abs(diff)))
        # count from here n_periods until after JD2
        n_periods=int((JD2-ec[n[0][0]])/Period)+1

        # make new lists of E (Er), Q1 (Q1r) and Q2 (Qr)
        En=ec[n[0][0]]
        Er=np.linspace(En,En+(n_periods*Period),n_periods+1)
        Q1r=np.linspace(En+(args.phase1*Period),En+(args.phase1*Period)+(n_periods*Period),n_periods+1)
        Q2r=np.linspace(En+(args.phase2*Period),En+(args.phase2*Period)+(n_periods*Period),n_periods+1)

        # convert the times to something usable mathematically
        Ert=Time(Er,format='jd',scale='utc')
        Q1rt=Time(Q1r,format='jd',scale='utc')
        Q2rt=Time(Q2r,format='jd',scale='utc')

        # ok so now we have a list of Q1 and Q2
        # cycle over them and check them for:
        # Moon, Twilight and Altitude

        # loop over Q1s
        if objects[obj][4] == 0 or args.phase1 != 0.25:
            for k in range(0, len(Q1rt)):
                if Q1rt[k].jd > JD1 and Q1rt[k].jd < JD2:
                    # moon + target location
                    m=ephem.Moon(Q1rt[k].iso)
                    telescope.date=Q1rt[k].iso
                    star.compute(telescope)
                    # separation
                    moon_sep=ephem.separation(m,star)
                    ms_q1=Deg(str(moon_sep).split(":"))
                    # phase
                    phase_q1=m.phase
                    # object altitude
                    alt_q1=Deg(str(star.alt).split(':'))
                    tc=Time(int(Q1rt[k].jd),format='jd',scale='utc')
                    telescope.date=tc.iso
                    e_twi1=telescope.next_setting(ephem.Sun(),use_center=True)
                    m_twi1=telescope.next_rising(ephem.Sun(),use_center=True)


                    # check moon and elevation limits
                    if ms_q1 >= MOON_ANGLE_LIMIT and alt_q1 > ELEVATION_LIMIT:
                        # if the pass check we are withing twilights
                        if Q1rt[k].jd > float(ephem.julian_date(e_twi1)) and Q1rt[k].jd < float(ephem.julian_date(m_twi1)):
                            obs[(Q1rt[k]-((exptime/2.)*u.second)).iso] = "Ph{} {} alt={} V={} Texp={}x{} MoonSep={} MoonIll={}%".format(round(args.phase1, 2), obj, int(alt_q1), round(Vmag, 1), exptime, n_spectra, int(ms_q1), int(phase_q1))
                            # append this object to a Q1 list
                            Q1_sched.append(obj)
                            # append the object to the catalogue file
                            if obj not in catalogue:
                                catalogue.append(obj)
        else:
            print('Ph{} done for {}, skipping'.format(round(args.phase1, 2), obj))

        # loop over the Q2s
        if objects[obj][5] == 0 or args.phase2 != 0.75:
            for k in range(0, len(Q2rt)):
                if Q2rt[k].jd > JD1 and Q2rt[k].jd < JD2:
                    # moon + target location
                    m=ephem.Moon(Q2rt[k].iso)
                    telescope.date=Q2rt[k].iso
                    star.compute(telescope)
                    # separation
                    moon_sep=ephem.separation(m,star)
                    ms_q2=Deg(str(moon_sep).split(":"))
                    # phase
                    phase_q2=m.phase
                    # object altitude
                    alt_q2=Deg(str(star.alt).split(':'))
                    tc=Time(int(Q2rt[k].jd),format='jd',scale='utc')
                    telescope.date=tc.iso
                    e_twi1=telescope.next_setting(ephem.Sun(),use_center=True)
                    m_twi1=telescope.next_rising(ephem.Sun(),use_center=True)

                    # check moon and elevation limits
                    if ms_q2 >= MOON_ANGLE_LIMIT and alt_q2 > ELEVATION_LIMIT:
                        # if the pass check we are withing twilights
                        if Q2rt[k].jd > float(ephem.julian_date(e_twi1)) and Q2rt[k].jd < float(ephem.julian_date(m_twi1)):
                            obs[(Q2rt[k]-((exptime/2.)*u.second)).iso] = "Ph{} {} alt={} V={} Texp={}x{} MoonSep={} MoonIll={}%".format(round(args.phase2, 2), obj, int(alt_q2), round(Vmag, 1), exptime, n_spectra, int(ms_q2), int(phase_q2))
                            Q2_sched.append(obj)
                            # append the object to the catalogue file
                            if obj not in catalogue:
                                catalogue.append(obj)
        else:
            print('Ph{} done for {}, skipping'.format(round(args.phase2, 2), obj))

    # print the final plan
    for i in sorted(obs):
        print(i, obs[i])

    print('Catalogue:')
    for c in catalogue:
        printCatalogueLine(c)
