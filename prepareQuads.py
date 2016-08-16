"""
Prepare quads for CAFE/IDS obs

To do:
"""
import math
import sys
import argparse as ap
from collections import defaultdict
import pymysql
import ephem
from astropy.time import Time
import astropy.units as u
import numpy as np
from numpy.polynomial import Polynomial as P

# connect to database
db = pymysql.connect(host='localhost', db='eblm')

def argParse():
    description = "Code to calculate QUADS for CAFE/IDS objects"
    status_choices = ["CONTINUE-CAFE",
                      "CONTINUE-STABILIZED",
                      "CONTINUE",
                      "OBSERVE"]
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
    parser.add_argument('--snr', help='Target SNR for obs', type=int, default=30)
    return parser.parse_args()

# observatory set up
telescope = ephem.Observer()
if args.instrument == 'CAFE':
    telescope.lat = str(37. + (13./60.) + (25./3600.))
    telescope.lon = str(-2. - (32./60.) - (46./3600.))
    telescope.elev = 2168
else:
    telescope.lon=str(-17.-(52./60.))
    telescope.lat=str(28.+(40./60.))
    telescope.elev=2326
# some observing constraints
MOON_ANGLE_LIMIT= 30
ELEVATION_LIMIT = 30

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
    predicted = math.ceil(predicted/60)*60
    print('Rounding up to nearest 60s {0:d}'.format(predicted))

    # check for texp>2700, scale to right number of
    # spectra to get required SNR
    if predicted > 2700:
        print('Predicted time > 2700s, looking for good combo of spectra')

        snr_max = snr[V[n]][-1]
        print('SNR_max @ 2700 = {0:.2f}'.format(snr_max))

        n_spectra = (args.snr/snr[V[n]][-1])**2
        print('This needs {0:.2f} spectra @ 2700'. format(n_spectra))

        n_spectra = math.ceil(n_spectra)
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

    print("%s %.2f %d x %.2fs" % (swasp_id,vmag,n_spectra,predicted))
    return n_spectra, predicted

def getIdsExptime():
    pass

def getObjects(status):
    objects = {}
    qry = """
        SELECT
        swasp_id, period, epoch, Vmag, paramfit_spec_type, q1, q2
        FROM
        eblm_parameters
        WHERE
        current_status = '{0:s}'
        """.format(status)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            objects[row[0]] = (row[1], row[2]+2450000, row[3], row[4], row[5], row[6])
    return objects

if __name__ == '__main__':
    args = argParse()
    if len(args.night1.split('-')) != 3 or len(args.night2.split('-')) != 3:
        print('Dates must have YYYY-MM-DD format')
        sys.exit(1)

    JD1 = (Time(args.night1, format='isot', scale='utc', in_subfmt='date') + 0.5*u.day).jd
    JD2 = (Time(args.night2, format='isot', scale='utc', in_subfmt='date') + 1.5*u.day).jd

    # dictionary to hold the plan
    obs = {}

    # list to hold quads
    Q1_sched, Q2_sched = [],[]

    # get the objects from the database
    objects = getObjects(args.status_flag)
    for obj in objects:
        Epoch = objects[obj][1]
        Period = objects[obj][0]
        Vmag = objects[obj][2]

        # get expected exposure time for args.snr
        if args.instrument == 'CAFE':
            n_spectra, exptime = getCafeExptime(obj, Vmag)
        # fluff IDS for now
        else:
            n_spectra = 1
            exptime = 900

        # generate an ephem star
        ra = "{0:s}:{1:s}:{2:s}".format(obj[7:9],obj[9:11],obj[11:16])
        dec = "{0:s}:{1:s}:{2:s}".format(obj[16:19],obj[19:21],obj[21:25])
        star=ephem.FixedBody()
        star._ra=ephem.hours(ra)
        star._dec=ephem.degrees(dec)
        print(obj, ra, dec,  objects[obj][0], objects[obj][1])

        # generate list of quads
        # can do this better. cycle the 
        # epoch forward or backwards until the E before the run
        # then simply scale it through the run from there
        # no need for a billion epochs for nothing... fix later
        nepoch = 20000
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
        Q1r=np.linspace(En+(0.25*Period),En+(0.25*Period)+(n_periods*Period),n_periods+1)
        Q2r=np.linspace(En+(0.75*Period),En+(0.75*Period)+(n_periods*Period),n_periods+1)

        # convert the times to something usable mathematically
        Ert=Time(Er,format='jd',scale='utc')
        Q1rt=Time(Q1r,format='jd',scale='utc')
        Q2rt=Time(Q2r,format='jd',scale='utc')

        # ok so now we have a list of Q1 and Q2
        # cycle over them and check them for:
        # Moon, Twilight and Altitude

        # loop over Q1s
        if obj[5] != 1:
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
                            obs[Q1rt[k].iso] = "Q1 {0:s} alt={1:.2f} V={2:.2f} Texp={3:d}x{4:d} MoonSep={5:03d} MoonIll={6:d}%".format(obj, alt_q1, Vmag, exptime, n_spectra, int(ms_q1), int(phase_q1))
                            # append this object to a Q1 list
                            Q1_sched.append(obj)
        else:
            print('Q1 done for {}, skipping'.format(obj[0]))

        # loop over the Q2s
        if obj[6] != 1:
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
                            obs[Q2rt[k].iso] = "Q2 {0:s} alt={1:.2f} V={2:.2f} Texp={3:d}x{4:d} MoonSep={5:03d} MoonIll={6:d}%".format(obj, alt_q2, Vmag, exptime, n_spectra, int(ms_q2), int(phase_q2))
                            Q2_sched.append(obj)

    # print the final plan
    for i in sorted(obs):
        print(i, obs[i])

