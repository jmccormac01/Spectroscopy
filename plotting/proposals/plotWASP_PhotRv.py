# script to plot WASP photometry and RVs from database
import matplotlib.pyplot as pl
from math import ceil
import argparse as ap
import numpy as np
import glob as g
import pymysql
import seaborn

def argParse():
	parser = ap.ArgumentParser()
	parser.add_argument('plot', choices=['save', 'show'], help = 'how to plot objects')
	parser.add_argument('plottype', choices=['full', 'zoom'], help = 'full plot or zoom on transit?')
	return parser.parse_args()

args = argParse()
seaborn.axes_style("darkgrid")
db = pymysql.connect(host='localhost',db='eblm')

# list of swasp objects to plot
swasp_ids = ["1SWASPJ230829.73+172019.6",
			"1SWASPJ174448.25+230056.0",
			"1SWASPJ004936.10+354240.2",
			"1SWASPJ184200.52+435734.9",
			"1SWASPJ221004.68+201121.4",
			"1SWASPJ235048.80+440127.4"]

height=int(ceil(len(swasp_ids)/2.))*2
fig,ax=pl.subplots(height,2,figsize=(10,12))

# plot counters
n,m=0,0
# different binning factors for each plot
# depends on the amount of raw data
bin_facs=np.array([10,10,15,15,20,10])
for swasp_id,bin_fac in zip(swasp_ids,bin_facs):
	# get the object parameters
	qry="SELECT period,epoch,Vmag,paramfit_spec_type FROM eblm_parameters WHERE swasp_id='%s'" % (swasp_id)
	with db.cursor() as cur:
		cur.execute(qry)
		for row in cur:
			period = row[0]
			epoch = row[1]
			vmag = row[2]
			spectype = row[3]

	# get the RV data
	hjd_rv, rv, rv_err = [], [], []
	qry2="SELECT hjd_mid,relative_velocity_kms,fwhm_peak_kms FROM eblm_ids_new WHERE swasp_id='%s' AND current_status != 'IGNORE' AND n_traces=1" % (swasp_id)
	with db.cursor() as cur:
		cur.execute(qry2)
		for row in cur:
			hjd_rv.append(row[0])
			rv.append(row[1])
			rv_err.append(row[2])
	hjd_rv = np.array(hjd_rv)
	rv = np.array(rv)
	rv_err = np.array(rv_err)
	
	hjd_phot, mag, mag_err = [], [], []
	qry3="SELECT hjd_mid,mag,mag_err FROM eblm_phot WHERE swasp_id='%s' AND instrument = 'wasp'" % (swasp_id)
	with db.cursor() as cur:
		cur.execute(qry3)
		for row in cur:
			hjd_phot.append(row[0])
			mag.append(row[1])
			mag_err.append(row[2])
		hjd_phot = np.array(hjd_phot)
		mag = np.array(mag)
		mag_err = np.array(mag_err)	

	# phase the RVs
	phase_rv=((hjd_rv-epoch)/period)%1
	
	hjd_phot_b=np.empty(len(hjd_phot)/bin_fac)
	mag_b=np.empty(len(hjd_phot)/bin_fac)
	q=0
	for j in range(0,len(hjd_phot)/bin_fac):
		hjd_phot_b[j]=((sum(hjd_phot[q:q+bin_fac]))/bin_fac)
		mag_b[j]=((sum(mag[q:q+bin_fac]))/bin_fac)
		q=q+bin_fac
	phase_phot_b=((hjd_phot_b-epoch)/period)%1

	spread=np.std(mag_b)
	print "%s\t%.2f\t%s\t%05d\t%d\t%.4f" % (swasp_id,vmag,spectype,len(hjd_phot),bin_fac,spread)


	# HERE NEED TO FIX STUPID PLOT!!

	#title="%s P=%.4f E=%.4f" % (swasp_id,period,epoch)
	#ax[n,m].set_title(title)
	ax[n,m].plot(phase_phot_b,mag_b,'k.',phase_phot_b-1,mag_b,'k.',phase_phot_b+1,mag_b,'k.')
	if args.plottype == 'full':
		ax[n,m].set_xlim(-0.5,1.5)
		ax[n,m].set_ylim(5*spread,-3*spread)
	else:
		ax[n,m].set_xlim(-0.25,0.25)
		ax[n,m].set_ylim(0.05,-0.05)

	
	ax[n+1,m].errorbar(phase_rv,rv,yerr=rv_err,fmt='k.')
	ax[n+1,m].set_xlim(0,1)

	# y labels
	if m == 0:
		ax[n,m].set_ylabel('Diff mag')
		ax[n+1,m].set_ylabel('RV (km/s)')	
	# incremenet the plotting area
	m+=1
	if m > 1:
		m=0
		n+=2

# x labels
ax[-1,0].set_xlabel('Orbital Phase')
ax[-1,1].set_xlabel('Orbital Phase')

if args.plot == 'show':
	pl.show()
elif args.plot == 'save':
	outname="EBLM_ProposalCAFE.png"
	fig.savefig(outname,dpi=200,bbox_inches='tight')


