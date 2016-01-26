# script to measure the RVs for a given object

import pymysql, sys
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as pl
import argparse as ap

def argParse():
	parser=ap.ArgumentParser()
	parser.add_argument('--track',type=int,help='starting id')
	return parser.parse_args()

args=argParse()

# function to generate a list of unique objects
def getAllInfo(object_in):
	image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
	qry="SELECT image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment FROM eblm_ids WHERE object_name='%s'" % (object_in)
	print qry
	cur.execute(qry)	
	for row in cur:
		image_id.append(row[0])             
		object_name.append(row[1])          
		template_image_id.append(row[2])    
		ref_star_name.append(row[3])        
		hjd_mid.append(row[4])              
		utmiddle.append(row[5])             
		n_traces.append(row[6])             
		peak_shift_pix.append(row[7])       
		correlation_height.append(row[8])   
		fwhm_peak_pix.append(row[9])        
		fwhm_peak_kms.append(row[10])        
		relative_velocity_kms.append(row[11])
		observed_velocity_kms.append(row[12])
		helio_velocity_kms.append(row[13])   
		comment.append(row[14]) 
	
	# move numbers to np.array as they are nicer to work on later
	hjd_mid=np.array(hjd_mid)
	n_traces=np.array(n_traces)
	sky_pa=np.array(sky_pa)
	peak_shift_pix=np.array(peak_shift_pix)
	correlation_height=np.array(correlation_height)
	fwhm_peak_pix=np.array(fwhm_peak_pix)
	fwhm_peak_kms=np.array(fwhm_peak_kms)
	relative_velocity_kms=np.array(relative_velocity_kms)
	observed_velocity_kms=np.array(observed_velocity_kms)
	helio_velocity_kms=np.array(helio_velocity_kms)
	
	qry="SELECT eblm_parameters.swasp_id,period,epoch,vmag,paramfit_spec_type,eblm_ids.object_name FROM eblm_parameters INNER JOIN eblm_ids ON eblm_parameters.swasp_id=eblm_ids.swasp_id WHERE object_name='%s'" % (object_in)	
	print qry
	cur.execute(qry)	
	for row in cur:
		swasp_id=row[0] 
		period=row[1]
		epoch=row[2]
		vmag=row[3]
		spectype=row[4]
	
	try:
		swasp_id
		return image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment,swasp_id,period,epoch,vmag,spectype
	except NameError:
		print "\n\nNO SWASPID FOUND FOR %s\n\n" % (object_in)
		return image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment,"",1.0,5000.0000,15.0,"XX"
	
def getUniqueObjectList():
	qry="SELECT object_name,n_traces FROM eblm_ids"
	cur.execute(qry)
	object_list=defaultdict(list)
	
	for row in cur:
		object_list[row[0]].append(row[1])
	return object_list

def phasePlot(swasp_id,object_name,epoch,period,hjd_mid,relative_velocity_kms,fwhm_peak_kms):
	t=((hjd_mid-epoch)/period)%1
	pl.figure(1)
	pl.errorbar(t,relative_velocity_kms,yerr=fwhm_peak_kms,fmt='ro')
	pl.xlim(0.0,1.0)
	pl.title('%s - %s' % (swasp_id,object_name))
	pl.xlabel('Phase')
	pl.ylabel('RV (km/s)')
	pl.show()


# set up the database
db=pymysql.connect(host='localhost',db='eblm')
cur=db.cursor()
	
blends = 'n'
object_list=getUniqueObjectList()
track=1
for i in object_list.keys():
	sum=0
	for j in object_list[i]:
		sum+=j
	if sum==len(object_list[i]):
		# if selecting the track ID to startr on, skip to it
		if args.track and track < args.track:
			track+=1
			continue
		# if this is true there are no blends
		# get all info about object X	
		image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment,swasp_id,period,epoch,vmag,spectype=getAllInfo(i)
		
		# check for unmatched objects and skip them
		if swasp_id=="":
			continue
		
		p2p=max(relative_velocity_kms)-min(relative_velocity_kms)
		
		# print a summary of the objects data
		print "\nTracker: %d" % (track) 
		print "SWASPID: %s" % (swasp_id)
		print "Target:  %s, n_RV: %d" % (i, len(image_id))
		print "Period:  %.6f d, Epoch: %.6f, V: %.2f SpecType: %s" % (period,epoch,vmag,spectype)
		print "Peak to Peak: %.4f km/s, Error: %.4f km/s " % (p2p,np.average(fwhm_peak_kms))
		print "\t HJD\t\tVhelio\t\tObserved\tRV\t\tRV_e\t\tComments"
		for k in range(0,len(image_id)):
			print "\t %.6f\t%010.4f\t%010.4f\t%010.4f\t%010.4f\t%s" % (hjd_mid[k],helio_velocity_kms[k],observed_velocity_kms[k],relative_velocity_kms[k],fwhm_peak_kms[k],comment[k])
		
		plt=raw_input('Phase & Plot?: ')
		if plt.lower()=='y':
			phasePlot(swasp_id,i,epoch,period,hjd_mid,relative_velocity_kms,fwhm_peak_kms)
		elif plt.lower()=='x':
			sys.exit(1)
		elif plt.lower().startswith('np'):
			while plt.startswith('np'):
				phasePlot(swasp_id,i,epoch,float(plt.split()[1]),hjd_mid,relative_velocity_kms,fwhm_peak_kms)
				plt=raw_input('Phase & Plot?: ')	

		track+=1

			