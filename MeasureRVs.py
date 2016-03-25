# script to measure the RVs for a given object

import pymysql, sys
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as pl
import argparse as ap

# set up the database
db=pymysql.connect(host='localhost',db='eblm')

def argParse():
	parser=ap.ArgumentParser()
	parser.add_argument('blends',choices=['yblends','nblends'],help='analyse blended objects?')
	parser.add_argument('--track',type=int,help='starting id')
	parser.add_argument('--swasp_id',help="swasp_id of the object you want")
	return parser.parse_args()

args=argParse()

# function to generate a list of unique objects
def getAllInfo(object_in,trace_id=None):
	image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
	if trace_id:
		qry="SELECT image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment FROM eblm_ids_new WHERE object_name='%s' AND n_traces=%d AND current_status != 'IGNORE'" % (object_in,trace_id)
	else:
		qry="SELECT image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment FROM eblm_ids_new WHERE object_name='%s' AND current_status != 'IGNORE'" % (object_in)
	print qry
	with db.cursor() as cur:
		cur.execute(qry)	
		for row in cur:
			image_id.append(row[0])             
			object_name.append(row[1])          
			template_image_id.append(row[2])    
			ref_star_name.append(row[3])        
			hjd_mid.append(row[4])              
			utmiddle.append(row[5])             
			n_traces.append(row[6]) 
			sky_pa.append(row[7])            
			peak_shift_pix.append(row[8])       
			correlation_height.append(row[9])   
			fwhm_peak_pix.append(row[10])        
			fwhm_peak_kms.append(row[11])        
			relative_velocity_kms.append(row[12])
			observed_velocity_kms.append(row[13])
			helio_velocity_kms.append(row[14])   
			comment.append(row[15]) 
	
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
	
	qry="SELECT eblm_parameters.swasp_id,period,epoch,vmag,paramfit_spec_type,eblm_ids_new.object_name FROM eblm_parameters INNER JOIN eblm_ids_new ON eblm_parameters.swasp_id=eblm_ids_new.swasp_id WHERE object_name='%s'" % (object_in)	
	print qry
	with db.cursor() as cur:
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
	
def getUniqueObjectList(blends):
	if blends == 'nblends':
		qry="SELECT object_name,n_traces FROM eblm_ids_new WHERE n_traces = 1 AND current_status != 'IGNORE'"
	else:
		qry="SELECT object_name,n_traces FROM eblm_ids_new WHERE n_traces < 0 AND current_status != 'IGNORE'"
	with db.cursor() as cur:
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

def runAnalysis(track,i,image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment,swasp_id,period,epoch,vmag,spectype):
	p2p=max(relative_velocity_kms)-min(relative_velocity_kms)
	
	phase=((hjd_mid-epoch)/period)%1

	# print a summary of the objects data
	print "\nTracker: %d" % (track) 
	print "SWASPID: %s" % (swasp_id)
	print "Target:  %s, n_RV: %d" % (i, len(image_id))
	print "Period:  %.6f d, Epoch: %.6f, V: %.2f SpecType: %s" % (period,epoch,vmag,spectype)
	print "Peak to Peak: %.4f km/s, Error: %.4f km/s " % (p2p,np.average(fwhm_peak_kms))
	print "\t ImageID\t\t\tPA\tHJD\t\tPhase\tVhelio\t\tObserved\tRV\t\tRV_e\t\tCCF_Height\tCCF_FWHM\tComments"
	for k in range(0,len(image_id)):
		print "\t %s\t%.2f\t%.6f\t%.3f\t%010.4f\t%010.4f\t%010.4f\t%010.4f\t%.3f\t\t%.2f\t\t%s" % (image_id[k],sky_pa[k],hjd_mid[k],phase[k],helio_velocity_kms[k],observed_velocity_kms[k],relative_velocity_kms[k],fwhm_peak_kms[k],correlation_height[k],fwhm_peak_pix[k],comment[k])
	
	plt=raw_input('Phase & Plot?: ')
	if plt.lower()=='y':
		phasePlot(swasp_id,i,epoch,period,hjd_mid,relative_velocity_kms,fwhm_peak_kms)
	elif plt.lower()=='x':
		sys.exit(1)
	elif plt.lower().startswith('np'):
		while plt.startswith('np'):
			phasePlot(swasp_id,i,epoch,float(plt.split()[1]),hjd_mid,relative_velocity_kms,fwhm_peak_kms)
			plt=raw_input('Phase & Plot?: ')

	
object_list=getUniqueObjectList(args.blends)

track=1
for i in sorted(object_list.keys()):
	# if selecting the track ID to start on, skip to it
	if args.track and track < args.track:
		track+=1
		continue
	if args.blends == 'nblends':
		image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment,swasp_id,period,epoch,vmag,spectype=getAllInfo(i,1)
		# check for unmatched objects and skip them
		if swasp_id=="":
			continue
		if args.swasp_id and args.swasp_id != swasp_id:
			continue
		runAnalysis(track,i,image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment,swasp_id,period,epoch,vmag,spectype)
	else:
		nt=abs(min(object_list[i]))
		for j in range(1,nt+1):
			image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment,swasp_id,period,epoch,vmag,spectype=getAllInfo(i,-j)
			if swasp_id=="":
				continue
			runAnalysis(track,i,image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,sky_pa,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment,swasp_id,period,epoch,vmag,spectype)
	track+=1


				