# script to measure the RVs for a given object

import pymysql
from collections import defaultdict
import numpy as np

# set up the database
db=pymysql.connect(host='localhost',db='eblm')
cur=db.cursor()

# function to generate a list of unique objects
def getAllInfo(object):
	image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
	qry="SELECT image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment FROM eblm_ids WHERE object_name='%s'" % (object)
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
	return image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment
	
def getUniqueObjectList():
	qry="SELECT object_name,n_traces FROM eblm_ids"
	cur.execute(qry)
	object_list=defaultdict(list)
	
	for row in cur:
		object_list[row[0]].append(row[1])
	return object_list


#		swid=raw_input('SWASPID: ')
#		period=raw_input('Period: ')
#		epoch=raw_input('Epoch: ')
#		spectype=raw_input('Spec Type: ')

	
	
blends = 'n'
object_list=getUniqueObjectList()

for i in object_list.keys():
	sum=0
	for j in object_list[i]:
		sum+=j
	if sum==len(object_list[i]):
		# if this is true there are no blends
		# get all info about object X	
		image_id,object_name,template_image_id,ref_star_name,hjd_mid,utmiddle,n_traces,peak_shift_pix,correlation_height,fwhm_peak_pix,fwhm_peak_kms,relative_velocity_kms,observed_velocity_kms,helio_velocity_kms,comment=getAllInfo(i)
		
		# print a summary of the objects data
		print "Target: %s, n_RV: %d" % (i, len(image_id))
		print "\tHJD\t\tVhelio\t\tObserved\tRV\t\tRV_e\t\tComments"
		for k in range(0,len(image_id)):
			print "\t%.6f\t%010.4f\t%010.4f\t%010.4f\t%010.4f\t%s" % (hjd_mid[k],helio_velocity_kms[k],observed_velocity_kms[k],relative_velocity_kms[k],fwhm_peak_kms[k],comment[k])
		
		

		
		