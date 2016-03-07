# function to cycle through the EBLMs and 
# reorder their aperture numbers according to the PA
# assume the first PA is the reference angle
# update the n_traces keyword to be -TRACE_ID for each trace, 
# after applying the flip
#
from collections import defaultdict
import numpy as np
import pymysql

db=pymysql.connect(host='localhost',db='eblm')

def updateTraceOrders(image_ids,trace_orders):
	for i in range(0,len(image_ids)):
		qry="UPDATE eblm_ids SET n_traces=%d WHERE image_id='%s' LIMIT 1" % (trace_orders[i],image_ids[i])
		print qry
		with db.cursor() as cur:
			cur.execute(qry)
			db.commit()

image_ids=defaultdict(list)
n_spectra=defaultdict(list)
n_traces=defaultdict(list)
pa=defaultdict(list)
qry="SELECT image_id,object_name,sky_pa,n_traces FROM eblm_ids WHERE n_traces<0 AND current_status != 'IGNORE' ORDER BY object_name"
with db.cursor() as cur:
	cur.execute(qry)
	for row in cur:
		image_ids[row[1]].append(row[0])
		pa[row[1]].append(float(row[2]))
		n_traces[row[1]]=int(row[3])
for i in image_ids.keys():
	n=[]
	for j in image_ids[i]:
		n.append(int(j.split('-')[-1]))
	n_traces[i]=max(n)
	n_spectra[i]=len(image_ids[i])/n_traces[i]


# now loop over the image_ids and check the SkyPAs
for i in pa.keys():
	if len(pa[i]) == n_traces[i]:
		print "\nThere is only one image for %s" % (i)
		updateTraceOrders(image_ids[i],range(-1,-n_traces[i]-1,-1))
		continue
	else:
		print "\nMultiple images for %s - checking PAs" % (i)
		print image_ids[i]
	
		print "Setting order based on 1st image"
		print image_ids[i][:n_traces[i]]
		for k in range(0,n_traces[i]):
			print i,image_ids[i][k],pa[i][k]
		updateTraceOrders(image_ids[i][:n_traces[i]],range(-1,-n_traces[i]-1,-1))

		# loop over the rest and flip any that might be needed
		for j in range(n_traces[i],len(pa[i]),n_traces[i]):
			invert=False
			if abs(pa[i][j]-pa[i][0]) >10 and abs(pa[i][j]-pa[i][0]) < 200:
				print "Need to invert traces of image %s for object %s" % (image_ids[i][j],i) 
				invert=True
			else:
				print "Traces aligned already"
			for k in range(j,j+n_traces[i]):
				print i,image_ids[i][k],pa[i][k]
			
			if invert==True:
				updateTraceOrders(image_ids[i][j:j+n_traces[i]],range(-1,-n_traces[i]-1,-1)[::-1])
			else:
				updateTraceOrders(image_ids[i][j:j+n_traces[i]],range(-1,-n_traces[i]-1,-1))
	print "\n"


