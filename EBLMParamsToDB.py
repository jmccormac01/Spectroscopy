# script to load combined mcmc and orion eblm table 
# into the database EBLM params table
# also to grab the spectral types from the Targets Description file
# and cross match the object_names with swasp_ids for easier analysis later

import pymysql
import argparse as ap
from collections import defaultdict

def argParse():
	parser=ap.ArgumentParser()
	parser.add_argument('--debug',help='run in debugging mode',action='store_true')
	parser.add_argument('--populate',help='populate the eblm parameters table',action='store_true')
	parser.add_argument('--swaspids',help='cross match the swasp_ids',action='store_true')
	return parser.parse_args()

args=argParse()

# set up the database
if not args.debug:
	db=pymysql.connect(host='localhost',db='eblm')
	cur=db.cursor()
	cur2=db.cursor()

# populate the eblm parameters table
if args.populate:
	base_dir="/Users/James/Documents/Observing/INTObs/EBLM/201508_09"
	orca_file="%s/EBLMS_orca_combined.txt" % (base_dir)
	spectype_file="%s/TargetsDescription.txt" % (base_dir) 
	
	# populate all the EBLMS
	f=open(orca_file,'r').readlines()
	for i in f:
		x=i.split()	
		qry="INSERT INTO eblm_parameters (swasp_id,period,epoch,Vmag) VALUES ('%s',%.6f,%.6f,%.2f)" % (x[0],float(x[1]),float(x[2]),float(x[3]))
		print qry
		if not args.debug:
			cur.execute(qry)
			db.commit()
		
	# then use the targets description list to give the spec types 
	# for the ones we looked at
	f2=open(spectype_file,'r.').readlines()
	for i in f2:
		x=i.split()
		qry="UPDATE eblm_parameters set paramfit_spec_type='%s' WHERE swasp_id='%s'" % (x[1],x[0])
		print qry
		if not args.debug:
			cur.execute(qry)
			db.commit()
	

# match object names to swasp_ids
if args.swaspids:
	# script to then check the data table for object names and match any swasp_ids
	# if there are more than one swasp_ids that match, skip and flag it for checking
	swasp_id=[]
	qry="SELECT swasp_id FROM eblm_parameters"
	if not args.debug:
		cur.execute(qry)
		for row in cur:
			swasp_id.append(row[0])
	
	qry="SELECT image_id,object_name from eblm_ids_new"
	multiples=defaultdict(list)
	if not args.debug:
		cur.execute(qry)
		
		for row in cur:
			matches=[]
			for j in swasp_id:
				if row[1] in j:
					matches.append(j)
			
			if len(matches) ==1:
				print "object_name: %s --> swasp_id: %s" % (row[1],matches[0])
				qry="UPDATE eblm_ids_new SET swasp_id='%s' WHERE image_id='%s'" % (matches[0],row[0])
				#print qry
				cur2.execute(qry)
				db.commit()
				continue
			elif len(matches) > 1:
				print "object_name: %s --> MULTIPLE clashes" % (row[1])
				multiples[row[1]]=matches
			else:
				print "object_name: %s NOT SWASP OBJECT?" % (row[1])
		
		print "DO THE MULTIPLES BY HAND!"
			
	
if not args.debug:
	cur.close()
	cur2.close()
	db.close()
