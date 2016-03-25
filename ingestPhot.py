'''
Script to ingest phot for interesting objects

Process:
	Check there is no ingested phot already
	If so, stop (remove later manuall?)
	If not, ingest it
'''
import os
import sys
import pymysql
import glob as g
import numpy as np
import argparse as ap

def argParse():
	description = '''Script for ingesting photometry into database
					Select the instrument e.g. WASP | NITES
					Data should be in diff_mags for ingesting
					You can use the options 
					swasp_id = 'find' and 
					infile = 'find' 
					to ingest many light curves of the same instrument set up
					'''	
	filters=["U", "B", "V", "R", "I", "clear", "u", "g", "r", "i", "z"]
	parser=ap.ArgumentParser(description=description)
	parser.add_argument('swasp_id',help='swasp_id')
	parser.add_argument('infile',help='path/to/phot/file')
	parser.add_argument('instrument',choices=['wasp','nites'],help='instrument')
	parser.add_argument('filter',choices=filters,help='filter used for observations')
	return parser.parse_args()

def ingest(t,m,e,swasp_id,instrument,filt):
	db=pymysql.connect(host='localhost',db='eblm')
	with db.cursor() as cur:
		for i,j,k in zip(t,m,e):
			qry="INSERT INTO eblm_phot (swasp_id,hjd_mid,mag,mag_err,instrument,filter) VALUES ('%s',%f,%f,%f,'%s','%s')" % (swasp_id,i,j,k,instrument,filt)
			cur.execute(qry)
			db.commit()

args=argParse()

# ingest many of the same instrument setup
if args.infile == 'find' and args.swasp_id == 'find':
	lis=g.glob('1SWASP*')
	for i in lis:
		swasp_id=i.split('_')[0]
		print swasp_id
		t,m,e=np.loadtxt(i,usecols=[0,1,2],unpack=True)
		ingest(t,m,e,swasp_id,args.instrument,args.filter)
		print "Ingested %s" % (i)
else:
	try:
		t,m,e=np.loadtxt(args.infile,usecols=[0,1,2],unpack=True)
		ingest(t,m,e,args.swasp_id,args.instrument,args.filter)
		print "Ingested %s" % (args.infile)
	except:
		print('Cannot find %s' % (args.infile))
		sys.exit()



