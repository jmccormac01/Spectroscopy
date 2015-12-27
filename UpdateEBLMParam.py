# tool to update the EBLM parameters database table

import argparse as ap
import pymysql

# function to parse the command line
def argParse():
	parser=ap.ArgumentParser()
	parser.add_argument('swaspid',help='SWASPID to update')
	parser.add_argument('param',help='parameter to update')
	parser.add_argument('value',help='value of update')
	parser.add_argument('--debug',help='run in debugging mode',action='store_true')
	return parser.parse_args()
	
	
args=argParse()

# set up the database
db=pymysql.connect(host='localhost',db='eblm')
cur=db.cursor()

qry="UPDATE eblm_parameters SET %s='%s' WHERE swasp_id='%s'" % (args.param,args.value,args.swaspid)
cur.execute(qry)
db.commit()

cur.close()
db.close()
