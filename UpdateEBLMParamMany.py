import pymysql
db=pymysql.connect(host='localhost',db='eblm')

status='CONTINUE-CAFE'

swasp_ids=["1SWASPJ230829.73+172019.6",
		"1SWASPJ172302.43+200000.6",
		"1SWASPJ174448.25+230056.0",
		"1SWASPJ004936.10+354240.2",
		"1SWASPJ184200.52+435734.9",
		"1SWASPJ221004.68+201121.4",
		"1SWASPJ235048.80+440127.4"]

for i in swasp_ids:
	qry="UPDATE eblm_parameters SET current_status='%s' WHERE swasp_id='%s'" % (status,i)	
	print qry
	with db.cursor() as cur:
		cur.execute(qry)
		db.commit()
