import pymysql
db=pymysql.connect(host='localhost',db='eblm', password='mysqlpassword')

status='ANALYSE_BLENDS'

swasp_ids=[""]

for i in swasp_ids:
	qry="UPDATE eblm_parameters SET current_status='%s' WHERE swasp_id='%s'" % (status,i)	
	print qry
	with db.cursor() as cur:
		cur.execute(qry)
		db.commit()
