"""
This script does a fresh ingest of the Mamajek stellar parameters table
"""
import pymysql
infile = 'mamajek_table.txt'
f = open(infile).readlines()
table = [line for line in f if not line.startswith('#')]

header = table[0].rstrip().split()[:-1]
header = [i.replace('-', '_') for i in header]
body = [line.rstrip() for line in table[1:]]

for row in body:
    cols = row.split()[:-1]
    args = []
    for col in cols:
        if '...' in col or '....' in col:
            args.append(None)
        else:
            # remove some spurious characters from entries
            col = col.replace('*', '')
            col = col.replace(':', '')
            try:
                args.append(round(float(col), 5))
            except ValueError:
                args.append(col)

    qry = """
        REPLACE INTO dwarf_star_properties
        ({})
        VALUES 
        """.format(','.join(header))

    qry = qry + """
        (%s, %s, %s, %s, %s, %s, %s,
         %s, %s, %s, %s, %s, %s, %s,
         %s, %s, %s, %s, %s, %s, %s,
         %s, %s, %s, %s, %s, %s, %s)
        """
    with pymysql.connect(host='localhost', db='eblm', password='mysqlpassword') as cur:
        cur.execute(qry, tuple(args))


