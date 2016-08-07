import glob as g
import os

t=g.glob('*.fits')
for i in t:
    nn = 'r%s' % (i.split('_00')[1])
    print i, nn
    os.system('mv %s %s' % (i, nn))
