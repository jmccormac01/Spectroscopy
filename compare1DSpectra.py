from astropy.io import fits
import matplotlib.pyplot as pl
import glob as g

t=g.glob('*.fits')

fig,ax=pl.subplots(len(t),1,figsize=(15,10),sharex=True)
for i in range(0,len(t)):
	h=fits.open(t[i])[0].data
	ax[i].plot(h,'k-')
	ax[i].set_xlim(91000,110000)
	ax[i].set_ylim(0.000,0.0005)
pl.show()
