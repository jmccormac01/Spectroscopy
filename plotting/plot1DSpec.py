from astropy.io import fits
import matplotlib.pyplot as plt
import glob as g
import numpy as np
from scipy.ndimage.interpolation import shift
from scipy.fftpack import fft, ifft
from scipy import conjugate, polyfit
from astropy.time import Time

t = g.glob('*.fits')
loc = np.empty(len(t))
halpha = 6562.8 # angstroms

dateobs, vhelio = [], []
# set up the figure
fig, ax = plt.subplots(1, figsize=(10, 5))
h = fits.open(t[0])
ref = h[0].data[40000:130000].astype(float)
crval = h[0].header['CRVAL1']
crpix = h[0].header['CRPIX1']
cdelt = h[0].header['CDELT1']
dateobs.append(Time(h[0].header['DATE-OBS'], format='isot', scale='utc'))
vhelio.append(float(h[0].header['VHELIO']))
wave = np.arange(ref.shape[0])*cdelt + crval
ax.plot(wave, ref, 'k-')
print(dateobs[0].jd, vhelio[0])

# calculate offset
for i in range(1, 2): #len(t)):
    h = fits.open(t[i])
    spec = h[0].data[40000:130000].astype(float)
    crval = h[0].header['CRVAL1']
    crpix = h[0].header['CRPIX1']
    cdelt = h[0].header['CDELT1']
    dateobs.append(Time(h[0].header['DATE-OBS'], format='isot', scale='utc'))
    vhelio.append(float(h[0].header['VHELIO']))
    print(dateobs[i].jd, vhelio[i])
    wave = np.arange(spec.shape[0])*cdelt + crval + 0.182756
    ax.plot(wave, spec, 'b-')

plt.show()
