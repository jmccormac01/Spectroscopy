import matplotlib.pyplot as pl
import argparse as ap
from math import sqrt, ceil
import numpy as np
import glob as g
import seaborn 
import sys

seaborn.axes_style("darkgrid")

period=16.9530955
epoch=3253.68317
rv_file='J234318_41_Combined.rv'
wasp_phot='1SWASPJ234318.41+295556.5_OF2328+3237_100_ORCA_TAMTFA.lc'
follow_up_phot=g.glob('*MCMCunbinned*')

plot_array_height=int(ceil(sqrt(len(follow_up_phot))))+1
fig,ax=pl.subplots(plot_array_height,2,figsize=(15,10))

if rv_file:
	t_rv,rv,err_rv=np.loadtxt(rv_file,usecols=[0,1,2],unpack=True)
	phase_rv=((t_rv-epoch)/period)%1
	ax[0,0].errorbar(phase_rv,rv,yerr=err_rv,fmt='ro')
	ax[0,0].errorbar(phase_rv-1,rv,yerr=err_rv,fmt='ro')
	ax[0,0].errorbar(phase_rv+1,rv,yerr=err_rv,fmt='ro')
	ax[0,0].set_ylabel('RV (km/s)')
	ax[0,0].set_xlim(-0.5,1.5)
if wasp_phot:
	t_phot,phot,err_phot=np.loadtxt(wasp_phot,usecols=[0,1,2],unpack=True)
	bin_fac=10
	t_phot_b=np.empty(len(t_phot)/bin_fac)
	phot_b=np.empty(len(t_phot)/bin_fac)
	q=0
	for j in range(0,len(t_phot)/bin_fac):
		t_phot_b[j]=((sum(t_phot[q:q+bin_fac]))/bin_fac)
		phot_b[j]=((sum(phot[q:q+bin_fac]))/bin_fac)
		q=q+bin_fac
	phase_phot_b=((t_phot_b-epoch)/period)%1

	ax[0,1].plot(phase_phot_b,phot_b,'r.')
	ax[0,1].plot(phase_phot_b-1,phot_b,'r.')
	ax[0,1].plot(phase_phot_b+1,phot_b,'r.')
	ax[0,1].set_ylabel('Diff mag')
	ax[0,1].set_ylabel('Phase')
	ax[0,1].set_ylim(0.05,-0.03)
	ax[0,1].set_xlim(-0.5,1.5)

if len(follow_up_phot) > 0:
	j,k=1,0
	for i in follow_up_phot:
		t_fphot,fphot,err_fphot=np.loadtxt(i,usecols=[0,1,2],unpack=True)
		phase_fphot=((t_fphot-epoch)/period)%1
		ax[j,k].plot(phase_fphot,fphot,'r.')
		ax[j,k].plot(phase_fphot-1,fphot,'r.')
		ax[j,k].plot(phase_fphot+1,fphot,'r.')
		ax[j,k].set_xlim(-0.005,0.015)
		ax[j,k].set_ylim(0.05,-0.02)

		# logic for cycling the plots
		k=k+1
		if k > 1: 
			k=0
			j=j+1

pl.show()

