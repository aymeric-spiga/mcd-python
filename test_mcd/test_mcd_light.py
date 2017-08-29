#! /usr/bin/env python

from fmcd import call_mcd
import numpy as np

lon = float(raw_input('Longitude?'))
lat = float(raw_input('Latitude?'))
xdate = float(raw_input('Ls?'))
loct = float(raw_input('Local time?'))
xz = float(raw_input('Altitude?'))

#dset = '/home/aymeric/Science/MCD_v4.3/data/'
dset = './MCD_DATA/'

zkey      = 3  # specify that xz is the altitude above surface (m)
datekey   = 1  # 1 = "Mars date": xdate is the value of Ls
dust      = 2  #our best guess MY24 scenario, with solar average conditions
hrkey     = 1  #set high resolution mode on (hrkey=0 to set high resolution off)
perturkey = 0  #integer perturkey ! perturbation type (0: none)
seedin    = 1  #random number generator seed (unused if perturkey=0)
gwlength  = 0. #gravity Wave wavelength (unused if perturkey=0)
#extvarkey = 1
#extvarkeys = np.ones(100)
extvarkeys = np.zeros(100)

(pres, dens, temp, zonwind, merwind, \
 meanvar, extvar, seedout, ierr) \
 = \
call_mcd(zkey,xz,lon,lat,hrkey, \
 datekey,xdate,loct,dset,dust, \
 perturkey,seedin,gwlength,extvarkeys )

print "temperature is %.0f K, pressure is %.0f Pa, density is %5.3e kg/m3, zonal wind is %.1f m/s, meridional wind is %.1f m/s" % (temp,pres,dens,zonwind,merwind)

