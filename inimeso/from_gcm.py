#! /usr/bin/env python
from string import split 
import numpy as np 
from ppclass import pp
import ppplot


shiftday = 258 + 27
shiftday = 258

### GET COORDINATES
lines = open("input_coord", "r").readlines()
lon   = float(split(lines[0])[0]) ; lat   = float(split(lines[1])[0])
xdate = float(split(lines[2])[0]) ; loct  = float(split(lines[3])[0])
utc = loct - lon/15.

### GCM
req = pp()
req.file = "diagfi.nc"
req.var = "vmr_h2ovap" 
req.var = "temp"
req.x = lon
req.y = lat
req.t = shiftday + utc/24.
req.verbose = True
req.changetime = "mars_dayini"
#prof,xx,yy,zz,tt = req.getfd()

### MCD
z,tpot,q,u,v = np.loadtxt("input_sounding",skiprows=1,unpack=True)
r,cp,p,rho,t = np.loadtxt("input_therm",unpack=True)
hgt,tsurf = np.loadtxt("input_more",unpack=True)
q = 0.001*q

### PLOT

nnn=60
dd=10
ppplot.rainbow(cb="jet",num=1+(nnn/dd))

sdg = ppplot.plot1d()
sdg.f = z-hgt
sdg.x = q 
sdg.x = t
sdg.legend = 'mcd'
sdg.marker = 's'
sdg.color = 'w'
sdg.make()


nnn=60
dd=10
ppplot.rainbow(cb="jet",num=nnn/dd)

ppp = ppplot.plot1d()
ppp.marker = '.'
ppp.ymax = 9000.
ppp.xlabel = 'temperature (K)'
ppp.ylabel = 'altitude (m)'

for iii in range(0,nnn,dd):
   if iii == 0: 
     sdg.verbose = True
   req.t = 258 + iii + utc/24.
   prof,xx,yy,zz,tt = req.getfd()
   ppp.x = prof
   ppp.f = zz*1000.
   ppp.legend = 'gcm sol 258+'+str(iii)
   ppp.make()

ppplot.save(mode="png",filename="comp")



