#! /usr/bin/env python
from string import split 
import numpy as np 
from ppclass import pp
import ppplot
from planets import Venus

############################
zefile = "yorgl.nc"
ttt = 1
############################

### GET COORDINATES
lines = open("input_coord", "r").readlines()
lon   = float(split(lines[0])[0]) ; lat   = float(split(lines[1])[0])
xdate = float(split(lines[2])[0]) ; loct  = float(split(lines[3])[0])

### GCM
temp  = pp(file=zefile,t=ttt,x=lon,y=lat,var="temp" ).getf()
pres  = pp(file=zefile,t=ttt,x=lon,y=lat,var="pres" ).getf()
vitu  = pp(file=zefile,t=ttt,x=lon,y=lat,var="vitu" ).getf()
vitv  = pp(file=zefile,t=ttt,x=lon,y=lat,var="vitv" ).getf()
dtlwr = pp(file=zefile,t=ttt,x=lon,y=lat,var="dtlwr").getf()
dtswr = pp(file=zefile,t=ttt,x=lon,y=lat,var="dtswr").getf()
geop  = pp(file=zefile,t=ttt,x=lon,y=lat,var="geop").getf()

### COMPUTATIONS
refpres = pres[0] ; surftpot = temp[0] ; surfh = 0.
### -- tpot
tpot = temp*(refpres/pres)**(Venus.R()/Venus.cp)
### -- alt
#alt = Venus.H(T0=temp)*np.log(refpres/pres)
alt = geop / Venus.g
### -- rho
rho = pres / Venus.R() / temp

### CHECK WITH A PLOT
pl = ppplot.plot1d(f=alt,x=tpot).makeshow()
pl = ppplot.plot1d(f=alt,x=rho).makeshow()
pl = ppplot.plot1d(f=alt,x=temp).makeshow()

### OPEN FILES TO BE WRITTEN
sounding = open("input_sounding", "w")
additional = open("input_therm", "w")
more = open("input_more", "w")
hr = open("input_hr", "w")

### GET and WRITE SURFACE VALUES
sounding.write( "%10.2f%12.2f%12.2f\n" % (refpres/100.,surftpot,0.) )
more.write( "%10.2f%10.2f" % (surfh,surftpot) )

### GET and WRITE VERTICAL PROFILE
nz = len(pres)
for iz in range(nz):
  sounding.write( "%10.2f%12.2f%12.2f%12.2f%12.2f\n" % (alt[iz],tpot[iz],0.,vitu[iz],vitv[iz]) )
  additional.write( "%12.2f%12.2f%18.6e%18.6e%12.2f\n" % (Venus.R(),Venus.cp,pres[iz],rho[iz],temp[iz]) )
  hr.write( "%18.6e%18.6e\n" % (dtswr[iz],dtlwr[iz]) )

### CLOSE FILES
sounding.close()
additional.close()
more.close()
hr.close()
