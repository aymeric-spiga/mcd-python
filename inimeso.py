#! /usr/bin/env python

### AS 10/05/2012. A python script to prepare initial state for idealized mesoscale runs.
###                use : ensure mcd class is working. fill in input_coord. execute inimeso.

from string import split ; import numpy as np ; import matplotlib.pyplot as mpl
from mcd import mcd

### MCD INSTANCE and SETTINGS (actually, default. but one never knows)
query = mcd() ; query.zkey = 3 ; query.dust = 2 ; query.hrkey = 1

### GET COORDINATES
lines = open("input_coord", "r").readlines()
query.lon   = float(split(lines[0])[0]) ; query.lat   = float(split(lines[1])[0]) 
query.xdate = float(split(lines[2])[0]) ; query.loct  = float(split(lines[3])[0])
query.printcoord()

### OPEN FILES TO BE WRITTEN
sounding = open("input_sounding", "w") ; additional = open("input_therm", "w") ; more = open("input_more", "w")

### GET and WRITE SURFACE VALUES
query.xz = 0. ; query.update() ; query.printmeanvar()
sounding.write( "%10.2f%12.2f%12.2f\n" % (query.pres/100.,query.temp*(610./query.pres)**(1.0/3.9),0.) )
more.write( "%10.2f%10.2f" % (query.extvar[1],query.extvar[14]) ) ; more.close()

### GET and WRITE VERTICAL PROFILE
query.profile( tabperso = np.append([0,1,5,10,20,50,100],np.linspace(200.,float(split(lines[4])[0])*1000.,float(split(lines[5])[0]))) )
for iz in range(len(query.prestab)):
    sounding.write(   "%10.2f%12.2f%12.2f%12.2f%12.2f\n" % ( \
                      query.extvartab[iz,2],query.temptab[iz]*(610./query.prestab[iz])**(1.0/3.9),\
                      0.,query.zonwindtab[iz],query.merwindtab[iz]) )
    additional.write( "%12.2f%12.2f%18.6e%18.6e%12.2f\n" % ( \
                      query.extvartab[iz,49],query.extvartab[iz,8],\
                      query.prestab[iz],query.denstab[iz],query.temptab[iz]) )

### FINISH
sounding.close() ; additional.close()
query.plot1d(["p","t","u","v"],vertplot=1) ; mpl.show()
