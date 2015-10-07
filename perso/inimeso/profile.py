#! /usr/bin/env python

# get a multicolumn file with profiles of z,P,T,rho

from string import split 
import numpy as np 
import matplotlib.pyplot as mpl
from mcd import mcd

### MCD INSTANCE and SETTINGS (actually, default. but one never knows)
query = mcd() ; query.zkey = 3 ; query.dust = 2 ; query.hrkey = 1
query.dust = 29
#query.dust = 24
query.toversion5(version="5.2")

### GET COORDINATES
lines = open("input_coord", "r").readlines()
query.lon   = float(split(lines[0])[0]) 
query.lat   = float(split(lines[1])[0]) 
query.xdate = float(split(lines[2])[0]) 
query.loct  = float(split(lines[3])[0])
query.printcoord()

### OPEN FILES TO BE WRITTEN
zefile = open("profile.txt", "w") 

### HEADER
query.gettitle()
zefile.write(query.ack + "\n")
zefile.write(query.title+"\n")
zefile.write("------------------------------\n")
zefile.write("altitude (m) // pressure (Pa) // density (kg/m3) // temperature (K)\n")
zefile.write("------------------------------\n")

### GET AND WRITE PROFILE
query.profile( tabperso = np.append([0.1,5,10,20,50,100],np.linspace(200.,float(split(lines[4])[0])*1000.,float(split(lines[5])[0]))) )


for iz in range(len(query.prestab)):
    zefile.write( "%15.4e%15.4e%15.4e%15.4e\n" % ( \
                  query.xcoord[iz],query.prestab[iz],query.denstab[iz],query.temptab[iz]) )

### FINISH
zefile.close()
query.plot1d(["p","t","rho"]) ; mpl.show()

