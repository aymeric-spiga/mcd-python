#! /usr/bin/env python

### AS 10/05/2012. A python script to prepare initial state for idealized mesoscale runs.
###                use : ensure mcd class is working. fill in input_coord. execute inimeso.

from string import split ; import numpy as np ; import matplotlib.pyplot as mpl
from mcd import mcd

rho_dust = 2500.  # Mars dust density (kg.m-3)
grav = 3.72
ksi = 3. / 4. / rho_dust / grav
nueff = 0.5

### MCD INSTANCE and SETTINGS (actually, default. but one never knows)
query = mcd() ; query.zkey = 3 ; query.dust = 2 ; query.hrkey = 1
query.dust = 29
#query.dust = 24
query.toversion5(version="5.2")

### GET COORDINATES
lines = open("input_coord", "r").readlines()
query.lon   = float(split(lines[0])[0]) ; query.lat   = float(split(lines[1])[0]) 
query.xdate = float(split(lines[2])[0]) ; query.loct  = float(split(lines[3])[0])
query.printcoord()

### OPEN FILES TO BE WRITTEN
sounding = open("input_sounding", "w") ; additional = open("input_therm", "w") ; more = open("input_more", "w")
dust = open("input_dust", "w") ; water = open("input_water", "w")

### GET and WRITE SURFACE VALUES
query.xz = 0.1 ; query.update() ; query.printmeanvar()
wvapor = query.extvar[42] #1.5*1e-3
wice = query.extvar[44] #0.0
sounding.write( "%10.2f%12.2f%12.2f\n" % (query.pres/100.,query.temp*(610./query.pres)**(1.0/3.9),(wvapor+wice)*1e3) )
more.write( "%10.2f%10.2f" % (query.extvar[1],query.extvar[14]) ) ; more.close()

### GET and WRITE VERTICAL PROFILE
closesurf = [0.1,5,10,20,50,100]
#closesurf = [0.1,0.5,1,2,5,10,20,50,100]
query.profile( tabperso = np.append(closesurf,np.linspace(200.,float(split(lines[4])[0])*1000.,float(split(lines[5])[0]))) )
for iz in range(len(query.prestab)):

    wvapor = query.extvartab[iz,42] #1.5*1e-3 
    wice = query.extvartab[iz,44] #0.0

    sounding.write(   "%10.2f%12.2f%12.2f%12.2f%12.2f\n" % ( \
                      query.extvartab[iz,2],query.temptab[iz]*(610./query.prestab[iz])**(1.0/3.9),\
                      (wvapor+wice)*1e3,\
                      query.zonwindtab[iz],query.merwindtab[iz]) )
    additional.write( "%12.2f%12.2f%18.6e%18.6e%12.2f\n" % ( \
                      query.extvartab[iz,53],query.extvartab[iz,8],\
                      query.prestab[iz],query.denstab[iz],query.temptab[iz]) )
    water.write( "%18.6e%18.6e\n" % (wvapor*1e3,wice*1e3) )

    ### DUST PROFILES
    q = query.extvartab[iz,38] # extvar(38)= Dust mass mixing ratio (kg/kg)
    reff = query.extvartab[iz,39] # extvar(39)= Dust effective radius (m)
    print q,reff
    N = (grav*ksi*((1+nueff)**3)/np.pi)*q/(reff**3)
    dust.write( "%18.6e%18.6e\n" % (q,N) )

### FINISH
sounding.close() ; additional.close() ; dust.close()
query.plot1d(["p","t","u","v","h2ovap","h2oice"]) ; mpl.show()






##### ESTIMATE CCNs
#dustq,dustn = np.loadtxt("input_dust",unpack=True)
#############################################################################
#############################################################################
##fqdust = np.squeeze(query.extvartab[:,38])
#fwice = np.squeeze(query.extvartab[:,44])
#####
##### I: fwice, fqdust --> O: finter
#####
##import numpy as np
#from scipy.interpolate import InterpolatedUnivariateSpline
### remove profile points with a cloud
#x = np.arange(fwice.size)
#w = np.where(fwice < np.mean(fwice)/10.)
### spline-interpolate with remaining points
#yyq = dustq[w]
#xxq = x[w]
#splq = InterpolatedUnivariateSpline(xxq, yyq)
#yyn = dustn[w]
#xxn = x[w]
#spln = InterpolatedUnivariateSpline(xxn, yyn)
### get interpolated profile
### -- a good match for qdust + qccn
### -- treat negative values because of spline
#finterq = splq(x) - dustq
#finterq[np.where(finterq < 0.)] = -finterq[np.where(finterq < 0.)]
#fintern = spln(x) - dustn
#fintern[np.where(fintern < 0.)] = -fintern[np.where(fintern < 0.)]
### plot
#mpl.figure()
#mpl.plot(yyn,xxn,'bo')
#mpl.plot(dustn,x)
#mpl.plot(fintern,x)
#mpl.figure()
#mpl.plot(yyq,xxq,'bo')
#mpl.plot(dustq,x)
#mpl.plot(finterq,x)
#mpl.figure()
#mpl.plot(fwice,x,'b.')
#mpl.show()
#############################################################################
#############################################################################
#
#ccn = open("input_ccn", "w")
#for iz in range(len(query.prestab)):
#  print iz,dustq[iz],finterq[iz]
#  #ccn.write( "%18.6e%18.6e\n" % (finterq[iz],fintern[iz]) )
#  ccn.write( "%18.6e%18.6e\n" % (np.max(finterq),np.max(fintern)) )
#ccn.close()
