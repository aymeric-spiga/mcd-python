#! /usr/bin/env python
import numpy as np

shift = -5.

# load data
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html
yy,yyy = np.loadtxt("input_more",unpack=True)
gg,ggg,gggg = np.loadtxt("input_sounding",unpack=True)
z,tpot,q,u,v = np.loadtxt("input_sounding",skiprows=1,unpack=True)
r,cp,p,rho,t = np.loadtxt("input_therm",unpack=True)


tpot_save = tpot
t = t + shift
tpot = t*(610./p)**(1.0/3.9)

print "------------------------------------"
print "shifts in potential temperature are:"
print tpot-tpot_save

more = open("input_more", "w")
more.write( "%10.2f%10.2f" % (yy,t[0]) ) 
more.close()
sounding = open("input_sounding", "w")
additional = open("input_therm", "w")
sounding.write( "%10.2f%12.2f%12.2f\n" % (gg[0],tpot[0],gggg[0]) )
for iz in range(len(p)):
    sounding.write( "%10.2f%12.2f%12.2f%12.2f%12.2f\n" % (z[iz],tpot[iz],q[iz],u[iz],v[iz]) )
    additional.write( "%12.2f%12.2f%18.6e%18.6e%12.2f\n" % ( r[iz],cp[iz],p[iz],rho[iz],t[iz] ) )
sounding.close()
additional.close()
